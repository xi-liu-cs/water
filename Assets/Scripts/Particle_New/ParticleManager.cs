using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;
using Unity.Mathematics;
using Random = UnityEngine.Random;
using UnityEditor;

public class ParticleManager : MonoBehaviour
{

    public static ParticleManager current;

    public enum NeighborSearchType {
        PrefixSum,
        CubeVolume
    }

    public struct Particle {
        public float3 position;
        public float speed;
        public int isBoid;
        public int touchedByBoid;
        public float3 boidInfluence;
    }

    private bool initialized = false;

    [Header("== REFERENCES ==")]
        public Boids3D boidManager;
        public PlaneObstacleManager obstacleManager;

    [Header("== WORLD CONFIGURATIONS ==")]
        #if UNITY_EDITOR
        [Help("These configurations can be used to adjust the simulation's world parameters. This includes how big the grid cells are, where the center of the simulation is in world space, etc.\n\n- `Bounds` is where the simulation is limited to\n- `Outer Bounds` is the world space including the buffer cells.\n- `Num Cells Per Axis` s is automatically calculated based on `Grid Cell Size` and `Outer Bounds`.\n- `Num Grid Cells` is the total number of grid cells populating `Outer Bounds` and is auto-calculated.\n- `G` is the gravitational force exerted on all particles in the simulation.", UnityEditor.MessageType.None)]
        #endif
        
        [Tooltip("How big are our grid cells? For Prefix Summation, recommended to match `smoothingRadius`")] 
        public float gridCellSize = 2f;
        
        [Tooltip("Where in world space are we centering the simulation around?")]
        public float[] origin = {0f,0f,0f};
        public Vector3 originVector3 { get => new Vector3(origin[0], origin[1], origin[2]); set {} }
        
        [Tooltip("What are the total world space length (per axis) is the simulation?")]
        public float[] bounds = {50f, 50f, 50f};
        public Vector3 boundsVector3 { get => new Vector3(bounds[0], bounds[1], bounds[2]); set {} }
        [ReadOnly, SerializeField, Tooltip("In our grid space, how many buffer cells are added to each axis? Prefix Sum: 2; CubeVolume: 1")]
        private int bufferCellsPerAxis;
        [ReadOnly, SerializeField, Tooltip("What is the world space size of the simulation, after taking into account grid cell size and buffer cells?")]
        private float[] _outerBounds = {60f, 60f, 60f};
        public float[] outerBounds { get => _outerBounds; set {} }
        [ReadOnly, SerializeField, Tooltip("Given `_outerBounds`, how many grid cells are along each axis?")]
        private int[] _numCellsPerAxis;
        public Vector3Int numCellsPerAxis { get => new Vector3Int(_numCellsPerAxis[0], _numCellsPerAxis[1], _numCellsPerAxis[2]); set {} }
        public int[] numCellsPerAxisInt { get => _numCellsPerAxis; set {} }
        [ReadOnly, SerializeField, Tooltip("How many grid cells do we have in total?")]
        private int numGridCells;
        [Tooltip("What's the gravitational force exerted on all particles?")]
        public float[] g = {0f, -9.81f, 0f};

    [Header("== PARTICLE CONFIGURATIONS ==")]
        #if UNITY_EDITOR
        [Help("These configurations adjust the particles themselves, such as how many particles are present and how big they are in world scale.", UnityEditor.MessageType.None)]
        #endif

        [Tooltip("How many particles will we use in the simulation?")]
        public int numParticles = 100000;
        [Tooltip("How big (visually only!) are each particle?")]
        public float particleRenderSize = 0.5f;
        public float particleRenderRadius { get => particleRenderSize / 2f; set { particleRenderSize = value * 2f; } }
        [Tooltip("The mesh used to render each particle in the simulation. Usually just the default `Sphere` mesh from Unity.")]
        public Mesh particle_mesh;
        public Material particle_material;
        // How many particles can we realistically fit into each grid cell? Calculated from particle render size. Intentional to use radius instead of size
        private int _numParticlesPerGridCell => Mathf.CeilToInt(Mathf.Pow(gridCellSize / particleRenderRadius, 3));
        public int numParticlesPerGridCell { get => _numParticlesPerGridCell; set {} }

    [Header("== FLUID MECHANICS ==")]
        #if UNITY_EDITOR
        [Help("These configurations adjust the behavior of the SPH fluid dynamics. These can get complicated, so make sure you know SPH inside out first. Learn about them here:\nhttps://www.cs.cornell.edu/~bindel/class/cs5220-f11/code/sph.pdf\nhttps://matthias-research.github.io/pages/publications/sca03.pdf\nThe current implementation is based off of the 1st reference.\n\n- `Dt` = The time step between frames. If set to any value < 0, then defaults to `Time.deltaTime`\n- `Smoothing Radius` = `K` - the radius used for the smoothing kernel.\n- `Particle Mass` = the mass of each individual particle.\n- `Viscosity Coefficient` = `mu`, how much particles stick to each other.\n- `Rest Density` = `p_0`, the default density the fluid is meant to embody when still.\n- `Damping` = How much the boundary walls reflect particles upon contact.\n - `Gas_constant` = how much the particles are excited. Normally dependent on temperature. Not used in this implementation.\n- `Bulk Modulus` = Every fluid has a modulus value - this can be looked up online. Try to pick values that aren't too large or small.", UnityEditor.MessageType.None)]
        #endif

        [Tooltip("The time difference between frames. If set to a negative number, will default to `Time.deltaTime`.")]
        public float dt = 0.0165f;
        [Tooltip("`h`: the smoothing kernel radius for SPH. Recommended to set to the same as `gridCellSize`.")]
        public float smoothingRadius = 1f;
        public float radius2 { get => Mathf.Pow(smoothingRadius, 2f); set {} }
        public float radius3 { get => Mathf.Pow(smoothingRadius, 3f); set {} }
        public float radius4 { get => Mathf.Pow(smoothingRadius, 4f); set {} }
        public float radius5 { get => Mathf.Pow(smoothingRadius, 5f); set {} }
        public float radius6 { get => Mathf.Pow(smoothingRadius, 6f); set {} }
        public float radius8 { get => Mathf.Pow(smoothingRadius, 8f); set {} }
        public float radius9 { get => Mathf.Pow(smoothingRadius, 9f); set {} }
        public float radius16 { get => Mathf.Pow(smoothingRadius, 16f); set {} }
        [Tooltip("How much mass does each particle have?")]
        public float particleMass = 1f;
        [Tooltip("The viscosity coefficient for SPH. The higher the value, the more particles are likely to clump together.")]
        public float viscosity_coefficient = 12.5f;
        [Tooltip("The resting density of particles.")]
        public float rest_density = 1f;
        [Tooltip("The amount of influence that boundaries have on particles when particles collide with boundaries. Must be a negative number b/w 0 and -1. Recommended: -0.5")]
        public float damping = -0.3f;
        [Tooltip("The gas constant of the particle liquid. The higher this value, the more particles vibrate and launch themselves in the air. Depends on temperature in the real world.")]
        public float gas_constant = 0f;
        [Tooltip("The bulk modulus of the fluid. Can be looked up on Google. Water = 300,000 psi")]
        public float bulkModulus = 1000f;

    [Header("== BOID CONFIGURATIONS ==")]
        #if UNITY_EDITOR
        [Help("The first `numBoids` particles are considered boids and will be constantly cemented to the boids from a separate boids manager. We just need to let the SPH compute shader know which of its particles are boids", UnityEditor.MessageType.None)]
        #endif
        [ReadOnly] public int numBoids = 100;
        public float boidMass = 5f;
    
    [Header("== OBSTACLE CONFIGURATIONS ==")]
        public PointCloudObstacle[] obstacles = new PointCloudObstacle[0];
        [SerializeField]
        private int numPointsOnObstacles = 0;

    [Header("== GPU SETTINGS ==")]
        #if UNITY_EDITOR
        [Help("GPU settings! Right now, the neighbor search type is set to `Cube Volume`, which limits the number of neighbors tracked per grid cell based on how many particles can realistically fit into a grid cell's volume based on `Grid Cell Size` and `Particle Render Size` defined above. `Prefix Sum` methodology does not work with the current implementation.", UnityEditor.MessageType.None)]
        #endif

        [Tooltip("Reference to the compute shader used to control particle movement")]
        public ComputeShader SPHComputeShader;
        [Tooltip("What kind of neighbor search and parsing should we use?")]
        public NeighborSearchType neighborSearchType = NeighborSearchType.CubeVolume;
        // number of threads running in each thread group that runs the simulation
        private const int _BLOCK_SIZE = 1024;
        public int BLOCK_SIZE { get => _BLOCK_SIZE; set {} }
        public int NUM_THREADS_FOR_PARTICLES { get => Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE); set {} }
        // number of thread groups that run on the GPU
        private int numBlocks;

    [Header("== GIZMOS CONFIGURATIONS ==")]
        [Tooltip("Toggle gizmos. Simple enough.")]
        public bool show_gizmos = false;
        [Tooltip("Show particles via Gizmos")]
        public bool show_particles = false;
        [Tooltip("Show boids")]
        public bool show_boids = false;
        [Tooltip("Show boid influence")]
        public bool show_boid_influence = false;
        [Tooltip("Show grid cells that have occupying particles via Gizmos")]
        public bool show_grid_cells = false;
        [Tooltip("Show the velocity vectors of each particle")]
        public bool show_velocities = false;
        [Tooltip("Show the force vectors of each particle")]
        public bool show_forces = false;
        [Tooltip("What color should the particles be via Gizmos?")]
        public Color gizmos_particle_color = Color.blue;

    [Header("== DEBUG CONFIGURATIONS ==")]
        public bool verbose_grid = false;
        public bool verbose_offsets = false;
        public bool verbose_prefix_sums = false;
        public bool verbose_rearrange = false;
        public bool verbose_num_neighbors = false;
        public bool verbose_density_pressure = false;
        public bool verbose_force = false;
        public bool verbose_velocity = false;
    
    // == KERNELS FOR GPU METHODS ==
    // Initialization and Resetting (Universal)
    int generate_particles_kernel;
    int clear_grid_kernel;
    int reset_num_neighbors_kernel;
    int update_grid_kernel;
    // Specific to prefix summation method
    int prefix_sum_kernel;
    int sum_blocks_kernel;
    int add_sums_kernel;
    int rearrange_particles_kernel;
    int count_num_neighbors_kernel;
    int ps_compute_density_kernel;
    int ps_compute_interact_acceleration_kernel;
    // Specific to cube volume method
    int cv_compute_density_kernel;
    int cv_compute_interact_acceleration_kernel;
    // Update particle positions (Universal)
    int compute_external_acceleration_kernel;
    int integrate_kernel;
    int dampen_by_bounds_kernel;
    int integrate_boid_particles_kernel;

    // public Vector3 velocity_initial = new Vector3(0, 10, 0);

    [Header("voxel")]
    public int n_point_per_axis = 50; 
    // z * n_point_per_axis * n_point_per_axis + y * n_point_per_axis + x, x ^ 3 + x ^ 2 + x = numParticles, (x = 50, numParticles = 130000)
    public float isolevel = 8;

    int size_property = Shader.PropertyToID("size");
    int particle_buffer_property = Shader.PropertyToID("particle_buffer");

    [Header("== COMPUTE BUFFERS ==")]
    public ComputeBuffer arg_buffer;
    
    // Stores how many particles are inside each grid cell. 
        // Length = numGridCells
        // To get grid cell count: 
        //  1) Find projected grid cell hash index
        //  2) grid[projected_index] <- the # of particles in that grid cell.
    private ComputeBuffer gridBuffer;
    // Stores all particle data in the system. 
        // Length = numParticles
    public ComputeBuffer particle_buffer;
    // Stores the "offsets" of each particle.
        // Length = numParticles
        // For example, each grid cell has `n` number of particles; `offset` is the index of the particle in their grid cell (max: `n-1`)
    public ComputeBuffer particleOffsetsBuffer;
    // Stores ID of particles in each grid cell. Length = numGridCells * _numParticlesPerGridCell
    // To iterate through neighbors of a particle's current cell:
    //  1) Get neighbor cells' projected indices. Can be done by getting current cell's XYZ indices, then iterating through 27 neighbor cells
    //  2) Iterate through neighbor cells. For each neighbor cell:
    //      2a) Get their starting index `j` for particleNeighborsBuffer (neighbor's hashed index * _numParticlesPerGridCell)
    //      2b) Get # of neighbors `n` in that neighbor cell (gridBuffer[<neighbor's hashed index>])
    //      2c) Loop `i` through `j` to `j+(n-1)`. Neighbor ID = particleNeighborsBuffer[i]
    public ComputeBuffer particleNeighborsBuffer;

    // The stuff below are specific to prfix summation
    private ComputeBuffer gridOffsetBuffer;
    private ComputeBuffer gridSumsBuffer1;
    private ComputeBuffer gridSumsBuffer2;
    private ComputeBuffer rearrangedParticlesBuffer;
    private ComputeBuffer numNeighborsBuffer;
    private ComputeBuffer obstacleParticlesBuffer;

    // Stores the masses of each particle
    public ComputeBuffer mass_buffer;
    // Stores the densities of each particle
    public ComputeBuffer density_buffer;
    // Stores the pressure exerted by each particle
    public ComputeBuffer pressure_buffer;
    // Stores the velocity of each particle
    public ComputeBuffer velocity_buffer;
    // Stores the force (acceleration) of each particle
    public ComputeBuffer force_buffer;
    // A temporary compute buffer for ayn thing debugging
    public ComputeBuffer temp_buffer;

    // These ones are ... older. Might use them later, but for now they're here
    public ComputeBuffer neighbor_tracker_buffer;
    public ComputeBuffer hash_grid_buffer;
    public ComputeBuffer hash_grid_tracker_buffer;
    public ComputeBuffer point_buffer;
    public ComputeBuffer noise_density_buffer;
    public ComputeBuffer triangle_buffer;
    public ComputeBuffer triangle_count_buffer;

    void OnDrawGizmos() {
        if (!show_gizmos) return;

        Gizmos.color = Color.white;
        Gizmos.DrawWireCube(originVector3, new Vector3(bounds[0], bounds[1], bounds[2]));
        Gizmos.color = Color.yellow;
        Gizmos.DrawWireCube(originVector3, new Vector3(_outerBounds[0], _outerBounds[1], _outerBounds[2]));
        Gizmos.color = Color.red;
        Gizmos.DrawSphere(originVector3,0.1f);

        if (!Application.isPlaying) return;

        Vector4 gridColor = new Vector4(1f,0f,0f,0.25f);
        Vector3 gridCellSize3D = new Vector3(gridCellSize, gridCellSize, gridCellSize);
        Particle[] temp_particles = new Particle[numParticles];
        particle_buffer.GetData(temp_particles);
        //int[] temp_num_neighbors = new int[numParticles];
        //numNeighborsBuffer.GetData(temp_num_neighbors);
        //float3[] temp_forces_array = new float3[numParticles];
        //force_buffer.GetData(temp_forces_array);
        float3[] temp_velocities_array = new float3[numParticles];
        velocity_buffer.GetData(temp_velocities_array);
        for(int i = 1; i < numParticles; i++) {
            if (show_particles) {
                Gizmos.color = gizmos_particle_color;
                Gizmos.DrawSphere(temp_particles[i].position, particleRenderRadius);
            }
            if (show_boid_influence) {
                if (temp_particles[i].isBoid > 0) {
                    Gizmos.color = Color.red;
                    Gizmos.DrawRay(temp_particles[i].position, temp_particles[i].boidInfluence);
                    Gizmos.color = Color.blue;
                    Gizmos.DrawRay(temp_particles[i].position, temp_velocities_array[i]);
                }
            }
            if (show_grid_cells) {
                Gizmos.color = gridColor;
                Gizmos.DrawCube(GetGridCellWorldPositionFromGivenPosition(temp_particles[i].position), gridCellSize3D);
            }
            /*
            if (show_forces) {
                Gizmos.color = Color.green;
                Gizmos.DrawLine(temp_particles[i].position, temp_particles[i].position + temp_forces_array[i]);
            }
            if (show_velocities) {
                Gizmos.color = Color.blue;
                Gizmos.DrawLine(temp_particles[i].position, temp_particles[i].position + temp_velocities_array[i]);
            }
            */
        }

        if (show_boids) {
            Gizmos.color = Color.yellow;
            for(int i = 0; i < numBoids; i++) {
                Gizmos.DrawSphere(temp_particles[i].position, particleRenderRadius);
            }
        } else {
            Gizmos.color = Color.yellow;
            Gizmos.DrawSphere(temp_particles[0].position, particleRenderRadius);
            Gizmos.color = gridColor;
            Gizmos.DrawCube(GetGridCellWorldPositionFromGivenPosition(temp_particles[0].position), gridCellSize3D);
            Gizmos.color = Color.yellow;
            Gizmos.DrawWireSphere(temp_particles[0].position, smoothingRadius);
        }
        /*
        Vector3 handlePos = new Vector3(
            temp_particles[0].position[0],
            temp_particles[0].position[1] + particleRenderSize*2f,
            temp_particles[0].position[2] + particleRenderSize*2f
        );
        Handles.Label(handlePos, temp_num_neighbors[0].ToString());
        */
    }

    public Vector3Int GetGridXYZIndices(Vector3 position) {
        return new Vector3Int(
            Mathf.FloorToInt((position.x - (origin[0] - (_numCellsPerAxis[0] * gridCellSize)/2f))/gridCellSize),
            Mathf.FloorToInt((position.y - (origin[1] - (_numCellsPerAxis[1] * gridCellSize)/2f))/gridCellSize),
            Mathf.FloorToInt((position.z - (origin[2] - (_numCellsPerAxis[2] * gridCellSize)/2f))/gridCellSize)
        );
    }
    public int GetProjectedGridIndexFromXYZ(int x, int y, int z) {
        return x + (_numCellsPerAxis[0] * y) + (_numCellsPerAxis[0] * _numCellsPerAxis[1] * z);
    }
    public int GetProjectedGridIndexFromXYZ(Vector3Int xyz) {
        return xyz.x + (_numCellsPerAxis[0] * xyz.y) + (_numCellsPerAxis[0] * _numCellsPerAxis[1] * xyz.z);
    }
    public Vector3 GetGridCellWorldPositionFromXYZIndices(int3 xyz) {
        return new Vector3(
            (origin[0] - ((_numCellsPerAxis[0] * gridCellSize)/2f)) + (xyz[0] * gridCellSize) + (gridCellSize/2f),
            (origin[1] - ((_numCellsPerAxis[1] * gridCellSize)/2f)) + (xyz[1] * gridCellSize) + (gridCellSize/2f),
            (origin[2] - ((_numCellsPerAxis[2] * gridCellSize)/2f)) + (xyz[2] * gridCellSize) + (gridCellSize/2f)
        );
    }
    public Vector3 GetGridCellWorldPositionFromXYZIndices(Vector3Int xyz) {
        return new Vector3(
            (origin[0] - ((_numCellsPerAxis[0] * gridCellSize)/2f)) + (xyz.x * gridCellSize) + (gridCellSize/2f),
            (origin[1] - ((_numCellsPerAxis[1] * gridCellSize)/2f)) + (xyz.y * gridCellSize) + (gridCellSize/2f),
            (origin[2] - ((_numCellsPerAxis[2] * gridCellSize)/2f)) + (xyz.z * gridCellSize) + (gridCellSize/2f)
        );
    }
    public Vector3 GetGridCellWorldPositionFromGivenPosition(Vector3 position) {
        Vector3Int xyz = GetGridXYZIndices(position);
        return GetGridCellWorldPositionFromXYZIndices(xyz);
    }

    private void Awake() {
        current = this;
        // If this game object is its own active gameobject, we run initialize;
        if (!initialized) Initialize();
    }

    // This is a public method that can be called from some external manager or the like.
    // If this is a standalone component without any managers, the `Awake` method will handle this by itself. No management needed.
    public void Initialize() {
        // Initialize boids first
        boidManager.Initialize();
        // Initialize obstacles second
        obstacleManager.Initialize();
        // Initialize key variables once
        InitializeVariables();
        // Determine the kernels from the GPU
        InitializeKernels();
        InitializeShaderVariables();
        InitializeBuffers();
        boidManager.AddParticleBuffer();
        boidManager.dt = dt;
        SPHComputeShader.Dispatch(clear_grid_kernel, numBlocks,1,1);
        SPHComputeShader.Dispatch(generate_particles_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE),1,1);
        initialized = true;
        Debug.Log("Particles Initialized");
    }

    // We need to initialize some variables for the simulation. Only runs once.
    // The primary purpose of this method is to determine:
    // - how many cells we can fit per axis
    // - the new bounds after adjusting for intricacies with float values and buffer cell appending
    // - how many grid cells we have in total
    // - how many thread groups we need
    private void InitializeVariables() {
        // Determine our process type
        bufferCellsPerAxis = 4;
        // Calculate how many cells will fit within the provided dimensions
        _numCellsPerAxis = new int[3];
        _numCellsPerAxis[0] = Mathf.CeilToInt(bounds[0] / gridCellSize) + bufferCellsPerAxis;
        _numCellsPerAxis[1] = Mathf.CeilToInt(bounds[1] / gridCellSize) + bufferCellsPerAxis;
        _numCellsPerAxis[2] = Mathf.CeilToInt(bounds[2] / gridCellSize) + bufferCellsPerAxis;
        // Re-Calculate the new bounds based on the number of cells per axis
        _outerBounds[0] = (float)_numCellsPerAxis[0] * gridCellSize;
        _outerBounds[1] = (float)_numCellsPerAxis[1] * gridCellSize;
        _outerBounds[2] = (float)_numCellsPerAxis[2] * gridCellSize;
        // How many grid cells are there in total?
        numGridCells = _numCellsPerAxis[0] * _numCellsPerAxis[1] * _numCellsPerAxis[2];
        // For volumetric neighbor search, the # of particles per grid cell is calculated automatically. No need to calculate it here.

        // However, we need to adjust the number of boids to ensure that it isn't larger than the # of particles in the current system.
        numBoids = Mathf.Min(boidManager.numBoids, Mathf.FloorToInt((float)numParticles*0.9f));

        // We also need to count number of particles attached to obstacles
        numPointsOnObstacles = 0;
        foreach(PointCloudObstacle obstacle in obstacles) {
            numPointsOnObstacles += obstacle.points.Count;
        }

        // == GPU SETTINGS ==
        numBlocks = Mathf.CeilToInt((float)numGridCells / (float)_BLOCK_SIZE);

        Debug.Log($"Number of particle grid cells per axis: ({_numCellsPerAxis[0]}, {_numCellsPerAxis[1]}, {_numCellsPerAxis[2]})");
        Debug.Log($"Total # of particle grid cells: {numGridCells}");
        Debug.Log($"# Particles per grid cell: {_numParticlesPerGridCell}");
        Debug.Log($"Size of particle neighbors: {numGridCells * _numParticlesPerGridCell}");
    }

    // We need to find out which kernels we need to invoke in the GPU for each process.
    // Keep in mind that the kernels will differ based on which neighbor search type we end up using.
    void InitializeKernels() {
        // Initialization and Resetting (Universal). Reset num neighbors is purely for debugging purposes
        clear_grid_kernel = SPHComputeShader.FindKernel("ClearGrid");
        generate_particles_kernel = SPHComputeShader.FindKernel("GenerateParticles");
        reset_num_neighbors_kernel = SPHComputeShader.FindKernel("ResetNumNeighbors");
        update_grid_kernel = SPHComputeShader.FindKernel("UpdateGridCellCounts");
        // Specific to Prefix Summation Method
        prefix_sum_kernel = SPHComputeShader.FindKernel("PrefixSum");
        sum_blocks_kernel = SPHComputeShader.FindKernel("SumBlocks");
        add_sums_kernel = SPHComputeShader.FindKernel("AddSums");
        rearrange_particles_kernel = SPHComputeShader.FindKernel("RearrangeParticles");
        count_num_neighbors_kernel = SPHComputeShader.FindKernel("CountNumNeighbors");
        ps_compute_density_kernel = SPHComputeShader.FindKernel("PS_ComputeDensity");
        ps_compute_interact_acceleration_kernel = SPHComputeShader.FindKernel("PS_ComputeInteractAcceleration");
        // Specific to Cube Volume Method
        cv_compute_density_kernel = SPHComputeShader.FindKernel("CV_ComputeDensity");
        cv_compute_interact_acceleration_kernel = SPHComputeShader.FindKernel("CV_ComputeInteractAcceleration");
        // This one is universal across process types
        compute_external_acceleration_kernel = SPHComputeShader.FindKernel("ComputeExternalAcceleration");
        integrate_kernel = SPHComputeShader.FindKernel("Integrate");
        dampen_by_bounds_kernel = SPHComputeShader.FindKernel("DampenByBounds");
        integrate_boid_particles_kernel = SPHComputeShader.FindKernel("IntegrateBoidsParticles");
    }

    void InitializeShaderVariables() {
        // == WORLD CONFIGURATIONS ==
        SPHComputeShader.SetFloat("gridCellSize", gridCellSize);
        SPHComputeShader.SetInts("numCellsPerAxis", _numCellsPerAxis);
        SPHComputeShader.SetInt("total_number_of_cells", numGridCells);
        SPHComputeShader.SetFloats("origin", origin);
        SPHComputeShader.SetFloats("bounds", bounds);
        SPHComputeShader.SetFloats("outerBounds", _outerBounds);
        SPHComputeShader.SetFloats("g", g);
        SPHComputeShader.SetFloat("epsilon", Mathf.Epsilon);
        SPHComputeShader.SetFloat("pi", Mathf.PI);

        // == PARTICLE CONFIGURATIONS ==
        SPHComputeShader.SetInt("numParticles", numParticles);
        SPHComputeShader.SetInt("numBoids", numBoids);
        SPHComputeShader.SetInt("numParticlesOnObstacles", numPointsOnObstacles);
        SPHComputeShader.SetFloat("particleRenderRadius", particleRenderRadius);
        SPHComputeShader.SetInt("numParticlesPerGridCell", _numParticlesPerGridCell);

        // == GPU SETTINGS
        SPHComputeShader.SetInt("numBlocks", numBlocks);
        SPHComputeShader.SetInt("randomSeed", Random.Range(0,int.MaxValue));

        // == OLD SETTINGS ==
        //SPHComputeShader.SetInt("n", 8);
        //SPHComputeShader.SetFloat("e", Mathf.Exp(1));
        //SPHComputeShader.SetVector("time", Shader.GetGlobalVector("_Time"));
        SPHComputeShader.SetInt("n_point_per_axis", n_point_per_axis);
        SPHComputeShader.SetFloat("isolevel", isolevel);

        // Update variables that may change over time due to modifying inspector values
        UpdateShaderVariables();
    }

    public void UpdateShaderVariables(bool updateDT = false) {
        // == FLUID MECHANICS ==
        SPHComputeShader.SetFloat("smoothingRadius", smoothingRadius);
        SPHComputeShader.SetFloat("radius2", radius2);
        SPHComputeShader.SetFloat("radius4", radius4);
        SPHComputeShader.SetFloat("radius6", radius6);
        SPHComputeShader.SetFloat("radius8", radius8);
        SPHComputeShader.SetFloat("radius9", radius9);
        SPHComputeShader.SetFloat("radius16", radius16);

        SPHComputeShader.SetFloat("particleMass", particleMass);
        SPHComputeShader.SetFloat("gas_constant", gas_constant);
        SPHComputeShader.SetFloat("rest_density", rest_density);
        SPHComputeShader.SetFloat("viscosity_coefficient", viscosity_coefficient);
        SPHComputeShader.SetFloat("damping", damping);
        SPHComputeShader.SetFloat("bulkModulus", bulkModulus);

        SPHComputeShader.SetFloats("g", g);
        if (updateDT) {
            if (dt > 0) SPHComputeShader.SetFloat("dt", dt);
            else SPHComputeShader.SetFloat("dt", Time.deltaTime);
        }
    }

    public unsafe void InitializeBuffers() {
        uint[] arg = {particle_mesh.GetIndexCount(0), (uint)numParticles, particle_mesh.GetIndexStart(0), particle_mesh.GetBaseVertex(0), 0};
        arg_buffer = new ComputeBuffer(1, arg.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        arg_buffer.SetData(arg);

        particle_buffer = new ComputeBuffer(numParticles, sizeof(float)*7 + sizeof(int)*2);
        particleOffsetsBuffer = new ComputeBuffer(numParticles, sizeof(int));
        gridBuffer = new ComputeBuffer(numGridCells, sizeof(int));
        particleNeighborsBuffer = new ComputeBuffer(numGridCells * _numParticlesPerGridCell, sizeof(int));

        
        //obstacleParticlesBuffer = new ComputeBuffer(numPointsOnObstacles, sizeof(Vector3));


        float[] particle_masses = new float[numParticles];
        for(int i = 0; i < numParticles; i++) particle_masses[i] = boidMass;
        for(int j = numBlocks; j < numParticles; j++) particle_masses[j] = particleMass;
        mass_buffer = new ComputeBuffer(numParticles, sizeof(float));
        mass_buffer.SetData(particle_masses);

        density_buffer = new ComputeBuffer(numParticles, sizeof(float));
        pressure_buffer = new ComputeBuffer(numParticles, sizeof(float));
        velocity_buffer = new ComputeBuffer(numParticles, 3 * sizeof(float));
        force_buffer = new ComputeBuffer(numParticles, 3 * sizeof(float));

        gridOffsetBuffer = new ComputeBuffer(numGridCells, sizeof(int));
        gridSumsBuffer1 = new ComputeBuffer(numBlocks, sizeof(int));
        gridSumsBuffer2 = new ComputeBuffer(numBlocks, sizeof(int));
        rearrangedParticlesBuffer = new ComputeBuffer(numParticles, sizeof(int));
        numNeighborsBuffer = new ComputeBuffer(numParticles, sizeof(int));

        temp_buffer = new ComputeBuffer(numParticles, sizeof(int));
        
        SPHComputeShader.SetBuffer(clear_grid_kernel, "grid", gridBuffer);
        SPHComputeShader.SetBuffer(reset_num_neighbors_kernel, "numNeighbors", numNeighborsBuffer);

        SPHComputeShader.SetBuffer(generate_particles_kernel, "particles", particle_buffer);
        SPHComputeShader.SetBuffer(generate_particles_kernel, "density", density_buffer);
        SPHComputeShader.SetBuffer(generate_particles_kernel, "pressure", pressure_buffer);
        SPHComputeShader.SetBuffer(generate_particles_kernel, "force", force_buffer);
        SPHComputeShader.SetBuffer(generate_particles_kernel, "velocity", velocity_buffer);

        SPHComputeShader.SetBuffer(update_grid_kernel, "particles", particle_buffer);
        SPHComputeShader.SetBuffer(update_grid_kernel, "grid", gridBuffer);
        SPHComputeShader.SetBuffer(update_grid_kernel, "particleOffsets", particleOffsetsBuffer);
        SPHComputeShader.SetBuffer(update_grid_kernel, "particleNeighbors", particleNeighborsBuffer);
        SPHComputeShader.SetBuffer(update_grid_kernel, "temp_buffer", temp_buffer);

        SPHComputeShader.SetBuffer(prefix_sum_kernel, "gridOffsetBufferIn", gridBuffer);
        SPHComputeShader.SetBuffer(prefix_sum_kernel, "gridOffsetBuffer", gridOffsetBuffer);
        SPHComputeShader.SetBuffer(prefix_sum_kernel, "gridSumsBuffer", gridSumsBuffer2);

        SPHComputeShader.SetBuffer(add_sums_kernel, "gridOffsetBuffer", gridOffsetBuffer);

        SPHComputeShader.SetBuffer(rearrange_particles_kernel, "particles", particle_buffer);
        SPHComputeShader.SetBuffer(rearrange_particles_kernel, "rearrangedParticles", rearrangedParticlesBuffer);
        SPHComputeShader.SetBuffer(rearrange_particles_kernel, "gridOffsetBuffer", gridOffsetBuffer);
        SPHComputeShader.SetBuffer(rearrange_particles_kernel, "particleOffsets", particleOffsetsBuffer);

        SPHComputeShader.SetBuffer(count_num_neighbors_kernel, "particles", particle_buffer);
        SPHComputeShader.SetBuffer(count_num_neighbors_kernel, "rearrangedParticles", rearrangedParticlesBuffer);
        SPHComputeShader.SetBuffer(count_num_neighbors_kernel, "gridOffsetBuffer", gridOffsetBuffer);
        SPHComputeShader.SetBuffer(count_num_neighbors_kernel, "numNeighbors", numNeighborsBuffer);

        SPHComputeShader.SetBuffer(ps_compute_density_kernel, "particles", particle_buffer);
        SPHComputeShader.SetBuffer(ps_compute_density_kernel, "rearrangedParticles", rearrangedParticlesBuffer);
        SPHComputeShader.SetBuffer(ps_compute_density_kernel, "gridOffsetBuffer", gridOffsetBuffer);
        SPHComputeShader.SetBuffer(ps_compute_density_kernel, "density", density_buffer);
    
        SPHComputeShader.SetBuffer(ps_compute_interact_acceleration_kernel, "particles", particle_buffer);
        SPHComputeShader.SetBuffer(ps_compute_interact_acceleration_kernel, "rearrangedParticles", rearrangedParticlesBuffer);
        SPHComputeShader.SetBuffer(ps_compute_interact_acceleration_kernel, "gridOffsetBuffer", gridOffsetBuffer);
        SPHComputeShader.SetBuffer(ps_compute_interact_acceleration_kernel, "density", density_buffer);
        SPHComputeShader.SetBuffer(ps_compute_interact_acceleration_kernel, "velocity", velocity_buffer);
        SPHComputeShader.SetBuffer(ps_compute_interact_acceleration_kernel, "force", force_buffer);

        SPHComputeShader.SetBuffer(cv_compute_density_kernel, "particles", particle_buffer);
        SPHComputeShader.SetBuffer(cv_compute_density_kernel, "grid", gridBuffer);
        SPHComputeShader.SetBuffer(cv_compute_density_kernel, "particleNeighbors", particleNeighborsBuffer);
        SPHComputeShader.SetBuffer(cv_compute_density_kernel, "density", density_buffer);

        SPHComputeShader.SetBuffer(cv_compute_interact_acceleration_kernel, "particles", particle_buffer);
        SPHComputeShader.SetBuffer(cv_compute_interact_acceleration_kernel, "grid", gridBuffer);
        SPHComputeShader.SetBuffer(cv_compute_interact_acceleration_kernel, "particleNeighbors", particleNeighborsBuffer);
        SPHComputeShader.SetBuffer(cv_compute_interact_acceleration_kernel, "density", density_buffer);
        SPHComputeShader.SetBuffer(cv_compute_interact_acceleration_kernel, "velocity", velocity_buffer);
        SPHComputeShader.SetBuffer(cv_compute_interact_acceleration_kernel, "force", force_buffer);
        
        SPHComputeShader.SetBuffer(compute_external_acceleration_kernel, "particles", particle_buffer);
        SPHComputeShader.SetBuffer(compute_external_acceleration_kernel, "force", force_buffer);
        
        SPHComputeShader.SetBuffer(integrate_kernel, "particles", particle_buffer);
        SPHComputeShader.SetBuffer(integrate_kernel, "velocity", velocity_buffer);
        SPHComputeShader.SetBuffer(integrate_kernel, "force", force_buffer);
        SPHComputeShader.SetBuffer(integrate_kernel, "density", density_buffer);

        SPHComputeShader.SetBuffer(dampen_by_bounds_kernel, "particles", particle_buffer);
        SPHComputeShader.SetBuffer(dampen_by_bounds_kernel, "velocity", velocity_buffer);

        SPHComputeShader.SetBuffer(integrate_boid_particles_kernel, "particles", particle_buffer);
        SPHComputeShader.SetBuffer(integrate_boid_particles_kernel, "velocity", velocity_buffer);
        SPHComputeShader.SetBuffer(integrate_boid_particles_kernel, "boids", boidManager.boidsBuffer);
        SPHComputeShader.SetBuffer(integrate_boid_particles_kernel, "boidVelocities", boidManager.boidVelocitiesBuffer);

    }

    private void Update() {
        if (initialized) UpdateParticles();
    }

    public void UpdateParticles() {

        int debugNumGridCells = Mathf.Min(1000, numGridCells);
        int debugNumParticles = Mathf.Min(500, numParticles);
        UpdateShaderVariables(true);

        // Reset the grid and num neighbors buffer, if we are debugging
        SPHComputeShader.Dispatch(clear_grid_kernel, numBlocks,1,1);
        SPHComputeShader.Dispatch(reset_num_neighbors_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE),1,1);
        
        // Calculate how many particles are in each grid cell, and get each particle's offset for their particular grid cell
        SPHComputeShader.Dispatch(update_grid_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE),1,1);
        if (verbose_grid) {
            DebugBufferInt("Grid", debugNumGridCells, gridBuffer, true);
            DebugBufferParticle("Particles", debugNumParticles, particle_buffer);
        }
        //DebugBufferInt("Projected indices", debugNumParticles, temp_buffer);
        if (verbose_offsets) {
            DebugBufferInt("Particle Offsets", debugNumParticles, particleOffsetsBuffer);
        }

        // We perform the updates for density, pressure, and force/acceleration. 
        // Which compute shader methods to run will be based on our chosen method.
        if (neighborSearchType == NeighborSearchType.PrefixSum) PrefixSumProcess(debugNumGridCells, debugNumParticles);
        else CubeVolumeProcess(debugNumGridCells, debugNumParticles);

        SPHComputeShader.Dispatch(compute_external_acceleration_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE), 1, 1);
        if (verbose_force) DebugBufferFloat3("Forces",debugNumParticles,force_buffer);

        // Integrate over particles, update their positions after taking all force calcualtions into account
        SPHComputeShader.Dispatch(integrate_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE), 1, 1);

        // Make sure that particles are within bounds, limit them if so
        SPHComputeShader.Dispatch(dampen_by_bounds_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE), 1, 1);
        if (verbose_velocity) DebugBufferFloat3("Velocities",debugNumParticles,velocity_buffer);

        // Integrate the boid particles specifically
        SPHComputeShader.Dispatch(integrate_boid_particles_kernel, Mathf.CeilToInt((float)numBoids / (float)_BLOCK_SIZE), 1, 1);
        
        // Render the particles
        particle_material.SetFloat(size_property, particleRenderSize);
        particle_material.SetBuffer(particle_buffer_property, particle_buffer);
        Graphics.DrawMeshInstancedIndirect(particle_mesh, 0, particle_material, new Bounds(Vector3.zero, new Vector3(1000f, 1000f, 1000f)), arg_buffer, castShadows: UnityEngine.Rendering.ShadowCastingMode.Off);
    }

    private void PrefixSumProcess(int debugNumGridCells, int debugNumParticles) {
        SPHComputeShader.Dispatch(prefix_sum_kernel, numBlocks, 1, 1);
        bool swap = false;
        for (int d = 1; d < numBlocks; d *= 2) {
            SPHComputeShader.SetBuffer(sum_blocks_kernel, "gridSumsBufferIn", swap ? gridSumsBuffer1 : gridSumsBuffer2);
            SPHComputeShader.SetBuffer(sum_blocks_kernel, "gridSumsBuffer", swap ? gridSumsBuffer2 : gridSumsBuffer1);
            SPHComputeShader.SetInt("d", d);
            SPHComputeShader.Dispatch(sum_blocks_kernel, Mathf.CeilToInt((float)numBlocks / (float)_BLOCK_SIZE), 1, 1);
            swap = !swap;
        }
        SPHComputeShader.SetBuffer(add_sums_kernel, "gridSumsBufferIn", swap ? gridSumsBuffer1 : gridSumsBuffer2);
        SPHComputeShader.Dispatch(add_sums_kernel, numBlocks, 1, 1);
        if (verbose_prefix_sums) DebugBufferInt("Prefix Sums", debugNumGridCells, gridOffsetBuffer);

        SPHComputeShader.Dispatch(rearrange_particles_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE), 1, 1);
        if (verbose_rearrange) DebugBufferInt("Rearranged particle", debugNumParticles, rearrangedParticlesBuffer);

        SPHComputeShader.Dispatch(count_num_neighbors_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE), 1, 1);
        if (verbose_num_neighbors) DebugBufferInt("Num Neighbors",debugNumParticles,numNeighborsBuffer);

        SPHComputeShader.Dispatch(ps_compute_density_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE), 1, 1);
        if (verbose_density_pressure) {
            DebugBufferFloat("Densities",debugNumParticles,density_buffer);
            DebugBufferFloat("Pressures",debugNumParticles,pressure_buffer);
        }

        SPHComputeShader.Dispatch(ps_compute_interact_acceleration_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE), 1, 1);
        if (verbose_force) DebugBufferFloat3("Forces",debugNumParticles,force_buffer);
    }

    private void CubeVolumeProcess(int debugNumGridCells, int debugNumParticles) {
        SPHComputeShader.Dispatch(cv_compute_density_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE), 1, 1);
        if (verbose_density_pressure) {
            //DebugBufferFloat("Densities",debugNumParticles,density_buffer);
            DebugBufferParticleTouched("Touched by Boid", 15000, particle_buffer);
        }
        SPHComputeShader.Dispatch(cv_compute_interact_acceleration_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE), 1, 1);
    }

    void OnDestroy() {
        particle_buffer.Dispose();
        particleOffsetsBuffer.Dispose();
        particleNeighborsBuffer.Dispose();
        arg_buffer.Dispose();
        density_buffer.Dispose();
        pressure_buffer.Dispose();
        velocity_buffer.Dispose();
        force_buffer.Dispose();     
        gridBuffer.Dispose();
        gridOffsetBuffer.Dispose();
        gridSumsBuffer1.Dispose();
        gridSumsBuffer2.Dispose();
        rearrangedParticlesBuffer.Dispose();
        numNeighborsBuffer.Dispose();
        temp_buffer.Dispose();
    }
    
    private void DebugBufferInt(string debugText, int debugSize, ComputeBuffer b, bool getTotal = false) {
        int[] temp = new int[debugSize];
        b.GetData(temp);
        string top = "", bottom = "";
        int count = 0;
        for(int i = 0; i < debugSize; i++) {
            top += $"{i}\t|";
            bottom += $"{temp[i]}\t|";
            count += temp[i];
        }
        Debug.Log($"{debugText}:\n{top}\n{bottom}");
        if (getTotal) Debug.Log($"Total: {count}");
        temp = null;
    }
    private void DebugBufferFloat(string debugText, int debugSize, ComputeBuffer b) {
        float[] temp = new float[debugSize];
        b.GetData(temp);
        string top = "", bottom = "";
        for(int i = 0; i < debugSize; i++) {
            top += $"{i}\t|";
            bottom += $"{temp[i]}\t|";
        }
        Debug.Log($"{debugText}:\n{top}\n{bottom}");
        temp = null;
    }
    private void DebugBufferFloat3(string debugText, int debugSize, ComputeBuffer b) {
        float3[] temp = new float3[debugSize];
        b.GetData(temp);
        string top = "", bottom = "";
        for(int i = 0; i < debugSize; i++) {
            top += $"{i}\t|";
            bottom += $"{temp[i]}\t|";
        }
        Debug.Log($"{debugText}:\n{top}\n{bottom}");
        temp = null;
    }
    private void DebugBufferParticle(string debugText, int debugSize, ComputeBuffer b) {
        Particle[] temp = new Particle[debugSize];
        b.GetData(temp);
        string top = "";
        string bottom = "";
        for(int i = 0; i < debugSize; i++) {
            top += $"{i}\t|";
            bottom += $"{temp[i].position}\t|";
        }
        Debug.Log($"{debugText} positions:\n{top}\n{bottom}");
        temp = null;
    }
    private void DebugBufferParticleTouched(string debugText, int debugSize, ComputeBuffer b) {
        Particle[] temp = new Particle[numBoids + debugSize];
        b.GetData(temp);
        string top = "";
        string bottom = "";
        for(int i = numBoids; i < numBoids + debugSize; i++) {
            if (temp[i].touchedByBoid > 0) {
                top += $"{i}\t|";
                bottom += $"{temp[i].touchedByBoid}\t|";
            }
        }
        Debug.Log($"{debugText} touched:\n{top}\n{bottom}");
        temp = null;
    }
}