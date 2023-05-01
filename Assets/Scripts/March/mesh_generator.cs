using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(MeshRenderer))]
public class mesh_generator : MonoBehaviour
{
    struct Triangle {
        public Vector3 a;
        public Vector3 b;
        public Vector3 c;

        public Vector3 this [int i] {
            get {
                switch (i) {
                    case 0:
                        return a;
                    case 1:
                        return b;
                    default:
                        return c;
                }
            }
        }
    }

    [Header("== REFERENCES ==")]
    [Tooltip("Reference to the script that handles all particles, aka our ParticleManager")]
    public fluid_gpu fluid_cs;

    [Header("== WORLD CONFIGURATIONS ==")]
        [Tooltip("How big are our grid cells? Recommended to match `smoothingRadius` from `fluid_gpu.cs")]
        public float gridCellSize = 1f;
        [Tooltip("Where in world space are we centering the simulation around?")]
        public float[] origin = {0f,0f,0f};
        [Tooltip("What are the total world space length (per axis) is the simulation?")]
        public float[] bounds = {50f, 50f, 50f};
        [Tooltip("In our grid space, how many buffer cells are added to each axis?")]
        public int[] bufferCells = {10,10,10};
        [ReadOnly, SerializeField, Tooltip("What is the world space size of the simulation, after taking into account grid cell size and buffer cells?")]
        private float[] outerBounds = {60f, 60f, 60f};
        [ReadOnly, SerializeField, Tooltip("Given `outerBounds`, how many grid cells are along each axis?")]
        private int[] _numCellsPerAxis;
        public Vector3Int numCellsPerAxis { get => new Vector3Int(_numCellsPerAxis[0], _numCellsPerAxis[1], _numCellsPerAxis[2]); set {} }
        [ReadOnly, SerializeField, Tooltip("How many grid cells do we have in total?")]
        private int numGridCells;
        [ReadOnly, SerializeField, Tooltip("Given the total number of grid cells, how many voxels do we have per axis?")]
        private int[] numVoxelsPerAxis;
        [ReadOnly, Tooltip("The total number of voxels, given `numVoxelsPerAxis")] public int numVoxels;

    [Header("== TRIANGLE CONFIGURATIONS ==")]
        [ReadOnly, SerializeField, Tooltip("How many triangles can we have, max?")]
        private int maxTriangleCount;
        [Tooltip("How many points will we render across the cube grid, per axis?")]
        public int numPointsPerAxis = 50;
        [ReadOnly, SerializeField] 
        private int numPoints;
        [Tooltip("ISO Level for the marching cubes algorithm")]
        public float isolevel;
        [ReadOnly, SerializeField, Tooltip("The calculated number of triangles outputted by the marching cubes shader.")]
        private int numTriangles = 0;
        // read-only: the triangles outputted by the marching cubes shader
        private Triangle[] triangles = new Triangle[0];

    [Header("== GPU CONFIGURATIOSN ==")]
        // how big each thread group is
        const int _BLOCK_SIZE = 1024;
        // # of threads per axis
        private int[] numThreadsPerAxis = new int[3];

    [Header("Debug Configurations")]
    [Tooltip("When active, prints out details about the triangles outputted by teh marching cubes shader")]
    public bool verbose_triangles = false;

    public density_generator density_gen;
    public ComputeShader shader;
    public Material material;
    public float boundsSize = 1;
    Mesh fluid;
    MeshRenderer fluid_mesh_renderer;
    public ComputeBuffer triangle_buffer;
    public ComputeBuffer point_buffer;
    public ComputeBuffer triangle_count_buffer;
    public ComputeBuffer voxel_density_buffer;

    private int size_property;
    private int particle_buffer_property; 

    private void Awake() {
        // Move to center of the world
        transform.position = Vector3.zero;

        // Initialize shader properties
        size_property = Shader.PropertyToID("size");
        particle_buffer_property = Shader.PropertyToID("particle_buffer");

        // Dunno what these do -_-
        fluid_mesh_renderer = GetComponent<MeshRenderer>();
        fluid = GetComponent<MeshFilter>().mesh;

        // Initialize Particles
        fluid_cs.Initialize();

        InitializeVariables();
        InitializeShaderVariables();
        InitializeBuffers();

        // Initialize Density
        density_gen.Initialize();
    }

    private void InitializeVariables() {
        _numCellsPerAxis = new int[3];
        _numCellsPerAxis[0] = Mathf.CeilToInt(bounds[0] / gridCellSize) + bufferCells[0];
        _numCellsPerAxis[1] = Mathf.CeilToInt(bounds[1] / gridCellSize) + bufferCells[1];
        _numCellsPerAxis[2] = Mathf.CeilToInt(bounds[2] / gridCellSize) + bufferCells[2];
        // Re-Calculate the new bounds based on the number of cells per axis
        outerBounds[0] = (float)_numCellsPerAxis[0] * gridCellSize;
        outerBounds[1] = (float)_numCellsPerAxis[1] * gridCellSize;
        outerBounds[2] = (float)_numCellsPerAxis[2] * gridCellSize;
        // How many grid cells are there in total?
        numGridCells = _numCellsPerAxis[0] * _numCellsPerAxis[1] * _numCellsPerAxis[2];
        // How many voxels are present per axis?
        numVoxelsPerAxis = new int[3] {
            _numCellsPerAxis[0] - 1,
            _numCellsPerAxis[1] - 1,
            _numCellsPerAxis[2] - 1
        };
        // How many voxels do we have in total?
        numVoxels = numVoxelsPerAxis[0] * numVoxelsPerAxis[1] * numVoxelsPerAxis[2];

        // How many triangles can we generate?
        maxTriangleCount = numVoxels * 5;

        // We gotta also set variables for GPU processing
        numThreadsPerAxis = new int[3] {
            Mathf.CeilToInt((float)numVoxelsPerAxis[0]/(float)_BLOCK_SIZE),
            Mathf.CeilToInt((float)numVoxelsPerAxis[1]/(float)_BLOCK_SIZE),
            Mathf.CeilToInt((float)numVoxelsPerAxis[2]/(float)_BLOCK_SIZE)
        };
    }

    private void InitializeShaderVariables() {
        shader.SetInt("n_point_per_axis", numPointsPerAxis);
        shader.SetFloat("isolevel", isolevel);
    }
    

    void Update() {
        // Tell our particle manager to update the particles
        fluid_cs.UpdateParticles();
        UpdateGridDensity();
        UpdateTriangles();
        UpdateMesh(fluid);
    }

    unsafe void InitializeBuffers() {
        triangle_buffer = new ComputeBuffer(maxTriangleCount, sizeof(Triangle), ComputeBufferType.Append);
        triangle_count_buffer = new ComputeBuffer(1, sizeof(int), ComputeBufferType.Raw);
        voxel_density_buffer = new ComputeBuffer(numPoints, sizeof(float));
        point_buffer = new ComputeBuffer(numPoints, sizeof(float) * 3);

        shader.SetBuffer(0, "triangles", triangle_buffer);
        shader.SetBuffer(0, "voxel_density", voxel_density_buffer);
        shader.SetBuffer(0, "particles", fluid_cs.particle_buffer);
        shader.SetBuffer(0, "points", point_buffer);
        
        material.SetFloat(size_property, fluid_cs.particleRenderRadius);
        material.SetBuffer(particle_buffer_property, fluid_cs.particle_buffer);
        fluid_mesh_renderer.material = material;
    }

    private void UpdateGridDensity() {
        Vector3 center = new Vector3(origin[0], origin[1], origin[2]);
        Vector3 worldBounds = new Vector3(outerBounds[0], outerBounds[1], outerBounds[2]);

        density_gen.generate(
            point_buffer, numPointsPerAxis, 
            boundsSize, 
            worldBounds, 
            center,
            Vector3.zero, 
            gridCellSize
        );
        
        /* 
        float[] a = new float[100];
        voxel_density_buffer.GetData(a);
        Debug.Log("voxel");
        for(int i = 0; i < 100; ++i) Debug.Log(a[i]); 
        */

    }

    private void UpdateTriangles() {
        triangle_buffer.SetCounterValue(0);
        shader.Dispatch (0, numThreadsPerAxis[0], numThreadsPerAxis[1], numThreadsPerAxis[2]);
        ComputeBuffer.CopyCount (triangle_buffer, triangle_count_buffer, 0);
        int[] triCountArray = { 0 };
        triangle_count_buffer.GetData (triCountArray);
        numTriangles = triCountArray[0];
        triangles = new Triangle[numTriangles];
        triangle_buffer.GetData(triangles, 0, 0, numTriangles);
        if (verbose_triangles) {
            Debug.Log($"Triangle count: {numTriangles}");
            for(int i = 0; i < numTriangles; i++) {
                Debug.Log($"a: {triangles[i].a} | b: {triangles[i].b} | c: {triangles[i].c}");
            }
        }
    }

    private void UpdateMesh(Mesh mesh) {
        mesh.Clear();
        var vertices = new Vector3[numTriangles * 3];
        var meshTriangles = new int[numTriangles * 3];

        for(int i = 0; i < numTriangles; i++) {
            for(int j = 0; j < 3; j++) {
                meshTriangles[i * 3 + j] = i * 3 + j; /* assign index */
                vertices[i * 3 + j] = triangles[i][j];
            }
        }
        mesh.vertices = vertices;
        mesh.triangles = meshTriangles;
        /* Color[] colors = new Color[vertices.Length];
        for(int i = 0; i < vertices.Length; ++i)
            colors[i] = Color.Lerp(Color.red, Color.green, vertices[i].y);
        mesh.colors = colors; */
        mesh.RecalculateNormals();
        /* Graphics.DrawMeshInstancedIndirect(mesh, 0, material, new Bounds(Vector3.zero, new Vector3(1000f, 1000f, 1000f)), fluid_cs.arg_buffer, castShadows: UnityEngine.Rendering.ShadowCastingMode.Off); */
    }

    void OnDestroy()
    {
        //particle_buffer.Release();
        voxel_density_buffer.Release();
        if(triangle_buffer != null) {
            triangle_buffer.Release();
            triangle_count_buffer.Release();
        }
    }
}