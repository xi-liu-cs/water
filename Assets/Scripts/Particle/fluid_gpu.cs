using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;
using Unity.Mathematics;
using Random = UnityEngine.Random;
using UnityEditor;

public class fluid_gpu : MonoBehaviour
{
    public struct Particle
    {
        public float3 position;
        public int3 grid_indices;
        public int projected_index;
        public int offset;
        public int index;
        public int numNeighbors;
    }

    [Header("== WORLD CONFIGURATIONS ==")]
    public float gridCellSize = 1f;
    public float[] origin = {0f,0f,0f};
    public float[] bounds = {50f, 50f, 50f};
    public int[] bufferCells = {10,10,10};
    private float[] outerBounds = {60f, 60f, 60f};
    private int[] numCellsPerAxis;
    private int numGridCells;
    public float[] g = {0f, -9.81f, 0f};

    [Header("== PARTICLE CONFIGURATIONS ==")]
    public int numParticles = 32;
    public float particleRenderSize = 8f;
    public float particleRenderRadius {
        get => particleRenderSize / 2f;
        set { particleRenderSize = value * 2f; }
    }
    public Mesh particle_mesh;

    [Header("== FLUID MECHANICS ==")]
    public float dt = 0.0008f;
    public float smoothingRadius = 8f;   /* h, smoothing length */
    public float particleMass = 1f;
    public float viscosity_coefficient = 0.01f;
    public float rest_density = 9f;
    public float damping = -1f;
    public float gas_constant = 2000f;

    [Header("== GPU SETTINGS ==")]
    public ComputeShader compute_shader;
    const int _BLOCK_SIZE = 1024;
    private int numBlocks;
    // -- Header Indices --
    int generate_grid_kernel;
    int generate_particles_kernel;
    int clear_grid_kernel;
    int reset_num_neighbors_kernel;
    int update_grid_kernel;
    int prefix_sum_kernel;
    int sum_blocks_kernel;
    int add_sums_kernel;
    int rearrange_particles_kernel;
    int count_num_neighbors_kernel;
    int compute_density_pressure_kernel;
    int compute_force_kernel;
    int integrate_kernel;
    int compute_collisions_kernel;
    int dampen_by_bounds_kernel;
    // int clear_hash_grid_kernel;
    // int compute_hash_grid_kernel;
    // int compute_neighbor_list_kernel;
    // int compute_sums_kernel;
    // int march_kernel;

    [Header("== DEBUG CONFIGURATIONS ==")]
    public bool show_gizmos = false;
    public bool show_velocities = false;
    public bool show_forces = false;
    public Color gizmos_particle_color = Color.blue;
    public bool verbose_grid = false;
    public bool verbose_prefix_sums = false;
    public bool verbose_rearrange = false;
    public bool verbose_num_neighbors = false;
    //public bool verbose_neighbors = false;
    public bool verbose_density_pressure = false;
    public bool verbose_force = false;
    public bool verbose_velocity = false;
    public bool verbose_collisions = false;





    [Header("particle")]
    public Mesh fluid_mesh;
    public Material material;
    //public float grid_size = 8f; /* 4 * smoothingRadius */
    
    public Vector3 position_offset = new Vector3(-140, -230, -60);
    public Vector3 velocity_initial = new Vector3(0, 10, 0);
    public float[] bound = {-160, 140, -250, -150, -100, 100};
    [HideInInspector]
    public float radius2,
    radius3,
    radius4,
    radius5,
    mass2;
    Vector3[] points;
    float[] noise_densities;
    int[] triangles;
    triangle[] march_triangles;

    [Header("fluid")]
    public int max_particles_per_grid = 12,
    intial_particles_per_grid = 5;
    public Vector3Int dimension;
    int dimension2,
    dimension3;
    public int[] dimension_array;

    [Header("voxel")]
    public int n_point_per_axis = 50; /* z * n_point_per_axis * n_point_per_axis + y * n_point_per_axis + x, x ^ 3 + x ^ 2 + x = numParticles, (x = 50, numParticles = 130000) */
    public float isolevel = 8;


    public struct triangle
    {
        public Vector3 vertex_a,
        vertex_b,
        vertex_c;
    }

    public struct grid_cell
    {
        public Vector3[] vertex;
        public float[] value;
        public grid_cell(Vector3[] vert, float[] val)
        {
            vertex = vert;
            value = val;
        }
    }

    int size_property = Shader.PropertyToID("size"),
    particle_buffer_property = Shader.PropertyToID("particle_buffer");

    int n = 8,
    n_bound = 6,
    thread_group_size,
    grid_size_over_2,
    n_debug = 64;
    Particle[] particles;
    int[] neighbor_list,
    neighbor_tracker,
    int_debug;
    uint[] hash_grid,
    hash_grid_tracker;
    float[] density,
    pressure,
    float_debug;
    Vector3[] velocity,
    force;

    [Header("== VARIABLES ==")]
    float[] kernel_sums;

    [Header("== COMPUTE BUFFERS ==")]
    public ComputeBuffer particle_buffer;
    public ComputeBuffer arg_buffer;
    public ComputeBuffer neighbor_list_buffer;
    public ComputeBuffer neighbor_tracker_buffer;
    public ComputeBuffer hash_grid_buffer;
    public ComputeBuffer hash_grid_tracker_buffer;
    public ComputeBuffer density_buffer;
    public ComputeBuffer pressure_buffer;
    public ComputeBuffer velocity_buffer;
    public ComputeBuffer force_buffer;
    public ComputeBuffer bound_buffer;
    public ComputeBuffer point_buffer;
    public ComputeBuffer noise_density_buffer;
    public ComputeBuffer triangle_buffer;
    public ComputeBuffer triangle_count_buffer;
    public ComputeBuffer int_debug_buffer;
    public ComputeBuffer float_debug_buffer;

    private ComputeBuffer gridBuffer;
    private ComputeBuffer gridOffsetBuffer;
    private ComputeBuffer gridSumsBuffer1;
    private ComputeBuffer gridSumsBuffer2;
    private ComputeBuffer rearrangedParticlesBuffer;
    private ComputeBuffer numNeighborsBuffer;
    private ComputeBuffer storedVelocitiesBuffer;

    public march_table table;

    public float max_density = 0;


    void OnDrawGizmos() {
        Vector3 _origin = new Vector3(origin[0], origin[1], origin[2]);
        Gizmos.color = Color.white;
        Gizmos.DrawWireCube(_origin, new Vector3(bounds[0], bounds[1], bounds[2]));
        Gizmos.color = Color.yellow;
        Gizmos.DrawWireCube(_origin, new Vector3(outerBounds[0], outerBounds[1], outerBounds[2]));
        Gizmos.color = Color.red;
        Gizmos.DrawSphere(_origin,0.1f);

        if (!Application.isPlaying || !show_gizmos) return;

        Vector4 gridColor = new Vector4(1f,0f,0f,0.25f);
        Vector3 gridCellSize3D = new Vector3(gridCellSize, gridCellSize, gridCellSize);
        Particle[] temp_particles = new Particle[numParticles];
        particle_buffer.GetData(temp_particles);
        int[] temp_num_neighbors = new int[numParticles];
        numNeighborsBuffer.GetData(temp_num_neighbors);
        float3[] temp_forces_array = new float3[numParticles];
        force_buffer.GetData(temp_forces_array);
        float3[] temp_velocities_array = new float3[numParticles];
        velocity_buffer.GetData(temp_velocities_array);
        for(int i = 1; i < numParticles; i++) {
            Gizmos.color = gizmos_particle_color;
            Gizmos.DrawSphere(temp_particles[i].position, particleRenderRadius);
            Gizmos.color = gridColor;
            Gizmos.DrawCube(GetGridCellWorldPositionFromXYZIndices(temp_particles[i].grid_indices), gridCellSize3D);
            if (show_forces) {
                Gizmos.color = Color.green;
                Gizmos.DrawLine(temp_particles[i].position, temp_particles[i].position + temp_forces_array[i]);
            }
            if (show_velocities) {
                Gizmos.color = Color.blue;
                Gizmos.DrawLine(temp_particles[i].position, temp_particles[i].position + temp_velocities_array[i]);
            }
        }
        Gizmos.color = Color.yellow;
        Gizmos.DrawSphere(temp_particles[0].position, particleRenderRadius);
        Gizmos.color = gridColor;
        Gizmos.DrawCube(GetGridCellWorldPositionFromXYZIndices(temp_particles[0].grid_indices), gridCellSize3D);
        Gizmos.color = Color.yellow;
        Gizmos.DrawWireSphere(temp_particles[0].position, smoothingRadius);
        Vector3 handlePos = new Vector3(
            temp_particles[0].position[0],
            temp_particles[0].position[1] + particleRenderSize*2f,
            temp_particles[0].position[2] + particleRenderSize*2f
        );
        Handles.Label(handlePos, temp_num_neighbors[0].ToString());

    }

    public void Awake() {
        InitializeVariables();

        //InitializeParticles();
        InitializeKernels();
        InitializeShaderVariables();
        InitializeBuffers();

        compute_shader.Dispatch(generate_grid_kernel, numBlocks,1,1);
        compute_shader.Dispatch(generate_particles_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE),1,1);
        /*
        Particle[] temp_particles = new Particle[numParticles];
        particle_buffer.GetData(temp_particles);
        string top = "";
        string bottom = "";
        for(int i = 0; i < numParticles; i++) {
            top += $"{i}\t|";
            bottom += $"{temp_particles[i].index}\t|";
        }
        Debug.Log($"Particle Indexes:\n{top}\n{bottom}");
        */

        /* Debug.Log("group_size " + thread_group_size);
        compute_shader.Dispatch(noise_density_kernel, thread_group_size, thread_group_size, thread_group_size);
        Vector4[] cpu = new Vector4[numParticles];
        point_buffer.GetData(cpu);
        Debug.Log("cpu");
        for(int i = 0; i < cpu.Length; ++i)
            Debug.Log(cpu[i]); */
    }

    private void InitializeVariables() {
        // Calculate how many cells will fit within the provided dimensions
        numCellsPerAxis = new int[3];
        numCellsPerAxis[0] = Mathf.CeilToInt(bounds[0] / gridCellSize) + bufferCells[0];
        numCellsPerAxis[1] = Mathf.CeilToInt(bounds[1] / gridCellSize) + bufferCells[1];
        numCellsPerAxis[2] = Mathf.CeilToInt(bounds[2] / gridCellSize) + bufferCells[2];
        // Re-Calculate the new bounds based on the number of cells per axis
        outerBounds[0] = (float)numCellsPerAxis[0] * gridCellSize;
        outerBounds[1] = (float)numCellsPerAxis[1] * gridCellSize;
        outerBounds[2] = (float)numCellsPerAxis[2] * gridCellSize;
        // How many grid cells are there in total?
        numGridCells = numCellsPerAxis[0] * numCellsPerAxis[1] * numCellsPerAxis[2];

        // == GPU SETTINGS ==
        numBlocks = Mathf.CeilToInt((float)numGridCells / (float)_BLOCK_SIZE);

        // == OLD SETTINGS ==
        grid_size_over_2 = (int)gridCellSize / 2;
        dimension = new Vector3Int((int)((bound[1] - bound[0]) / gridCellSize), (int)((bound[3] - bound[2]) / gridCellSize), (int)((bound[5] - bound[4]) / gridCellSize));
        dimension_array = new int[]{dimension.x, dimension.y, dimension.z};
        dimension2 = dimension.x * dimension.y;
        dimension3 = dimension2 * dimension.z;
        table = new march_table();
        radius2 = smoothingRadius * smoothingRadius;
        radius3 = radius2 * smoothingRadius;
        radius4 = radius3 * smoothingRadius;
        radius5 = radius4 * smoothingRadius;
        mass2 = particleMass * particleMass;
        //thread_group_size = numParticles / dimension.x;

        Debug.Log($"{numCellsPerAxis[0]},{numCellsPerAxis[1]},{numCellsPerAxis[2]}");
        Debug.Log($"{numGridCells}");
    }

    void InitializeKernels() {
        generate_grid_kernel = compute_shader.FindKernel("GenerateGrid");
        generate_particles_kernel = compute_shader.FindKernel("GenerateParticles");
        clear_grid_kernel = compute_shader.FindKernel("ClearGrid");
        reset_num_neighbors_kernel = compute_shader.FindKernel("ResetNumNeighbors");
        update_grid_kernel = compute_shader.FindKernel("UpdateGridCellCounts");
        prefix_sum_kernel = compute_shader.FindKernel("PrefixSum");
        sum_blocks_kernel = compute_shader.FindKernel("SumBlocks");
        add_sums_kernel = compute_shader.FindKernel("AddSums");
        rearrange_particles_kernel = compute_shader.FindKernel("RearrangeParticles");
        count_num_neighbors_kernel = compute_shader.FindKernel("CountNumNeighbors");
        compute_density_pressure_kernel = compute_shader.FindKernel("ComputeDensityPressure");
        compute_force_kernel = compute_shader.FindKernel("ComputeForce");
        integrate_kernel = compute_shader.FindKernel("Integrate");
        compute_collisions_kernel = compute_shader.FindKernel("ComputeCollisions");
        dampen_by_bounds_kernel = compute_shader.FindKernel("DampenByBounds");
    }

    void InitializeShaderVariables() {
        // == WORLD CONFIGURATIONS ==
        compute_shader.SetFloat("gridCellSize", gridCellSize);
        compute_shader.SetInts("numCellsPerAxis", numCellsPerAxis);
        compute_shader.SetInt("total_number_of_cells", numGridCells);
        compute_shader.SetFloats("origin", origin);
        compute_shader.SetFloats("bounds", bounds);
        compute_shader.SetFloats("outerBounds", outerBounds);
        compute_shader.SetFloats("g", g);

        // == PARTICLE CONFIGURATIONS ==
        compute_shader.SetInt("numParticles", numParticles);
        compute_shader.SetFloat("particleRenderRadius", particleRenderRadius);

        // == FLUID MECHANICS ==
        compute_shader.SetFloat("smoothingRadius", smoothingRadius);
        compute_shader.SetFloat("radius2", radius2);
        compute_shader.SetFloat("radius3", radius3);
        compute_shader.SetFloat("radius4", radius4);
        compute_shader.SetFloat("radius5", radius5);
        compute_shader.SetFloat("radius6", radius5 * smoothingRadius);
        compute_shader.SetFloat("radius9", radius5 * radius4);
        compute_shader.SetFloat("particleMass", particleMass);
        compute_shader.SetFloat("mass2", mass2);
        compute_shader.SetFloat("gas_constant", gas_constant);
        compute_shader.SetFloat("rest_density", rest_density);
        compute_shader.SetFloat("viscosity_coefficient", viscosity_coefficient);
        compute_shader.SetFloat("damping", damping);
        compute_shader.SetInts("dimension", dimension_array);

        // == GPU SETTINGS
        compute_shader.SetInt("numBlocks", numBlocks);
        compute_shader.SetInt("randomSeed", Random.Range(0,int.MaxValue));

        // == CONST SETTINGS ==
        compute_shader.SetFloat("epsilon", Mathf.Epsilon);
        compute_shader.SetFloat("pi", Mathf.PI);

        // == OLD SETTINGS ==
        compute_shader.SetInt("n", 8);
        compute_shader.SetInt("max_particles_per_grid", max_particles_per_grid);
        compute_shader.SetFloat("e", Mathf.Exp(1));
        compute_shader.SetVector("time", Shader.GetGlobalVector("_Time"));
        compute_shader.SetInt("n_point_per_axis", n_point_per_axis);
        compute_shader.SetFloat("isolevel", isolevel);
    }

    public unsafe void InitializeBuffers() {
        uint[] arg = {particle_mesh.GetIndexCount(0), (uint)numParticles, particle_mesh.GetIndexStart(0), particle_mesh.GetBaseVertex(0), 0};
        arg_buffer = new ComputeBuffer(1, arg.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        arg_buffer.SetData(arg);

        particle_buffer = new ComputeBuffer(numParticles, sizeof(float)*3 + sizeof(int)*7);
        density_buffer = new ComputeBuffer(numParticles, sizeof(float));
        pressure_buffer = new ComputeBuffer(numParticles, sizeof(float));
        velocity_buffer = new ComputeBuffer(numParticles, 3 * sizeof(float));
        force_buffer = new ComputeBuffer(numParticles, 3 * sizeof(float));
        gridBuffer = new ComputeBuffer(numGridCells, sizeof(int));
        gridOffsetBuffer = new ComputeBuffer(numGridCells, sizeof(int));
        gridSumsBuffer1 = new ComputeBuffer(numBlocks, sizeof(int));
        gridSumsBuffer2 = new ComputeBuffer(numBlocks, sizeof(int));
        rearrangedParticlesBuffer = new ComputeBuffer(numParticles, sizeof(float)*3 + sizeof(int)*7);
        numNeighborsBuffer = new ComputeBuffer(numParticles, sizeof(int));
        storedVelocitiesBuffer = new ComputeBuffer(numParticles, sizeof(float)*3);

        compute_shader.SetBuffer(generate_grid_kernel, "grid", gridBuffer);

        compute_shader.SetBuffer(reset_num_neighbors_kernel, "numNeighbors", numNeighborsBuffer);

        compute_shader.SetBuffer(generate_particles_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(generate_particles_kernel, "density", density_buffer);
        compute_shader.SetBuffer(generate_particles_kernel, "pressure", pressure_buffer);
        compute_shader.SetBuffer(generate_particles_kernel, "force", force_buffer);
        compute_shader.SetBuffer(generate_particles_kernel, "velocity", velocity_buffer);

        compute_shader.SetBuffer(clear_grid_kernel, "grid", gridBuffer);
        
        compute_shader.SetBuffer(update_grid_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(update_grid_kernel, "grid", gridBuffer);

        compute_shader.SetBuffer(prefix_sum_kernel, "gridOffsetBufferIn", gridBuffer);
        compute_shader.SetBuffer(prefix_sum_kernel, "gridOffsetBuffer", gridOffsetBuffer);
        compute_shader.SetBuffer(prefix_sum_kernel, "gridSumsBuffer", gridSumsBuffer2);

        compute_shader.SetBuffer(add_sums_kernel, "gridOffsetBuffer", gridOffsetBuffer);

        compute_shader.SetBuffer(rearrange_particles_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(rearrange_particles_kernel, "rearrangedParticles", rearrangedParticlesBuffer);
        compute_shader.SetBuffer(rearrange_particles_kernel, "gridOffsetBuffer", gridOffsetBuffer);

        compute_shader.SetBuffer(count_num_neighbors_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(count_num_neighbors_kernel, "rearrangedParticles", rearrangedParticlesBuffer);
        compute_shader.SetBuffer(count_num_neighbors_kernel, "gridOffsetBuffer", gridOffsetBuffer);
        compute_shader.SetBuffer(count_num_neighbors_kernel, "numNeighbors", numNeighborsBuffer);

        compute_shader.SetBuffer(compute_density_pressure_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_density_pressure_kernel, "rearrangedParticles", rearrangedParticlesBuffer);
        compute_shader.SetBuffer(compute_density_pressure_kernel, "gridOffsetBuffer", gridOffsetBuffer);
        compute_shader.SetBuffer(compute_density_pressure_kernel, "density", density_buffer);
        compute_shader.SetBuffer(compute_density_pressure_kernel, "pressure", pressure_buffer);

        compute_shader.SetBuffer(compute_force_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_force_kernel, "rearrangedParticles", rearrangedParticlesBuffer);
        compute_shader.SetBuffer(compute_force_kernel, "gridOffsetBuffer", gridOffsetBuffer);
        compute_shader.SetBuffer(compute_force_kernel, "density", density_buffer);
        compute_shader.SetBuffer(compute_force_kernel, "pressure", pressure_buffer);
        compute_shader.SetBuffer(compute_force_kernel, "velocity", velocity_buffer);
        compute_shader.SetBuffer(compute_force_kernel, "force", force_buffer);

        compute_shader.SetBuffer(integrate_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(integrate_kernel, "velocity", velocity_buffer);
        compute_shader.SetBuffer(integrate_kernel, "force", force_buffer);
        compute_shader.SetBuffer(integrate_kernel, "density", density_buffer);
        compute_shader.SetBuffer(integrate_kernel, "storedVelocities", storedVelocitiesBuffer);

        compute_shader.SetBuffer(compute_collisions_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_collisions_kernel, "rearrangedParticles", rearrangedParticlesBuffer);
        compute_shader.SetBuffer(compute_collisions_kernel, "gridOffsetBuffer", gridOffsetBuffer);
        compute_shader.SetBuffer(compute_collisions_kernel, "velocity", velocity_buffer);
        compute_shader.SetBuffer(compute_collisions_kernel, "storedVelocities", storedVelocitiesBuffer);

        compute_shader.SetBuffer(dampen_by_bounds_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(dampen_by_bounds_kernel, "velocity", velocity_buffer);
    }

    int calculate_cube_index(grid_cell cell)
    {
        int cube_index = 0;
        for(int i = 0; i < 8; ++i)
            if(cell.value[i] < isolevel)
                cube_index |= 1 << i;
        return cube_index;
    }

    Vector3 interpolate(Vector3 v1, float val1, Vector3 v2, float val2)
    {
        Vector3 interpolated;
        float mu = (isolevel - val1) / (val2 - val1);
        interpolated.x = mu * (v2.x - v1.x) + v1.x;
        interpolated.y = mu * (v2.y - v1.y) + v1.y;
        interpolated.z = mu * (v2.z - v1.z) + v1.z;
        return interpolated;
    }

    Vector3[] get_intersection_coordinates(grid_cell cell)
    {
        Vector3[] intersections = new Vector3[12];
        int cube_index = calculate_cube_index(cell);
        int intersection_key = table.edge_table[cube_index];
        int i = 0;
        while(intersection_key != 0)
        {
            if((intersection_key & 1) != 0)
            {
                int v1 = table.edge_to_vertex[i, 0], v2 = table.edge_to_vertex[i, 1];
                Vector3 intersection_point = interpolate(cell.vertex[v1], cell.value[v1], cell.vertex[v2], cell.value[v2]);
                intersections[i] = intersection_point;
            }
            ++i;
            intersection_key >>= 1;
        }
        return intersections;
    }

    List<triangle> get_triangles(Vector3[] intersections, int cube_index)
    {
        List<triangle> a = new List<triangle>();
        for(int i = 0; table.triangle_table[cube_index, i] != -1; i += 3)
        {
            triangle t;
            t.vertex_a = intersections[table.triangle_table[cube_index, i]];
            t.vertex_b = intersections[table.triangle_table[cube_index, i + 1]];
            t.vertex_c = intersections[table.triangle_table[cube_index, i + 2]];
            a.Add(t);
        }
        return a;
    }

    void print_triangles(List<triangle> a)
    {
        for(int i = 0; i < a.Count; ++i)
        {
            Debug.Log(a[i].vertex_a);
            Debug.Log(a[i].vertex_b);
            Debug.Log(a[i].vertex_c);
        }
    }

    List<triangle> triangulate_cell(grid_cell cell)
    {
        int cube_index = calculate_cube_index(cell);
        Vector3[] intersections = get_intersection_coordinates(cell);
        List<triangle> a = get_triangles(intersections, cube_index);
        return a;
    }

    int index_from_coordinate(int x, int y, int z)
    {
        return z * n_point_per_axis * n_point_per_axis + y * n_point_per_axis + x;
    }

    int[] triangulate_field(Vector3[] points)
    {
        int[] a = new int[numParticles];
        int idx = 0;
        for(int i = 0; i + 1 < n_point_per_axis; ++i)
        {
            for(int j = 0; j + 1 < n_point_per_axis; ++j)
            {
                for(int k = 0; k + 1 < n_point_per_axis; ++k)
                {
                    float x = i, y = j, z = k;
                    if(index_from_coordinate(i + 1, j + 1, k + 1) >= numParticles)
                        break;
                    grid_cell cell = new grid_cell
                    (
                        new Vector3[]
                        {
                            new Vector3(x, y, z), new Vector3(x + 1f, y, z),
                            new Vector3(x + 1f, y, z + 1f), new Vector3(x, y, z + 1f),
                            new Vector3(x, y + 1f, z), new Vector3(x + 1f, y + 1f, z),
                            new Vector3(x + 1f, y + 1f, z + 1f), new Vector3(x, y + 1f, z + 1f)
                        },
                        new float[]
                        {
                            noise_densities[index_from_coordinate(i, j, k)], noise_densities[index_from_coordinate(i + 1, j, k)],
                            noise_densities[index_from_coordinate(i + 1, j, k + 1)], noise_densities[index_from_coordinate(i, j, k + 1)],
                            noise_densities[index_from_coordinate(i, j + 1, k)], noise_densities[index_from_coordinate(i + 1, j + 1, k)],
                            noise_densities[index_from_coordinate(i + 1, j + 1, k + 1)], noise_densities[index_from_coordinate(i, j + 1, k + 1)]
                        }
                    );
                    List<triangle> cell_triangles = triangulate_cell(cell);
                    for(int l = 0; l < cell_triangles.Count; ++l)
                    {
                        triangle tri = cell_triangles[l];
                        a[idx++] = (int)tri.vertex_a[0];
                        a[idx++] = (int)tri.vertex_a[1];
                        a[idx++] = (int)tri.vertex_a[2];
                        a[idx++] = (int)tri.vertex_b[0];
                        a[idx++] = (int)tri.vertex_b[1];
                        a[idx++] = (int)tri.vertex_b[2];
                        a[idx++] = (int)tri.vertex_c[0];
                        a[idx++] = (int)tri.vertex_c[1];
                        a[idx++] = (int)tri.vertex_c[2];
                    }
                }
            }
        }
        return a;
    }

    void InitializeParticles()
    {
        particles = new Particle[numParticles];
        density = new float[numParticles];
        pressure = new float[numParticles];
        velocity = new Vector3[numParticles];
        force = new Vector3[numParticles];
        int i = 0;
        while(i < numParticles)
        {
            for(int x = 0; x < dimension.x; ++x)
                for(int y = 0; y < dimension.y; ++y)
                    for(int z = 0; z < dimension.z; ++z)
                    {
                        Vector3 grid_pos = new Vector3(
                            bound[0] + x * gridCellSize, 
                            bound[2] + y * gridCellSize, 
                            bound[4] + z * gridCellSize
                        );
                        for(int a = 0; a < intial_particles_per_grid; ++a)
                        {
                            //Vector3 pos = grid_pos - new Vector3(gridCellSize_over_2, gridCellSize_over_2, gridCellSize_over_2) - new Vector3(UnityEngine.Random.Range(0f, gridCellSize_over_2), UnityEngine.Random.Range(0f, gridCellSize_over_2), UnityEngine.Random.Range(0f, gridCellSize_over_2));
                            Vector3 pos = new Vector3(
                                Random.Range(bound[0],bound[1]),
                                Random.Range(bound[2],bound[3]),
                                Random.Range(bound[4],bound[5])
                            );
                            Vector3Int grid_indices = GetGridXYZIndices(pos);
                            int projected_index = GetProjectedGridIndexFromXYZ(grid_indices);
                            particles[i] = new Particle{
                                position = new(pos.x, pos.y, pos.z), 
                                grid_indices = new(grid_indices.x, grid_indices.y, grid_indices.z),
                                projected_index = projected_index,
                                index = i
                            };
                            density[i] = -1;
                            pressure[i] = 0;
                            force[i] = Vector3.zero;
                            velocity[i] = velocity_initial;
                            if(++i == numParticles) return;
                        }
                    }
        }
    }
    public Vector3Int GetGridXYZIndices(Vector3 position) {
        return new Vector3Int(
            Mathf.FloorToInt((position.x - bound[0])/gridCellSize),
            Mathf.FloorToInt((position.y - bound[2])/gridCellSize),
            Mathf.FloorToInt((position.z - bound[4])/gridCellSize)
        );
    }
    public int GetProjectedGridIndexFromXYZ(int x, int y, int z) {
        return x + (dimension_array[1] * y) + (dimension_array[1] * dimension_array[2] * z);
    }
    public int GetProjectedGridIndexFromXYZ(Vector3Int xyz) {
        return xyz.x + (dimension_array[1] * xyz.y) + (dimension_array[1] * dimension_array[2] * xyz.z);
    }
    public Vector3 GetGridCellWorldPositionFromXYZIndices(int3 xyz) {
        return new Vector3(
            (origin[0] - ((numCellsPerAxis[0] * gridCellSize)/2f)) + (xyz[0] * gridCellSize) + (gridCellSize/2f),
            (origin[1] - ((numCellsPerAxis[1] * gridCellSize)/2f)) + (xyz[1] * gridCellSize) + (gridCellSize/2f),
            (origin[2] - ((numCellsPerAxis[2] * gridCellSize)/2f)) + (xyz[2] * gridCellSize) + (gridCellSize/2f)
        );
    }

    public void Update() {

        if (dt > 0) compute_shader.SetFloat("dt", dt);
        else compute_shader.SetFloat("dt", Time.deltaTime);
        compute_shader.SetFloat("particleRenderRadius", particleRenderRadius);
        compute_shader.SetFloats("g", g);

        string top = "";
        string bottom = "";

        compute_shader.Dispatch(clear_grid_kernel, numBlocks,1,1);
        compute_shader.Dispatch(reset_num_neighbors_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE),1,1);
        compute_shader.Dispatch(update_grid_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE),1,1);

        if (verbose_grid) {
            uint[] temp_grid = new uint[numGridCells];
            gridBuffer.GetData(temp_grid);
            for(int i = 0; i < numGridCells; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_grid[i]}\t|";
            }
            Debug.Log($"Grid\n{top}\n{bottom}");
            Particle[] temp_particles = new Particle[numParticles];
            particle_buffer.GetData(temp_particles);
            top = "";
            bottom = "";
            for(int i = 0; i < numParticles; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_particles[i].position}\t|";
            }
            Debug.Log($"Particles World Positions:\n{top}\n{bottom}");
            top = "";
            bottom = "";
            for(int i = 0; i < numParticles; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_particles[i].grid_indices}\t|";
            }
            Debug.Log($"Particles Grid Indices:\n{top}\n{bottom}");
            top = "";
            bottom = "";
            for(int i = 0; i < numParticles; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_particles[i].projected_index}\t|";
            }
            Debug.Log($"Particles Projected Indices:\n{top}\n{bottom}");
            top = "";
            bottom = "";
            for(int i = 0; i < numParticles; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_particles[i].offset}\t|";
            }
            Debug.Log($"Particles Offsets:\n{top}\n{bottom}");
        }

        compute_shader.Dispatch(prefix_sum_kernel, numBlocks, 1, 1);
        bool swap = false;
        for (int d = 1; d < numBlocks; d *= 2) {
            compute_shader.SetBuffer(sum_blocks_kernel, "gridSumsBufferIn", swap ? gridSumsBuffer1 : gridSumsBuffer2);
            compute_shader.SetBuffer(sum_blocks_kernel, "gridSumsBuffer", swap ? gridSumsBuffer2 : gridSumsBuffer1);
            compute_shader.SetInt("d", d);
            compute_shader.Dispatch(sum_blocks_kernel, Mathf.CeilToInt((float)numBlocks / (float)_BLOCK_SIZE), 1, 1);
            swap = !swap;
        }
        compute_shader.SetBuffer(add_sums_kernel, "gridSumsBufferIn", swap ? gridSumsBuffer1 : gridSumsBuffer2);
        compute_shader.Dispatch(add_sums_kernel, numBlocks, 1, 1);
        if (verbose_prefix_sums) {
            int[] temp_offset_grid = new int[numGridCells];
            gridOffsetBuffer.GetData(temp_offset_grid);
            top = "";
            bottom="";
            for(int k = 0; k < numGridCells; k++) {
                top += $"{k}\t|";
                bottom += $"{temp_offset_grid[k]}\t|";
            }
            Debug.Log("PREFIX SUMS:\n"+top+"\n"+bottom);
        }

        compute_shader.Dispatch(rearrange_particles_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE), 1, 1);
        if (verbose_rearrange) {
            Particle[] temp_rearranged = new Particle[numParticles];
            rearrangedParticlesBuffer.GetData(temp_rearranged);
            top = "";
            bottom="";
            for(int k = 0; k < numParticles; k++) {
                top += $"{k}\t|";
                bottom += $"{temp_rearranged[k].index}\t|";
            }
            Debug.Log("REARRANGED PARTICLES:\n"+top+"\n"+bottom);
        }

        compute_shader.Dispatch(count_num_neighbors_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE), 1, 1);
        if (verbose_num_neighbors) {
            Particle[] temp_particles = new Particle[numParticles];
            particle_buffer.GetData(temp_particles);
            top = "";
            bottom="";
            for(int k = 0; k < numParticles; k++) {
                if (temp_particles[k].numNeighbors == 0) continue;
                top += $"{k}\t|";
                bottom += $"{temp_particles[k].numNeighbors}\t|";
            }
            Debug.Log("NUM PARTICLE NEIGHBORS:\n"+top+"\n"+bottom);
            /*
            int[] temp_num_neighbors = new int[numParticles];
            numNeighborsBuffer.GetData(temp_num_neighbors);
            top = "";
            bottom="";
            for(int k = 0; k < numParticles; k++) {
                top += $"{k}\t|";
                bottom += $"{temp_num_neighbors[k]}\t|";
            }
            Debug.Log("NUM PARTICLE NEIGHBORS 2:\n"+top+"\n"+bottom);
            */
        }

        compute_shader.Dispatch(compute_density_pressure_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE), 1, 1);
        if (verbose_density_pressure) {
            float[] temp_densities = new float[numParticles];
            float[] temp_pressures = new float[numParticles];
            density_buffer.GetData(temp_densities);
            top = "";
            bottom = "";
            for(int i = 0; i < numParticles; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_densities[i]}\t|";
            }
            Debug.Log($"Densities:\n{top}\n{bottom}");
            pressure_buffer.GetData(temp_pressures);
            top = "";
            bottom = "";
            for(int i = 0; i < numParticles; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_pressures[i]}\t|";
            }
            Debug.Log($"Pressures:\n{top}\n{bottom}");
        }

        compute_shader.Dispatch(compute_force_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE), 1, 1);
        if (verbose_force) {
            float3[] temp_forces = new float3[numParticles];
            force_buffer.GetData(temp_forces);
            top = "";
            bottom = "";
            for(int i = 0; i < numParticles; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_forces[i]}\t|";
            }
            Debug.Log($"Forces:\n{top}\n{bottom}");
        }

        compute_shader.Dispatch(integrate_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE), 1, 1);
        if (verbose_velocity) {
            float3[] temp_velocities = new float3[numParticles];
            velocity_buffer.GetData(temp_velocities);
            top = "";
            bottom = "";
            for(int i = 0; i < numParticles; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_velocities[i]}\t|";
            }
            Debug.Log($"Velocities:\n{top}\n{bottom}");
        }

        /*
        compute_shader.Dispatch(compute_collisions_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE), 1, 1);
        if (verbose_collisions) {
            float3[] temp_velocities = new float3[numParticles];
            velocity_buffer.GetData(temp_velocities);
            top = "";
            bottom = "";
            for(int i = 0; i < numParticles; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_velocities[i]}\t|";
            }
            Debug.Log($"Velocities After Collisions:\n{top}\n{bottom}");
        }
        */

        compute_shader.Dispatch(dampen_by_bounds_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE), 1, 1);

        /*
        if (verbose_hash_grid) {
            uint[] temp_hash_grid = new uint[dimension3 * max_particles_per_grid];
            uint[] temp_hash_grid_tracker = new uint[dimension3];        
            hash_grid_buffer.GetData(temp_hash_grid);
            hash_grid_tracker_buffer.GetData(temp_hash_grid_tracker);
            for(int i = 0; i < dimension3 * max_particles_per_grid; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_hash_grid[i]}\t|";
            }
            Debug.Log($"Hash Grid\n{top}\n{bottom}");
            top = "";
            bottom = "";
            for(int i = 0; i < dimension3; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_hash_grid_tracker[i]}\t|";
            }
            Debug.Log($"Hash Grid Tracker\n{top}\n{bottom}");
            Particle[] temp_particles = new Particle[numParticles];
            particle_buffer.GetData(temp_particles);
            top = "";
            bottom = "";
            for(int i = 0; i < numParticles; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_particles[i].projected_index}\t|";
            }
            Debug.Log($"Particles Projected Indices:\n{top}\n{bottom}");
            top = "";
            bottom = "";
            for(int i = 0; i < numParticles; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_particles[i].offset}\t|";
            }
            Debug.Log($"Particles Offsets:\n{top}\n{bottom}");
        }
        */

        /*
        compute_shader.Dispatch(compute_neighbor_list_kernel, thread_group_size, 1, 1);
        
        if (verbose_neighbors) {
            int[] temp_neighbor_list = new int[numParticles * max_particles_per_grid * n];
            int[] temp_neighbor_tracker = new int[numParticles];
            neighbor_list_buffer.GetData(temp_neighbor_list);
            neighbor_tracker_buffer.GetData(temp_neighbor_tracker);
            top = "";
            bottom = "";
            for(int i = 0; i < numParticles * max_particles_per_grid * n; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_neighbor_list[i]}\t|";
            }
            Debug.Log($"Neighbor list\n{top}\n{bottom}");
            top = "";
            bottom = "";
            for(int i = 0; i < numParticles; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_neighbor_tracker[i]}\t|";
            }
            Debug.Log($"Neighbor Tracker\n{top}\n{bottom}");
        }
        compute_shader.Dispatch(compute_density_pressure_kernel, thread_group_size, 1, 1);
        if (verbose_density_pressure) {
            float[] temp_densities = new float[numParticles];
            float[] temp_pressures = new float[numParticles];
            density_buffer.GetData(temp_densities);
            top = "";
            bottom = "";
            for(int i = 0; i < numParticles; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_densities[i]}\t|";
            }
            Debug.Log($"Densities:\n{top}\n{bottom}");
            pressure_buffer.GetData(temp_pressures);
            top = "";
            bottom = "";
            for(int i = 0; i < numParticles; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_pressures[i]}\t|";
            }
            Debug.Log($"Pressures:\n{top}\n{bottom}");
        }

        // density_buffer.GetData(density);
        // for(int i = 0; i < density.Length; ++i)
        //     if(density[i] > max_density)
        //         max_density = density[i];
        // Debug.Log("max den = " + max_density);
        // compute_shader.SetFloat("max_density_multiplier", 1 / max_density);
    
        compute_shader.Dispatch(compute_force_kernel, thread_group_size, 1, 1);
        if (verbose_force) {
            float3[] temp_forces = new float3[numParticles];
            force_buffer.GetData(temp_forces);
            top = "";
            bottom = "";
            for(int i = 0; i < numParticles; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_forces[i]}\t|";
            }
            Debug.Log($"Forces:\n{top}\n{bottom}");
        }

        compute_shader.Dispatch(integrate_kernel, thread_group_size, 1, 1);
        if (verbose_integrate) {
            float[] temp_bounds = new float[n_bound];
            bound_buffer.GetData(temp_bounds);
            top = "";
            bottom = "";
            for(int i = 0; i < n_bound; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_bounds[i]}\t|";
            }
            Debug.Log($"Bounds:\n{top}\n{bottom}");
        }
        */
        
        if (show_gizmos) return;
        material.SetFloat(size_property, particleRenderSize);
        material.SetBuffer(particle_buffer_property, particle_buffer);
        Graphics.DrawMeshInstancedIndirect(particle_mesh, 0, material, new Bounds(Vector3.zero, new Vector3(1000f, 1000f, 1000f)), arg_buffer, castShadows: UnityEngine.Rendering.ShadowCastingMode.Off);
        
    }

    void OnDestroy()
    {
        free();
    }

    void free()
    {
        particle_buffer.Dispose();
        arg_buffer.Dispose();
        //neighbor_list_buffer.Dispose();
        //neighbor_tracker_buffer.Dispose();
        //hash_grid_buffer.Dispose();
        //hash_grid_tracker_buffer.Dispose();
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
        storedVelocitiesBuffer.Dispose();
        /* point_buffer.Dispose();
        triangle_buffer.Dispose();
        int_debug_buffer.Dispose();
        float_debug_buffer.Dispose(); */
    }
}