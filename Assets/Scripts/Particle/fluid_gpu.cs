using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;
using Unity.Mathematics;
using Random = UnityEngine.Random;

public class fluid_gpu : MonoBehaviour
{
    [Header("particle")]
    public Mesh particle_mesh,
    fluid_mesh;
    public Material material;
    public float mass = 1f;
    public float viscosity_coefficient = 0.01f;
    public float particle_size = 8f;
    public float radius = 8f; /* h, smoothing length */
    public float grid_size = 8f; /* 4 * radius */
    public float gas_constant = 2000f;
    public float dt = 0.0008f;
    public float rest_density = 9f;
    public float damping = -1f;
    
    public Vector3 position_offset = new Vector3(-140, -230, -60);
    public Vector3 velocity_initial = new Vector3(0, 10, 0);
    public float[] g = {0f, -9.81f, 0f};
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
    public int n_particle = 700000, /* max = 1000000 */
    max_particles_per_grid = 12,
    intial_particles_per_grid = 5;
    public Vector3Int dimension;
    int dimension2,
    dimension3;
    public int[] dimension_array;

    [Header("voxel")]
    public int n_point_per_axis = 50; /* z * n_point_per_axis * n_point_per_axis + y * n_point_per_axis + x, x ^ 3 + x ^ 2 + x = n_particle, (x = 50, n_particle = 130000) */
    public float isolevel = 8;

    public struct particle
    {
        public float3 position;
        public int3 grid_indices;
        public int projected_index;
        public int offset;
        public int index;
    }

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
    ceil_grid_size,
    n_debug = 64;
    particle[] particles;
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
    private ComputeBuffer kernel_sums_buffer;

    private ComputeBuffer gridBuffer;
    private ComputeBuffer gridOffsetBuffer;
    private ComputeBuffer gridSumsBuffer1;
    private ComputeBuffer gridSumsBuffer2;
    private ComputeBuffer rearrangedParticlesBuffer;

    public ComputeShader compute_shader;

    [Header("== KERNELS ==")]
    int clear_hash_grid_kernel;
    int compute_hash_grid_kernel;
    int compute_neighbor_list_kernel;
    int compute_sums_kernel;
    int compute_density_pressure_kernel;
    int compute_force_kernel;
    int integrate_kernel;
    int compute_density_kernel;
    int march_kernel;
    
    int generate_grid_kernel;
    int generate_particles_kernel;
    int clear_grid_kernel;
    int update_grid_kernel;
    int prefix_sum_kernel;
    int sum_blocks_kernel;
    int add_sums_kernel;
    int rearrange_particles_kernel;

    public march_table table;

    public float max_density = 0;

    public bool verbose_hash_grid = false;
    public bool verbose_neighbors = false;
    public bool verbose_density_pressure = false;
    public bool verbose_force = false;
    public bool verbose_integrate = false;
    public bool verbose_velocity = false;
    public bool verbose_prefix_sums = false;

    int _BLOCK_SIZE = 512;
    int numBlocks;

    void OnDrawGizmos() {
        Gizmos.color = Color.white;
        Gizmos.DrawWireCube(Vector3.zero, new Vector3(
            Mathf.Abs((float)bound[1] - (float)bound[0]),
            Mathf.Abs((float)bound[3] - (float)bound[2]),
            Mathf.Abs((float)bound[5] - (float)bound[4])
        ));
    }

    public void Awake()
    {
        dimension = new Vector3Int((int)((bound[1] - bound[0]) / grid_size), (int)((bound[3] - bound[2]) / grid_size), (int)((bound[5] - bound[4]) / grid_size));
        dimension_array = new int[]{dimension.x, dimension.y, dimension.z};
        dimension2 = dimension.x * dimension.y;
        dimension3 = dimension2 * dimension.z;
        grid_size_over_2 = (int)grid_size / 2;
        ceil_grid_size = (int)Math.Ceiling(grid_size);
        table = new march_table();
        radius2 = radius * radius;
        radius3 = radius2 * radius;
        radius4 = radius3 * radius;
        radius5 = radius4 * radius;
        mass2 = mass * mass;
        thread_group_size = n_particle / dimension.x;

        numBlocks = Mathf.CeilToInt((float)dimension3 / (float)_BLOCK_SIZE);

        //malloc_particle();
        find_kernel();
        compute_shader_init();
        compute_buffer_init();

        compute_shader.Dispatch(generate_grid_kernel, numBlocks,1,1);
        compute_shader.Dispatch(generate_particles_kernel, Mathf.CeilToInt((float)n_particle / (float)_BLOCK_SIZE),1,1);
        /*
        particle[] temp_particles = new particle[n_particle];
        particle_buffer.GetData(temp_particles);
        string top = "";
        string bottom = "";
        for(int i = 0; i < n_particle; i++) {
            top += $"{i}\t|";
            bottom += $"{temp_particles[i].index}\t|";
        }
        Debug.Log($"Particle Indexes:\n{top}\n{bottom}");
        */

        /* Debug.Log("group_size " + thread_group_size);
        compute_shader.Dispatch(noise_density_kernel, thread_group_size, thread_group_size, thread_group_size);
        Vector4[] cpu = new Vector4[n_particle];
        point_buffer.GetData(cpu);
        Debug.Log("cpu");
        for(int i = 0; i < cpu.Length; ++i)
            Debug.Log(cpu[i]); */
    }


    void debug_function()
    {
        points = new Vector3[n_particle];
        point_buffer.GetData(points);
        Debug.Log("point");
        Debug.Log(points[0]);

        noise_densities = new float[n_particle];
        noise_density_buffer.GetData(noise_densities);
        Debug.Log("noise density");
        Debug.Log(noise_densities[0]);

        triangles = triangulate_field(points);
        Debug.Log("triangles");
        for(int i = 0; i < 10; ++i)
            Debug.Log(triangles[i]);
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
        int[] a = new int[n_particle];
        int idx = 0;
        for(int i = 0; i + 1 < n_point_per_axis; ++i)
        {
            for(int j = 0; j + 1 < n_point_per_axis; ++j)
            {
                for(int k = 0; k + 1 < n_point_per_axis; ++k)
                {
                    float x = i, y = j, z = k;
                    if(index_from_coordinate(i + 1, j + 1, k + 1) >= n_particle)
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

    void malloc_particle()
    {
        particles = new particle[n_particle];
        density = new float[n_particle];
        pressure = new float[n_particle];
        velocity = new Vector3[n_particle];
        force = new Vector3[n_particle];
        int i = 0;
        while(i < n_particle)
        {
            for(int x = 0; x < dimension.x; ++x)
                for(int y = 0; y < dimension.y; ++y)
                    for(int z = 0; z < dimension.z; ++z)
                    {
                        Vector3 grid_pos = new Vector3(
                            bound[0] + x * grid_size, 
                            bound[2] + y * grid_size, 
                            bound[4] + z * grid_size
                        );
                        for(int a = 0; a < intial_particles_per_grid; ++a)
                        {
                            //Vector3 pos = grid_pos - new Vector3(grid_size_over_2, grid_size_over_2, grid_size_over_2) - new Vector3(UnityEngine.Random.Range(0f, grid_size_over_2), UnityEngine.Random.Range(0f, grid_size_over_2), UnityEngine.Random.Range(0f, grid_size_over_2));
                            Vector3 pos = new Vector3(
                                Random.Range(bound[0],bound[1]),
                                Random.Range(bound[2],bound[3]),
                                Random.Range(bound[4],bound[5])
                            );
                            Vector3Int grid_indices = GetGridXYZIndices(pos);
                            int projected_index = GetProjectedGridIndexFromXYZ(grid_indices);
                            particles[i] = new particle{
                                position = new(pos.x, pos.y, pos.z), 
                                grid_indices = new(grid_indices.x, grid_indices.y, grid_indices.z),
                                projected_index = projected_index,
                                index = i
                            };
                            density[i] = -1;
                            pressure[i] = 0;
                            force[i] = Vector3.zero;
                            velocity[i] = velocity_initial;
                            if(++i == n_particle) return;
                        }
                    }
        }
    }
    public Vector3Int GetGridXYZIndices(Vector3 position) {
        return new Vector3Int(
            Mathf.FloorToInt((position.x - bound[0])/grid_size),
            Mathf.FloorToInt((position.y - bound[2])/grid_size),
            Mathf.FloorToInt((position.z - bound[4])/grid_size)
        );
    }
    public int GetProjectedGridIndexFromXYZ(int x, int y, int z) {
        return x + (dimension_array[1] * y) + (dimension_array[1] * dimension_array[2] * z);
    }
    public int GetProjectedGridIndexFromXYZ(Vector3Int xyz) {
        return xyz.x + (dimension_array[1] * xyz.y) + (dimension_array[1] * dimension_array[2] * xyz.z);
    }

    void find_kernel() {
        generate_grid_kernel = compute_shader.FindKernel("GenerateGrid");
        generate_particles_kernel = compute_shader.FindKernel("GenerateParticles");
        clear_grid_kernel = compute_shader.FindKernel("ClearGrid");
        update_grid_kernel = compute_shader.FindKernel("UpdateGridCellCounts");
        prefix_sum_kernel = compute_shader.FindKernel("PrefixSum");
        sum_blocks_kernel = compute_shader.FindKernel("SumBlocks");
        add_sums_kernel = compute_shader.FindKernel("AddSums");
        rearrange_particles_kernel = compute_shader.FindKernel("RearrangeParticles");
        compute_density_pressure_kernel = compute_shader.FindKernel("ComputeDensityPressure");
        compute_force_kernel = compute_shader.FindKernel("ComputeForce");
        integrate_kernel = compute_shader.FindKernel("Integrate");

        /*
        // clear_hash_grid_kernel = compute_shader.FindKernel("clear_hash_grid");
        // compute_hash_grid_kernel =  compute_shader.FindKernel("compute_hash_grid");
        compute_neighbor_list_kernel = compute_shader.FindKernel("compute_neighbor_list");
        compute_sums_kernel = compute_shader.FindKernel("compute_sums");
        compute_density_pressure_kernel = compute_shader.FindKernel("compute_density_pressure");
        compute_force_kernel = compute_shader.FindKernel("compute_force");
        integrate_kernel = compute_shader.FindKernel("integrate");
        */
    }

    void compute_shader_init()
    {
        compute_shader.SetInt("n", 8);
        compute_shader.SetInt("n_particle", n_particle);
        compute_shader.SetFloat("grid_size", grid_size);
        compute_shader.SetInt("max_particles_per_grid", max_particles_per_grid);
        compute_shader.SetFloat("radius", radius);
        compute_shader.SetFloat("radius2", radius2);
        compute_shader.SetFloat("radius3", radius3);
        compute_shader.SetFloat("radius4", radius4);
        compute_shader.SetFloat("radius5", radius5);
        compute_shader.SetFloat("mass", mass);
        compute_shader.SetFloat("mass2", mass2);
        compute_shader.SetFloat("gas_constant", gas_constant);
        compute_shader.SetFloat("rest_density", rest_density);
        compute_shader.SetFloat("viscosity_coefficient", viscosity_coefficient);
        compute_shader.SetFloat("damping", damping);
        //compute_shader.SetFloat("dt", dt);
        compute_shader.SetFloats("g", g);
        compute_shader.SetFloat("epsilon", Mathf.Epsilon);
        compute_shader.SetFloat("e", Mathf.Exp(1));
        compute_shader.SetFloat("pi", Mathf.PI);
        compute_shader.SetVector("time", Shader.GetGlobalVector("_Time"));
        compute_shader.SetInt("n_point_per_axis", n_point_per_axis);
        compute_shader.SetFloat("isolevel", isolevel);
        compute_shader.SetInts("dimension", dimension_array);

        compute_shader.SetInt("total_number_of_cells", dimension3);
        compute_shader.SetInt("numBlocks", numBlocks);
        compute_shader.SetInt("randomSeed", Random.Range(0,int.MaxValue));
    }

    public unsafe void compute_buffer_init()
    {
        uint[] arg = {particle_mesh.GetIndexCount(0), (uint)n_particle, particle_mesh.GetIndexStart(0), particle_mesh.GetBaseVertex(0), 0};
        particle_buffer = new ComputeBuffer(n_particle, sizeof(float)*3 + sizeof(int)*6);
        //particle_buffer.SetData(particles);
        arg_buffer = new ComputeBuffer(1, arg.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        arg_buffer.SetData(arg);
        neighbor_list = new int[n_particle * max_particles_per_grid * n];
        neighbor_tracker = new int[n_particle];
        hash_grid = new uint[dimension3 * max_particles_per_grid];
        hash_grid_tracker = new uint[dimension3];
        /* march_triangles = new triangle[n_particle]; */
        
        /*
        neighbor_list_buffer = new ComputeBuffer(n_particle * max_particles_per_grid * n, sizeof(int));
        neighbor_list_buffer.SetData(neighbor_list);
        neighbor_tracker_buffer = new ComputeBuffer(n_particle, sizeof(int));
        neighbor_tracker_buffer.SetData(neighbor_tracker);
        //hash_grid_buffer = new ComputeBuffer(dimension3 * max_particles_per_grid, sizeof(uint));
        //hash_grid_buffer.SetData(hash_grid);
        hash_grid_tracker_buffer = new ComputeBuffer(dimension3, sizeof(uint));
        hash_grid_tracker_buffer.SetData(hash_grid_tracker);
        */
        density_buffer = new ComputeBuffer(n_particle, sizeof(float));
        //density_buffer.SetData(density);
        pressure_buffer = new ComputeBuffer(n_particle, sizeof(float));
        //pressure_buffer.SetData(pressure);
        velocity_buffer = new ComputeBuffer(n_particle, 3 * sizeof(float));
        //velocity_buffer.SetData(velocity);
        force_buffer = new ComputeBuffer(n_particle, 3 * sizeof(float));
        //force_buffer.SetData(force);

        bound_buffer = new ComputeBuffer(n_bound, sizeof(float));
        bound_buffer.SetData(bound);

        gridBuffer = new ComputeBuffer(dimension3, sizeof(int));
        gridOffsetBuffer = new ComputeBuffer(dimension3, sizeof(int));
        gridSumsBuffer1 = new ComputeBuffer(numBlocks, sizeof(int));
        gridSumsBuffer2 = new ComputeBuffer(numBlocks, sizeof(int));
        rearrangedParticlesBuffer = new ComputeBuffer(n_particle, sizeof(float)*3 + sizeof(int)*6);

        compute_shader.SetBuffer(generate_grid_kernel, "grid", gridBuffer);

        compute_shader.SetBuffer(generate_particles_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(generate_particles_kernel, "density", density_buffer);
        compute_shader.SetBuffer(generate_particles_kernel, "pressure", pressure_buffer);
        compute_shader.SetBuffer(generate_particles_kernel, "force", force_buffer);
        compute_shader.SetBuffer(generate_particles_kernel, "velocity", velocity_buffer);
        compute_shader.SetBuffer(generate_particles_kernel, "bound", bound_buffer);

        // compute_shader.SetBuffer(clear_hash_grid_kernel, "hash_grid_tracker", hash_grid_tracker_buffer);
        compute_shader.SetBuffer(clear_grid_kernel, "grid", gridBuffer);
        
        // compute_shader.SetBuffer(compute_hash_grid_kernel, "particles", particle_buffer);
        // compute_shader.SetBuffer(compute_hash_grid_kern1el, "hash_grid", hash_grid_buffer);
        // compute_shader.SetBuffer(compute_hash_grid_kernel, "hash_grid_tracker", hash_grid_tracker_buffer);
        // compute_shader.SetBuffer(compute_hash_grid_kernel, "bound", bound_buffer);
        compute_shader.SetBuffer(update_grid_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(update_grid_kernel, "grid", gridBuffer);
        compute_shader.SetBuffer(update_grid_kernel, "bound", bound_buffer);

        compute_shader.SetBuffer(prefix_sum_kernel, "gridOffsetBufferIn", gridBuffer);
        compute_shader.SetBuffer(prefix_sum_kernel, "gridOffsetBuffer", gridOffsetBuffer);
        compute_shader.SetBuffer(prefix_sum_kernel, "gridSumsBuffer", gridSumsBuffer2);

        compute_shader.SetBuffer(add_sums_kernel, "gridOffsetBuffer", gridOffsetBuffer);

        compute_shader.SetBuffer(rearrange_particles_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(rearrange_particles_kernel, "rearrangedParticles", rearrangedParticlesBuffer);
        compute_shader.SetBuffer(rearrange_particles_kernel, "gridOffsetBuffer", gridOffsetBuffer);

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
        compute_shader.SetBuffer(integrate_kernel, "bound", bound_buffer);

        /*
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "hash_grid", hash_grid_buffer);
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "hash_grid_tracker", hash_grid_tracker_buffer);
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "neighbor_list", neighbor_list_buffer);
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "neighbor_tracker", neighbor_tracker_buffer);
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "bound", bound_buffer);

        compute_shader.SetBuffer(compute_density_pressure_kernel, "neighbor_list", neighbor_list_buffer);
        compute_shader.SetBuffer(compute_density_pressure_kernel, "neighbor_tracker", neighbor_tracker_buffer);
        compute_shader.SetBuffer(compute_density_pressure_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_density_pressure_kernel, "density", density_buffer);
        compute_shader.SetBuffer(compute_density_pressure_kernel, "pressure", pressure_buffer);

        compute_shader.SetBuffer(compute_force_kernel, "neighbor_list", neighbor_list_buffer);
        compute_shader.SetBuffer(compute_force_kernel, "neighbor_tracker", neighbor_tracker_buffer);
        compute_shader.SetBuffer(compute_force_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_force_kernel, "density", density_buffer);
        compute_shader.SetBuffer(compute_force_kernel, "pressure", pressure_buffer);
        compute_shader.SetBuffer(compute_force_kernel, "velocity", velocity_buffer);
        compute_shader.SetBuffer(compute_force_kernel, "force", force_buffer);

        compute_shader.SetBuffer(integrate_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(integrate_kernel, "velocity", velocity_buffer);
        compute_shader.SetBuffer(integrate_kernel, "force", force_buffer);
        compute_shader.SetBuffer(integrate_kernel, "bound", bound_buffer);
        */
    }

    public void Update()
    {

        compute_shader.SetFloat("dt", Time.deltaTime);

        string top = "";
        string bottom = "";

        //compute_shader.Dispatch(clear_hash_grid_kernel, dimension2, 1, 1);
        compute_shader.Dispatch(clear_grid_kernel, numBlocks,1,1);
        //compute_shader.Dispatch(compute_hash_grid_kernel, thread_group_size, 1, 1);
        compute_shader.Dispatch(update_grid_kernel, Mathf.CeilToInt((float)n_particle / (float)_BLOCK_SIZE),1,1);

        if (verbose_hash_grid) {
            uint[] temp_grid = new uint[dimension3];
            gridBuffer.GetData(temp_grid);
            for(int i = 0; i < dimension3; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_grid[i]}\t|";
            }
            Debug.Log($"Grid\n{top}\n{bottom}");
            particle[] temp_particles = new particle[n_particle];
            particle_buffer.GetData(temp_particles);
            top = "";
            bottom = "";
            for(int i = 0; i < n_particle; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_particles[i].position}\t|";
            }
            Debug.Log($"Particles World Positions:\n{top}\n{bottom}");
            top = "";
            bottom = "";
            for(int i = 0; i < n_particle; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_particles[i].grid_indices}\t|";
            }
            Debug.Log($"Particles Grid Indices:\n{top}\n{bottom}");
            top = "";
            bottom = "";
            for(int i = 0; i < n_particle; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_particles[i].projected_index}\t|";
            }
            Debug.Log($"Particles Projected Indices:\n{top}\n{bottom}");
            top = "";
            bottom = "";
            for(int i = 0; i < n_particle; i++) {
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
            int[] temp_offset_grid = new int[dimension3];
            gridOffsetBuffer.GetData(temp_offset_grid);
            top = "";
            bottom="";
            for(int k = 0; k < dimension3; k++) {
                top += $"{k}\t|";
                bottom += $"{temp_offset_grid[k]}\t|";
            }
            Debug.Log("PREFIX SUMS:\n"+top+"\n"+bottom);
        }

        compute_shader.Dispatch(rearrange_particles_kernel, Mathf.CeilToInt((float)n_particle / (float)_BLOCK_SIZE), 1, 1);

        compute_shader.Dispatch(compute_density_pressure_kernel, Mathf.CeilToInt((float)n_particle / (float)_BLOCK_SIZE), 1, 1);
        if (verbose_density_pressure) {
            float[] temp_densities = new float[n_particle];
            float[] temp_pressures = new float[n_particle];
            density_buffer.GetData(temp_densities);
            top = "";
            bottom = "";
            for(int i = 0; i < n_particle; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_densities[i]}\t|";
            }
            Debug.Log($"Densities:\n{top}\n{bottom}");
            pressure_buffer.GetData(temp_pressures);
            top = "";
            bottom = "";
            for(int i = 0; i < n_particle; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_pressures[i]}\t|";
            }
            Debug.Log($"Pressures:\n{top}\n{bottom}");
        }

        compute_shader.Dispatch(compute_force_kernel, Mathf.CeilToInt((float)n_particle / (float)_BLOCK_SIZE), 1, 1);
        if (verbose_force) {
            float3[] temp_forces = new float3[n_particle];
            force_buffer.GetData(temp_forces);
            top = "";
            bottom = "";
            for(int i = 0; i < n_particle; i++) {
                if (
                    Mathf.Abs(temp_forces[i][0]) < Mathf.Infinity
                    && Mathf.Abs(temp_forces[i][1]) < Mathf.Infinity
                    && Mathf.Abs(temp_forces[i][2]) < Mathf.Infinity
                ) continue;
                top += $"{i}\t|";
                bottom += $"{temp_forces[i]}\t|";
            }
            Debug.Log($"Forces:\n{top}\n{bottom}");
        }

        compute_shader.Dispatch(integrate_kernel, Mathf.CeilToInt((float)n_particle / (float)_BLOCK_SIZE), 1, 1);
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
        if (verbose_velocity) {
            float3[] temp_velocities = new float3[n_particle];
            velocity_buffer.GetData(temp_velocities);
            top = "";
            bottom = "";
            for(int i = 0; i < n_particle; i++) {
                if (
                    Mathf.Abs(temp_velocities[i][0]) < Mathf.Infinity
                    && Mathf.Abs(temp_velocities[i][1]) < Mathf.Infinity
                    && Mathf.Abs(temp_velocities[i][2]) < Mathf.Infinity
                ) continue;
                top += $"{i}\t|";
                bottom += $"{temp_velocities[i]}\t|";
            }
            Debug.Log($"Velocities:\n{top}\n{bottom}");
        }

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
            particle[] temp_particles = new particle[n_particle];
            particle_buffer.GetData(temp_particles);
            top = "";
            bottom = "";
            for(int i = 0; i < n_particle; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_particles[i].projected_index}\t|";
            }
            Debug.Log($"Particles Projected Indices:\n{top}\n{bottom}");
            top = "";
            bottom = "";
            for(int i = 0; i < n_particle; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_particles[i].offset}\t|";
            }
            Debug.Log($"Particles Offsets:\n{top}\n{bottom}");
        }
        */

        /*
        compute_shader.Dispatch(compute_neighbor_list_kernel, thread_group_size, 1, 1);
        
        if (verbose_neighbors) {
            int[] temp_neighbor_list = new int[n_particle * max_particles_per_grid * n];
            int[] temp_neighbor_tracker = new int[n_particle];
            neighbor_list_buffer.GetData(temp_neighbor_list);
            neighbor_tracker_buffer.GetData(temp_neighbor_tracker);
            top = "";
            bottom = "";
            for(int i = 0; i < n_particle * max_particles_per_grid * n; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_neighbor_list[i]}\t|";
            }
            Debug.Log($"Neighbor list\n{top}\n{bottom}");
            top = "";
            bottom = "";
            for(int i = 0; i < n_particle; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_neighbor_tracker[i]}\t|";
            }
            Debug.Log($"Neighbor Tracker\n{top}\n{bottom}");
        }
        compute_shader.Dispatch(compute_density_pressure_kernel, thread_group_size, 1, 1);
        if (verbose_density_pressure) {
            float[] temp_densities = new float[n_particle];
            float[] temp_pressures = new float[n_particle];
            density_buffer.GetData(temp_densities);
            top = "";
            bottom = "";
            for(int i = 0; i < n_particle; i++) {
                top += $"{i}\t|";
                bottom += $"{temp_densities[i]}\t|";
            }
            Debug.Log($"Densities:\n{top}\n{bottom}");
            pressure_buffer.GetData(temp_pressures);
            top = "";
            bottom = "";
            for(int i = 0; i < n_particle; i++) {
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
            float3[] temp_forces = new float3[n_particle];
            force_buffer.GetData(temp_forces);
            top = "";
            bottom = "";
            for(int i = 0; i < n_particle; i++) {
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

        material.SetFloat(size_property, particle_size);
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
        bound_buffer.Dispose();
        force_buffer.Dispose();
        
        gridBuffer.Dispose();
        gridOffsetBuffer.Dispose();
        gridSumsBuffer1.Dispose();
        gridSumsBuffer2.Dispose();
        rearrangedParticlesBuffer.Dispose();
        /* point_buffer.Dispose();
        triangle_buffer.Dispose();
        int_debug_buffer.Dispose();
        float_debug_buffer.Dispose(); */
    }
}