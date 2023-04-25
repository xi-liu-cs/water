using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;

public class fluid : MonoBehaviour
{
    [Header("particle")]
    public Mesh particle_mesh,
    fluid_mesh;
    public Material material;
    public float mass = 1f,
    viscosity_coefficient = 0.01f,
    particle_size = 8f,
    radius = 2f,
    grid_size = 8f, /* 4 * radius */
    gas_constant = 2000f,
    dt = 0.0008f,
    rest_density = 9f,
    damping = -1f;
    public Vector3 position_offset,
    velocity_initial = new Vector3(0, 500, 0);
    public float[] g = {0f, -9.81f * 2000f, 0f},
    bound = {-100, 200, -300, -100, -100, 100}; /* {-200, 200, -300, -12, -150, 150} */
    float radius2,
    radius3,
    radius4,
    radius5,
    mass2;
    int dimension2,
    dimension3;
    Vector3[] points;
    float[] noise_densities;
    int[] triangles;
    triangle[] march_triangles;

    [Header("fluid")]
    public int n_particle = 50000, /* max = 130000 */
    dimension = 100,
    max_particles_per_grid = 500;

    [Header("voxel")]
    public int n_point_per_axis = 50; /* z * n_point_per_axis * n_point_per_axis + y * n_point_per_axis + x, x ^ 3 + x ^ 2 + x = n_particle, (x = 50, n_particle = 130000) */
    public float isolevel = 8;

    public struct particle
    {
        public Vector3 position;
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

    public ComputeBuffer particle_buffer,
    arg_buffer,
    neighbor_list_buffer,
    neighbor_tracker_buffer,
    hash_grid_buffer,
    hash_grid_tracker_buffer,
    density_buffer,
    pressure_buffer,
    velocity_buffer,
    force_buffer,
    bound_buffer,
    point_buffer,
    noise_density_buffer,
    triangle_buffer,
    triangle_count_buffer,
    int_debug_buffer,
    float_debug_buffer;

    public ComputeShader compute_shader;
    int clear_hash_grid_kernel,
    compute_hash_grid_kernel,
    compute_neighbor_list_kernel,
    compute_density_pressure_kernel,
    compute_force_kernel,
    integrate_kernel,
    noise_density_kernel,
    march_kernel;

    public march_table table;

    public void Awake()
    {
        table = new march_table();
        radius2 = radius * radius;
        radius3 = radius2 * radius;
        radius4 = radius3 * radius;
        radius5 = radius4 * radius;
        mass2 = mass * mass;
        dimension2 = dimension * dimension;
        dimension3 = dimension * dimension2;
        thread_group_size = n_particle / dimension;

        /* gameObject.AddComponent<MeshFilter>();
        gameObject.AddComponent<MeshRenderer>(); */

        malloc_particle();
        find_kernel();
        compute_shader_init();
        compute_buffer_init();
        
        /* Debug.Log("group_size " + thread_group_size);
        compute_shader.Dispatch(noise_density_kernel, thread_group_size, thread_group_size, thread_group_size);
        Vector4[] cpu = new Vector4[n_particle];
        point_buffer.GetData(cpu);
        Debug.Log("cpu");
        for(int i = 0; i < cpu.Length; ++i)
            Debug.Log(cpu[i]); */
    }

    public void Update()
    {
        compute_shader.Dispatch(clear_hash_grid_kernel, dimension2, 1, 1);
        compute_shader.Dispatch(compute_hash_grid_kernel, thread_group_size, 1, 1);
        compute_shader.Dispatch(compute_neighbor_list_kernel, thread_group_size, 1, 1);
        compute_shader.Dispatch(compute_density_pressure_kernel, thread_group_size, 1, 1);
        compute_shader.Dispatch(compute_force_kernel, thread_group_size, 1, 1);
        compute_shader.Dispatch(integrate_kernel, thread_group_size, 1, 1);

        material.SetFloat(size_property, particle_size);
        material.SetBuffer(particle_buffer_property, particle_buffer);
        Graphics.DrawMeshInstancedIndirect(particle_mesh, 0, material, new Bounds(Vector3.zero, new Vector3(1000f, 1000f, 1000f)), arg_buffer, castShadows: UnityEngine.Rendering.ShadowCastingMode.Off);

        /* compute_shader.Dispatch(noise_density_kernel, thread_group_size, 1, 1); */

        /* int[] int_array_in_update_function = new int[n_debug];
        hash_grid_buffer.GetData(int_array_in_update_function);
        hash_grid_tracker_buffer.GetData(int_array_in_update_function);
        neighbor_list_buffer.GetData(int_array_in_update_function);
        neighbor_tracker_buffer.GetData(int_array_in_update_function);
        density_buffer.GetData(int_array_in_update_function);
        for(int i = 0; i < n_debug; ++i) Debug.Log(int_array_in_update_function[i]); */

        /* int march_kernel_n_thread = 4,
        march_kernel_group_size = (int)(Mathf.Pow(n_particle, 1.0f / 3.0f) / march_kernel_n_thread);
        triangle_buffer.SetCounterValue(0);
        compute_shader.Dispatch(march_kernel, march_kernel_group_size, march_kernel_group_size, march_kernel_group_size);
        
        ComputeBuffer.CopyCount(triangle_buffer, triangle_count_buffer, 0);
        int[] triangle_count_array = {0};
        triangle_count_buffer.GetData(triangle_count_array);
        int n_triangle = triangle_count_array[0];

        march_triangles = new triangle[n_triangle];
        triangle_buffer.GetData(march_triangles, 0, 0, n_triangle);
        Debug.Log("n_triangle " + n_triangle);
        for(int i = 0; i < n_triangle; ++i)
        {
            Debug.Log(march_triangles[i].vertex_a);
            Debug.Log(march_triangles[i].vertex_b);
            Debug.Log(march_triangles[i].vertex_c);
        }

        int_debug = new int[n_debug];
        int_debug_buffer.GetData(int_debug);
        for(int i = 0; i < n_debug; ++i)
            Debug.Log(String.Format("debug[{0}] = {1}", i, int_debug[i])); */
        
        /* debug_function(); */
        
        /* fluid_mesh = GetComponent<MeshFilter>().mesh;
        fluid_mesh.vertices = points;
        fluid_mesh.triangles = triangles;
        Graphics.DrawMesh(fluid_mesh, Vector3.zero, Quaternion.identity, material, 0); */
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
        int particle_per_dimension = Mathf.CeilToInt(Mathf.Pow(n_particle, 1f / 3f)),
        i = 0;
        while(i < n_particle)
        {
            for(int x = 0; x < particle_per_dimension; ++x)
                for(int y = 0; y < particle_per_dimension; ++y)
                    for(int z = 0; z < particle_per_dimension; ++z)
                    {
                        float r = UnityEngine.Random.Range(0f, 5f);
                        Vector3 pos = new Vector3(dimension - 1, dimension - 1, dimension - 1) - new Vector3(x / 2f, y / 2f, z / 2f)  - new Vector3(UnityEngine.Random.Range(0f, 0.01f), UnityEngine.Random.Range(0f, 0.01f), UnityEngine.Random.Range(0f, 0.01f)); /* Vector3 pos = new Vector3(1.5f * r * x, r * y, r * z) + position_offset; */
                        particles[i] = new particle{position = pos};
                        density[i] = -1;
                        pressure[i] = 0;
                        force[i] = Vector3.zero;
                        velocity[i] = velocity_initial;
                        if(++i == n_particle) return;
                    }
        }
    }

    void find_kernel()
    {
        clear_hash_grid_kernel = compute_shader.FindKernel("clear_hash_grid");
        compute_hash_grid_kernel = compute_shader.FindKernel("compute_hash_grid");
        compute_neighbor_list_kernel = compute_shader.FindKernel("compute_neighbor_list");
        compute_density_pressure_kernel = compute_shader.FindKernel("compute_density_pressure");
        compute_force_kernel = compute_shader.FindKernel("compute_force");
        integrate_kernel = compute_shader.FindKernel("integrate");
        noise_density_kernel = compute_shader.FindKernel("noise_density");
        march_kernel = compute_shader.FindKernel("march");
    }

    void compute_shader_init()
    {
        compute_shader.SetInt("n", 8);
        compute_shader.SetFloat("grid_size", grid_size);
        compute_shader.SetInt("dimension", dimension);
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
        compute_shader.SetFloat("dt", dt);
        compute_shader.SetFloats("g", g);
        compute_shader.SetFloats("epsilon", Mathf.Epsilon);
        compute_shader.SetFloats("pi", Mathf.PI);
        compute_shader.SetVector("time", Shader.GetGlobalVector("_Time"));
        compute_shader.SetInt("n_point_per_axis", n_point_per_axis);
        compute_shader.SetFloat("isolevel", isolevel);
    }

    public unsafe void compute_buffer_init()
    {
        uint[] arg = {particle_mesh.GetIndexCount(0), (uint)n_particle, particle_mesh.GetIndexStart(0), particle_mesh.GetBaseVertex(0), 0};
        particle_buffer = new ComputeBuffer(n_particle, sizeof(particle));
        particle_buffer.SetData(particles);
        arg_buffer = new ComputeBuffer(1, arg.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        arg_buffer.SetData(arg);
        neighbor_list = new int[n_particle * max_particles_per_grid * n];
        neighbor_tracker = new int[n_particle];
        hash_grid = new uint[dimension3 * max_particles_per_grid];
        hash_grid_tracker = new uint[dimension3];
        /* march_triangles = new triangle[n_particle]; */
        
        neighbor_list_buffer = new ComputeBuffer(n_particle * max_particles_per_grid * n, sizeof(int));
        neighbor_list_buffer.SetData(neighbor_list);
        neighbor_tracker_buffer = new ComputeBuffer(n_particle, sizeof(int));
        neighbor_tracker_buffer.SetData(neighbor_tracker);
        hash_grid_buffer = new ComputeBuffer(dimension3 * max_particles_per_grid, sizeof(uint));
        hash_grid_buffer.SetData(hash_grid);
        hash_grid_tracker_buffer = new ComputeBuffer(dimension3, sizeof(uint));
        hash_grid_tracker_buffer.SetData(hash_grid_tracker);
        density_buffer = new ComputeBuffer(n_particle, sizeof(float));
        density_buffer.SetData(density);
        pressure_buffer = new ComputeBuffer(n_particle, sizeof(float));
        pressure_buffer.SetData(pressure);
        velocity_buffer = new ComputeBuffer(n_particle, 3 * sizeof(float));
        velocity_buffer.SetData(velocity);
        force_buffer = new ComputeBuffer(n_particle, 3 * sizeof(float));
        force_buffer.SetData(force);
        bound_buffer = new ComputeBuffer(n_bound, sizeof(float));
        bound_buffer.SetData(bound);
        /* point_buffer = new ComputeBuffer(n_particle, 3 * sizeof(float));
        noise_density_buffer = new ComputeBuffer(n_particle, sizeof(float));
        triangle_buffer = new ComputeBuffer(n_particle, sizeof(triangle), ComputeBufferType.Append);
        triangle_count_buffer = new ComputeBuffer(1, sizeof(int), ComputeBufferType.Raw);
        int_debug_buffer = new ComputeBuffer(n_debug, sizeof(int));
        float_debug_buffer = new ComputeBuffer(n_debug, sizeof(float)); */

        compute_shader.SetBuffer(clear_hash_grid_kernel, "hash_grid_tracker", hash_grid_tracker_buffer);

        compute_shader.SetBuffer(compute_hash_grid_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_hash_grid_kernel, "hash_grid", hash_grid_buffer);
        compute_shader.SetBuffer(compute_hash_grid_kernel, "hash_grid_tracker", hash_grid_tracker_buffer);

        compute_shader.SetBuffer(compute_neighbor_list_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "hash_grid", hash_grid_buffer);
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "hash_grid_tracker", hash_grid_tracker_buffer);
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "neighbor_list", neighbor_list_buffer);
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "neighbor_tracker", neighbor_tracker_buffer);
        /* compute_shader.SetBuffer(compute_neighbor_list_kernel, "int_debug", int_debug_buffer); */
        
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

        /* compute_shader.SetBuffer(noise_density_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(noise_density_kernel, "points", point_buffer);
        compute_shader.SetBuffer(noise_density_kernel, "noise_densities", noise_density_buffer);

        compute_shader.SetBuffer(march_kernel, "points", point_buffer);
        compute_shader.SetBuffer(march_kernel, "noise_densities", noise_density_buffer);
        compute_shader.SetBuffer(march_kernel, "triangles", triangle_buffer);
        compute_shader.SetBuffer(march_kernel, "int_debug", int_debug_buffer);
        compute_shader.SetBuffer(march_kernel, "float_debug", float_debug_buffer); */
    }

    void OnDestroy()
    {
        free();
    }

    void free()
    {
        particle_buffer.Dispose();
        arg_buffer.Dispose();
        neighbor_list_buffer.Dispose();
        neighbor_tracker_buffer.Dispose();
        hash_grid_buffer.Dispose();
        hash_grid_tracker_buffer.Dispose();
        density_buffer.Dispose();
        pressure_buffer.Dispose();
        velocity_buffer.Dispose();
        force_buffer.Dispose();
        point_buffer.Dispose();
        triangle_buffer.Dispose();
    }
}