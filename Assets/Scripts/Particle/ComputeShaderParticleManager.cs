using System.Runtime.InteropServices;
using UnityEngine;
using Random = UnityEngine.Random;

public class ComputeShaderParticleManager : MonoBehaviour
{
    [Header("particle")]
    public Mesh particle_mesh,
    fluid_mesh;
    public Material material;
    public float mass = 4f,
    viscosity_coefficient = 2.5f,
    particle_size = 40f,
    radius = 1f,
    grid_size = 4f, /* 4 * radius */
    gas_constant = 2000f,
    dt = 0.0008f,
    rest_density = 9f,
    damping = -1f;
    public Vector3 position_offset,
    velocity_initial = new Vector3(0, 500, 0);
    public float[] g = {0f, -9.81f * 2000f, 0f},
    bound = {0, 0, 0, -300, 0, 0};
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
        public Vector4 color;
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

    int size_property = Shader.PropertyToID("_size"),
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

    private void Awake()
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

    #region Initialisation
    
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
                        Vector3 pos = new Vector3(1.5f * r * x, r * y, r * z) + position_offset; /* Vector3 pos = new Vector3(dimension - 1, dimension - 1, dimension - 1) - new Vector3(x / 2f, y / 2f, z / 2f)  - new Vector3(UnityEngine.Random.Range(0f, 0.01f), UnityEngine.Random.Range(0f, 0.01f), UnityEngine.Random.Range(0f, 0.01f)); */
                        particles[i] = new particle
                        {
                            position = pos,
                            color = new Vector4(0.3f, 0.5f, 1f, 0.5f)
                        };
                        density[i] = -1;
                        pressure[i] = 0;
                        force[i] = Vector3.zero;
                        velocity[i] = velocity_initial;
                        if(++i == n_particle) return;
                    }
        }
    }

    private void find_kernel()
    {
        clear_hash_grid_kernel = compute_shader.FindKernel("clear_hash_grid");
        compute_hash_grid_kernel = compute_shader.FindKernel("compute_hash_grid");
        compute_neighbor_list_kernel = compute_shader.FindKernel("compute_neighbor_list");
        compute_density_pressure_kernel = compute_shader.FindKernel("compute_density_pressure");
        compute_force_kernel = compute_shader.FindKernel("compute_force");
        integrate_kernel = compute_shader.FindKernel("integrate");
    }

    private void compute_shader_init()
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

    unsafe void compute_buffer_init()
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
        march_triangles = new triangle[n_particle];
        
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
        point_buffer = new ComputeBuffer(n_particle, 3 * sizeof(float));
        noise_density_buffer = new ComputeBuffer(n_particle, sizeof(float));
        triangle_buffer = new ComputeBuffer(n_particle, sizeof(triangle), ComputeBufferType.Append);
        triangle_count_buffer = new ComputeBuffer(1, sizeof(int), ComputeBufferType.Raw);
        int_debug_buffer = new ComputeBuffer(n_debug, sizeof(int));
        float_debug_buffer = new ComputeBuffer(n_debug, sizeof(float));

        compute_shader.SetBuffer(clear_hash_grid_kernel, "hash_grid_tracker", hash_grid_tracker_buffer);

        compute_shader.SetBuffer(compute_hash_grid_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_hash_grid_kernel, "hash_grid", hash_grid_buffer);
        compute_shader.SetBuffer(compute_hash_grid_kernel, "hash_grid_tracker", hash_grid_tracker_buffer);

        compute_shader.SetBuffer(compute_neighbor_list_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "hash_grid", hash_grid_buffer);
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "hash_grid_tracker", hash_grid_tracker_buffer);
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "neighbor_list", neighbor_list_buffer);
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "neighbor_tracker", neighbor_tracker_buffer);
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "int_debug", int_debug_buffer);
        
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

        compute_shader.SetBuffer(noise_density_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(noise_density_kernel, "points", point_buffer);
        compute_shader.SetBuffer(noise_density_kernel, "noise_densities", noise_density_buffer);

        compute_shader.SetBuffer(march_kernel, "points", point_buffer);
        compute_shader.SetBuffer(march_kernel, "noise_densities", noise_density_buffer);
        compute_shader.SetBuffer(march_kernel, "triangles", triangle_buffer);
        compute_shader.SetBuffer(march_kernel, "int_debug", int_debug_buffer);
        compute_shader.SetBuffer(march_kernel, "float_debug", float_debug_buffer);
    }
    
    #endregion

    void Update()
    {
        compute_shader.Dispatch(clear_hash_grid_kernel, dimension * dimension * dimension / 100, 1, 1);
        compute_shader.Dispatch(compute_hash_grid_kernel, n_particle / 100, 1, 1);
        compute_shader.Dispatch(compute_neighbor_list_kernel, n_particle / 100, 1, 1);
        compute_shader.Dispatch(compute_density_pressure_kernel, n_particle / 100, 1, 1);
        compute_shader.Dispatch(compute_force_kernel, n_particle / 100, 1, 1);
        compute_shader.Dispatch(integrate_kernel, n_particle / 100, 1, 1);
        
        material.SetFloat(size_property, particle_size);
        material.SetBuffer(particle_buffer_property, particle_buffer);
        Graphics.DrawMeshInstancedIndirect(particle_mesh, 0, material, new Bounds(Vector3.zero, new Vector3(100.0f, 100.0f, 100.0f)), arg_buffer, castShadows: UnityEngine.Rendering.ShadowCastingMode.Off);
        int[] nei = new int[1200];
        // hash_grid_buffer.GetData(nei); // 1
        // hash_grid_tracker_buffer.GetData(nei); // 1
        neighbor_list_buffer.GetData(nei); // 1
        // neighbor_tracker_buffer.GetData(nei); // 1
        // density_buffer.GetData(nei); // 1
        // pressure_buffer.GetData(nei); // 1
        // velocity_buffer.GetData(nei); // 1
        // force_buffer.GetData(nei); // 1
        for(int i = 0; i < 1000; ++i) if(nei[i] != 0) Debug.Log(nei[i]);
    }

    private void OnDestroy()
    {
        ReleaseBuffers();
    }

    private void ReleaseBuffers()
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
    }
}
