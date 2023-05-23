using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;

public class fluid2 : MonoBehaviour
{
    [Header("particle")]
    public Mesh particle_mesh,
    fluid_mesh;
    public Material material;
    public float mass = 1f,
    viscosity_coefficient = 0.1f,
    particle_size = 2f,
    radius = 2f, /* h, smoothing length */
    grid_size = 2f, /* 4 * radius */
    gas_constant = 2000f,
    dt = 0.0001f,
    rest_density = 1000f,
    damping = -0.5f,
    bulk_modulus = 500f;
    public Vector3 position_offset = new Vector3(0, 0, 0),
    velocity_initial = new Vector3(0, 0, 0);
    public float[] g = {0f, -9.81f, 0f},
    bound = {-50, 50, -50, 50, -50, 50}; /* {-160, 140, -250, -150, -100, 100} */
    [HideInInspector]
    public float radius2,
    radius3,
    radius4,
    radius5,
    radius8,
    mass2;
    Vector3[] points;
    float[] noise_densities;
    int[] triangles;
    triangle[] march_triangles;

    [Header("fluid")]
    public int n_particle = 100000, /* max = 1000000 */
    max_particles_per_grid = 12,
    intial_particles_per_grid = 5;
    public Vector3Int dimension;
    int dimension2,
    dimension3;
    public int[] dimension_array;

    [Header("voxel")]
    public int n_point_per_axis = 84; /* z * n_point_per_axis * n_point_per_axis + y * n_point_per_axis + x, x ^ 3 + x ^ 2 + x = n_particle, (x = 50, n_particle = 130000) */
    public float isolevel = 0;

    /* need to match compute shader file */
    int thread_group_size = 512,
    dispatch_size,
    cell_dispatch_size;

    public struct particle /* struct alignment */
    {
        public Vector3 position,
        velocity,
        acceleration;
        public float density,
        pressure;
        public int neighbor;
    };

    public struct triangle
    {
        public Vector3 vertex_a,
        vertex_b,
        vertex_c;
    }

    int size_property = Shader.PropertyToID("size"),
    particle_buffer_property = Shader.PropertyToID("particle_buffer");

    int n = 8,
    n_bound = 6,
    grid_size_over_2,
    ceil_grid_size,
    n_debug = 64,
    n_cell;
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
    bound_buffer,
    density_buffer,
    hash_grid_buffer,
    hash_grid_tracker_buffer;

    public ComputeShader compute_shader;
    int malloc_particle_kernel,
    clear_hash_grid_kernel,
    compute_hash_grid_kernel,
    compute_density_kernel,
    compute_acceleration_kernel,
    integrate_kernel;

    public march_table table;

    public void Awake()
    {
        n_point_per_axis = (int)Mathf.Pow(n_particle, 1.0f / 3.0f);
        dimension = new Vector3Int((int)((bound[1] - bound[0]) / grid_size), (int)((bound[3] - bound[2]) / grid_size), (int)((bound[5] - bound[4]) / grid_size));
        dimension_array = new int[]{dimension.x, dimension.y, dimension.z};
        dimension2 = dimension.x * dimension.y;
        dimension3 = dimension2 * dimension.z;
        grid_size_over_2 = (int)grid_size / 2;
        ceil_grid_size = (int)Math.Ceiling(grid_size);
        radius2 = radius * radius;
        radius3 = radius2 * radius;
        radius4 = radius3 * radius;
        radius5 = radius4 * radius;
        radius8 = radius2 * radius4;
        mass2 = mass * mass;
        n_cell = dimension.x + dimension.x * (dimension.y + dimension.y * dimension.z);
        dispatch_size = Mathf.CeilToInt((float)n_particle / (float)thread_group_size);
        cell_dispatch_size = Mathf.CeilToInt((float)(n_cell + 1) / (float)thread_group_size);

        find_kernel();
        compute_shader_init();
        compute_buffer_init();
        compute_shader.Dispatch(malloc_particle_kernel, thread_group_size, 1, 1);
    }

    public void Update()
    {
        compute_shader.Dispatch(clear_hash_grid_kernel, cell_dispatch_size, 1, 1);
        compute_shader.Dispatch(compute_hash_grid_kernel, dispatch_size, 1, 1);
        compute_shader.Dispatch(compute_density_kernel, dispatch_size, 1, 1);
        compute_shader.Dispatch(compute_acceleration_kernel, dispatch_size, 1, 1);
        compute_shader.Dispatch(integrate_kernel, dispatch_size, 1, 1);
        material.SetFloat(size_property, particle_size);
        material.SetBuffer(particle_buffer_property, particle_buffer);
        Graphics.DrawMeshInstancedIndirect(particle_mesh, 0, material, new Bounds(Vector3.zero, new Vector3(1000f, 1000f, 1000f)), arg_buffer, castShadows: UnityEngine.Rendering.ShadowCastingMode.Off);
        particle[] debug_buf = new particle[100000];
        particle_buffer.GetData(debug_buf);
        for(int i = 0; i < 100; ++i) Debug.Log("position = " + debug_buf[i].position
        + ", velocity = " + debug_buf[i].velocity
        + ", acceleration = " + debug_buf[i].acceleration
        + ", density = " + debug_buf[i].density
        + ", pressure = " + debug_buf[i].pressure
        + ", neighbor = " + debug_buf[i].neighbor);
        /* int[] neighbor_debug_buf = new int[100];
        neighbor_buffer.GetData(neighbor_debug_buf);
        Debug.Log("neighbor: ");
        for(int i = 0; i < 100; ++i) Debug.Log(debug_buf[neighbor_debug_buf[i]].position); */
    }

    void find_kernel()
    {
        malloc_particle_kernel = compute_shader.FindKernel("malloc_particle");
        clear_hash_grid_kernel = compute_shader.FindKernel("clear_hash_grid");
        compute_hash_grid_kernel = compute_shader.FindKernel("compute_hash_grid");
        compute_density_kernel = compute_shader.FindKernel("compute_density");
        compute_acceleration_kernel = compute_shader.FindKernel("compute_acceleration");
        integrate_kernel = compute_shader.FindKernel("integrate");
    }

    void compute_shader_init()
    {
        compute_shader.SetInt("n", 8);
        compute_shader.SetFloat("grid_size", grid_size);
        compute_shader.SetInt("max_particles_per_grid", max_particles_per_grid);
        compute_shader.SetFloat("radius", radius);
        compute_shader.SetFloat("radius2", radius2);
        compute_shader.SetFloat("radius3", radius3);
        compute_shader.SetFloat("radius4", radius4);
        compute_shader.SetFloat("radius5", radius5);
        compute_shader.SetFloat("radius8", radius8);
        compute_shader.SetFloat("particle_size", particle_size);
        compute_shader.SetFloat("particle_size2", 2 * particle_size);
        compute_shader.SetFloat("mass", mass);
        compute_shader.SetFloat("mass2", mass2);
        compute_shader.SetFloat("gas_constant", gas_constant);
        compute_shader.SetFloat("rest_density", rest_density);
        compute_shader.SetFloat("viscosity_coefficient", viscosity_coefficient);
        compute_shader.SetFloat("damping", damping);
        compute_shader.SetFloat("dt", dt);
        compute_shader.SetFloats("g", g);
        compute_shader.SetFloat("epsilon", Mathf.Epsilon);
        compute_shader.SetFloat("e", Mathf.Exp(1));
        compute_shader.SetFloat("pi", Mathf.PI);
        compute_shader.SetFloat("bulk_modulus", bulk_modulus);
        compute_shader.SetVector("time", Shader.GetGlobalVector("_Time"));
        compute_shader.SetInt("n_point_per_axis", n_point_per_axis);
        compute_shader.SetFloat("isolevel", isolevel);
        compute_shader.SetInts("dimension", dimension_array);
        compute_shader.SetInt("n_particle", n_particle);
        compute_shader.SetInt("seed", UnityEngine.Random.Range(0, int.MaxValue));
    }

    public unsafe void compute_buffer_init()
    {
        uint[] arg = {particle_mesh.GetIndexCount(0), (uint)n_particle, particle_mesh.GetIndexStart(0), particle_mesh.GetBaseVertex(0), 0};
        particle_buffer = new ComputeBuffer(n_particle, sizeof(particle));
        arg_buffer = new ComputeBuffer(1, arg.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        arg_buffer.SetData(arg);
        bound_buffer = new ComputeBuffer(n_bound, sizeof(float));
        bound_buffer.SetData(bound);
        hash_grid_buffer = new ComputeBuffer((n_cell + 1) * max_particles_per_grid, sizeof(int));
        hash_grid_tracker_buffer = new ComputeBuffer(n_cell + 1, sizeof(int));

        compute_shader.SetBuffer(malloc_particle_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(malloc_particle_kernel, "bound", bound_buffer);

        compute_shader.SetBuffer(clear_hash_grid_kernel, "hash_grid_tracker", hash_grid_tracker_buffer);

        compute_shader.SetBuffer(compute_hash_grid_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_hash_grid_kernel, "hash_grid", hash_grid_buffer);
        compute_shader.SetBuffer(compute_hash_grid_kernel, "hash_grid_tracker", hash_grid_tracker_buffer);
        compute_shader.SetBuffer(compute_hash_grid_kernel, "bound", bound_buffer);

        compute_shader.SetBuffer(compute_density_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_density_kernel, "hash_grid", hash_grid_buffer);
        compute_shader.SetBuffer(compute_density_kernel, "hash_grid_tracker", hash_grid_tracker_buffer);
        compute_shader.SetBuffer(compute_density_kernel, "bound", bound_buffer);

        compute_shader.SetBuffer(compute_acceleration_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_acceleration_kernel, "hash_grid", hash_grid_buffer);
        compute_shader.SetBuffer(compute_acceleration_kernel, "hash_grid_tracker", hash_grid_tracker_buffer);
        compute_shader.SetBuffer(compute_acceleration_kernel, "bound", bound_buffer);

        compute_shader.SetBuffer(integrate_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(integrate_kernel, "bound", bound_buffer);
    }

    void malloc_particle()
    {
        particles = new particle[n_particle];
        int i = 0;
        while(i < n_particle)
        {
            for(int x = 0; x < dimension.x; ++x)
                for(int y = 0; y < dimension.y; ++y)
                    for(int z = 0; z < dimension.z; ++z)
                    {
                        Vector3 grid_pos = new Vector3(bound[0] + x * grid_size, bound[2] + y * grid_size, bound[4] + z * grid_size);
                        for(int a = 0; a < intial_particles_per_grid; ++a)
                        {
                            Vector3 pos = grid_pos - new Vector3(grid_size_over_2, grid_size_over_2, grid_size_over_2) - new Vector3(UnityEngine.Random.Range(0f, grid_size_over_2), UnityEngine.Random.Range(0f, grid_size_over_2), UnityEngine.Random.Range(0f, grid_size_over_2));
                            particles[i] = new particle{position = pos, velocity = velocity_initial, acceleration = Vector3.zero, density = 0, pressure = 0, neighbor = 0};
                            if(++i == n_particle) return;
                        }
                    }
        }
    }


    void OnDrawGizmos()
    {
        Vector3 bound_center = new Vector3(0.5f * (bound[0] + bound[1]), 0.5f * (bound[2] + bound[3]), 0.5f * (bound[4] + bound[5])),
        bound_size = new Vector3(bound[1] - bound[0], bound[3] - bound[2], bound[5] - bound[4]);
        Gizmos.color = new Color(0, 0, 0, 1);
        Gizmos.DrawWireCube(bound_center, bound_size);
    }

    public void OnDestroy()
    {
        free();
    }

    public void free()
    {
        particle_buffer.Dispose();
        arg_buffer.Dispose();
    }
}