using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;

public class fluid : MonoBehaviour
{
    [Header("particle")]
    public Mesh particle_mesh;
    public Material material;
    public float mass = 4f,
    viscosity_coefficient = 2.5f,
    particle_size = 2f,
    radius = 1f,
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
    int dimension3;

    [Header("fluid")]
    public int n_particle = 50000,
    dimension = 100,
    max_particle_per_grid = 500;

    [StructLayout(LayoutKind.Sequential, Size = 28)]
    public struct particle
    {
        public Vector3 position;
        public Vector4 color;
    }

    int size_property = Shader.PropertyToID("size"),
    particle_buffer_property = Shader.PropertyToID("particle_buffer");

    int n = 8, /* 8 grids */
    n_bound = 6;
    particle[] particles;
    int[] neighbor_list,
    neighbor_per_particle,
    hash_grid,
    particle_per_grid;
    float[] density,
    pressure;
    Vector3[] velocity,
    force;

    ComputeBuffer particle_buffer,
    arg_buffer,
    neighbor_list_buffer,
    neighbor_per_particle_buffer,
    hash_grid_buffer,
    particle_per_grid_buffer,
    density_buffer,
    pressure_buffer,
    velocity_buffer,
    force_buffer,
    bound_buffer;

    public ComputeShader compute_shader;
    int clear_hash_grid_kernel,
    compute_hash_grid_kernel,
    compute_neighbor_list_kernel,
    compute_density_pressure_kernel,
    compute_force_kernel,
    integrate_kernel;

    void Awake()
    {
        radius2 = radius * radius;
        radius3 = radius2 * radius;
        radius4 = radius3 * radius;
        radius5 = radius4 * radius;
        mass2 = mass * mass;
        dimension3 = dimension * dimension * dimension;

        malloc_particle();
        find_kernel();
        compute_shader_init();
        compute_buffer_init();
    }

    void Update()
    {
        compute_shader.Dispatch(clear_hash_grid_kernel, dimension3 / dimension, 1, 1);
        compute_shader.Dispatch(compute_hash_grid_kernel, n_particle / dimension, 1, 1);
        compute_shader.Dispatch(compute_neighbor_list_kernel, n_particle / dimension, 1, 1);
        compute_shader.Dispatch(compute_density_pressure_kernel, n_particle / dimension, 1, 1);
        compute_shader.Dispatch(compute_force_kernel, n_particle / dimension, 1, 1);
        compute_shader.Dispatch(integrate_kernel, n_particle / dimension, 1, 1);

        material.SetFloat(size_property, particle_size);
        material.SetBuffer(particle_buffer_property, particle_buffer);
        Graphics.DrawMeshInstancedIndirect(particle_mesh, 0, material, new Bounds(Vector3.zero, new Vector3(1000f, 1000f, 1000f)), arg_buffer, castShadows: UnityEngine.Rendering.ShadowCastingMode.Off);
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
                        float r = Random.Range(0f, 5f);
                        Vector3 pos = new Vector3(1.5f * r * x, r * y, r * z) + position_offset;
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

    void find_kernel()
    {
        clear_hash_grid_kernel = compute_shader.FindKernel("clear_hash_grid");
        compute_hash_grid_kernel = compute_shader.FindKernel("compute_hash_grid");
        compute_neighbor_list_kernel = compute_shader.FindKernel("compute_neighbor_list");
        compute_density_pressure_kernel = compute_shader.FindKernel("compute_density_pressure");
        compute_force_kernel = compute_shader.FindKernel("compute_force");
        integrate_kernel = compute_shader.FindKernel("integrate");
    }

    void compute_shader_init()
    {
        compute_shader.SetFloat("grid_size", 4 * radius);
        compute_shader.SetInt("dimension", dimension);
        compute_shader.SetInt("max_particle_per_grid", max_particle_per_grid);
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
    }

    public unsafe void compute_buffer_init()
    {
        uint[] arg = {particle_mesh.GetIndexCount(0), (uint)n_particle, particle_mesh.GetIndexStart(0), particle_mesh.GetBaseVertex(0), 0};
        particle_buffer = new ComputeBuffer(n_particle, sizeof(particle));
        particle_buffer.SetData(particles);
        arg_buffer = new ComputeBuffer(1, arg.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        arg_buffer.SetData(arg);
        neighbor_list = new int[n_particle * max_particle_per_grid * n];
        neighbor_per_particle = new int[n_particle];
        hash_grid = new int[dimension3 * max_particle_per_grid];
        particle_per_grid = new int[dimension3];
        
        neighbor_list_buffer = new ComputeBuffer(n_particle * max_particle_per_grid * n, sizeof(int));
        neighbor_list_buffer.SetData(neighbor_list);
        neighbor_per_particle_buffer = new ComputeBuffer(n_particle, sizeof(int));
        neighbor_per_particle_buffer.SetData(neighbor_per_particle);
        hash_grid_buffer = new ComputeBuffer(dimension3 * max_particle_per_grid, sizeof(int));
        hash_grid_buffer.SetData(hash_grid);
        particle_per_grid_buffer = new ComputeBuffer(dimension3, sizeof(int));
        particle_per_grid_buffer.SetData(particle_per_grid);
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

        compute_shader.SetBuffer(clear_hash_grid_kernel, "particle_per_grid", particle_per_grid_buffer);

        compute_shader.SetBuffer(compute_hash_grid_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_hash_grid_kernel, "hash_grid", hash_grid_buffer);
        compute_shader.SetBuffer(compute_hash_grid_kernel, "particle_per_grid", particle_per_grid_buffer);

        compute_shader.SetBuffer(compute_neighbor_list_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "hash_grid", hash_grid_buffer);
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "particle_per_grid", particle_per_grid_buffer);
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "neighbor_list", neighbor_list_buffer);
        compute_shader.SetBuffer(compute_neighbor_list_kernel, "neighbor_per_particle", neighbor_per_particle_buffer);
        
        compute_shader.SetBuffer(compute_density_pressure_kernel, "neighbor_list", neighbor_list_buffer);
        compute_shader.SetBuffer(compute_density_pressure_kernel, "neighbor_per_particle", neighbor_per_particle_buffer);
        compute_shader.SetBuffer(compute_density_pressure_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_density_pressure_kernel, "density", density_buffer);
        compute_shader.SetBuffer(compute_density_pressure_kernel, "pressure", pressure_buffer);

        compute_shader.SetBuffer(compute_force_kernel, "neighbor_list", neighbor_list_buffer);
        compute_shader.SetBuffer(compute_force_kernel, "neighbor_per_particle", neighbor_per_particle_buffer);
        compute_shader.SetBuffer(compute_force_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_force_kernel, "density", density_buffer);
        compute_shader.SetBuffer(compute_force_kernel, "pressure", pressure_buffer);
        compute_shader.SetBuffer(compute_force_kernel, "velocity", velocity_buffer);
        compute_shader.SetBuffer(compute_force_kernel, "force", force_buffer);

        compute_shader.SetBuffer(integrate_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(integrate_kernel, "velocity", velocity_buffer);
        compute_shader.SetBuffer(integrate_kernel, "force", force_buffer);
        compute_shader.SetBuffer(integrate_kernel, "bound", bound_buffer);
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
        neighbor_per_particle_buffer.Dispose();
        hash_grid_buffer.Dispose();
        particle_per_grid_buffer.Dispose();
        density_buffer.Dispose();
        pressure_buffer.Dispose();
        velocity_buffer.Dispose();
        force_buffer.Dispose();
    }
}