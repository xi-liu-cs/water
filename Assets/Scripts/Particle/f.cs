using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;

public class f : MonoBehaviour
{
    [Header("particle")]
    public Mesh particle_mesh;
    public Material material;
    public float mass = 1f,
    viscosity_coefficient = 0.1f,
    particle_size = 2f,
    radius = 2f, /* h, smoothing length */
    grid_size = 1f,
    gas_constant = 2000f,
    dt = 0.0001f,
    rest_density = 1000f,
    damping = -0.5f,
    g = -9.81f;
    public Vector3 position_offset = new Vector3(0, 0, 0),
    velocity_initial = new Vector3(0, 0, 0);
    public float[] bound = {-50, 50, -50, 50, -50, 50}; /* {-160, 140, -250, -150, -100, 100} */
    [HideInInspector]
    public float radius2,
    radius3,
    radius4,
    radius5,
    radius8,
    mass2;

    [Header("fluid")]
    public int n_particle = 100000; /* max = 1000000 */
    public Vector3Int dimension;
    int dimension2,
    dimension3;
    public int[] dimension_array;

    /* need to match compute shader file */
    int thread_per_group = 512,
    dispatch_size,
    cell_dispatch_size,
    neighbor_dispatch_size,
    n_bound = 6;
    bool exclusive = true;

    public struct particle /* struct alignment */
    {
        public Vector3 position,
        velocity,
        acceleration;
        public float mass,
        density,
        pressure;
        public int neighbor,
        cell_index,
        index_in_cell;
    };

    int size_property = Shader.PropertyToID("size"),
    particle_buffer_property = Shader.PropertyToID("particle_buffer");

    public ComputeBuffer particle_buffer,
    arg_buffer,
    bound_buffer,
    aux_buffer,
    cell_particle_count_buffer,
    cell_offset_buffer,
    sort_particle_index_buffer,
    neighbor_count_buffer,
    neighbor_offset_buffer,
    neighbor_buffer,
    total_neighbor_count_buffer;

    public ComputeShader compute_shader;
    int malloc_particle_kernel,
    clear_count_kernel,
    compute_count_kernel,
    scan_kernel,
    scan_inclusive_kernel,
    scan_exclusive_kernel,
    scan_result_kernel,
    scan_add_result_kernel,
    sequential_scan_kernel,
    count_sort_particle_index_kernel,
    compute_neighbor_count_kernel,
    compute_total_neighbor_count_kernel,
    compute_neighbor_kernel,
    debug_test_kernel,
    compute_density_kernel,
    reset_acceleration_kernel,
    compute_acceleration_kernel,
    integrate_kernel;

    public float max_density = 0,
    bulk_modulus = 1000f;

    particle[] particles;
    int[] sort_particle,
    cell_particle_count,
    cell_offset,
    neighbor_count,
    neighbor_offset,
    neighbors;
    int total_neighbor;

    public void Awake()
    {
        dimension = new Vector3Int((int)((bound[1] - bound[0]) / grid_size), (int)((bound[3] - bound[2]) / grid_size), (int)((bound[5] - bound[4]) / grid_size));
        dimension_array = new int[]{dimension.x, dimension.y, dimension.z};
        dimension3 = dimension.x * dimension.y * dimension.z;
        radius2 = radius * radius;
        radius3 = radius2 * radius;
        radius4 = radius3 * radius;
        radius5 = radius4 * radius;
        radius8 = radius2 * radius4;
        mass2 = mass * mass;
        dispatch_size = Mathf.CeilToInt((float)n_particle / (float)thread_per_group);
        cell_dispatch_size = Mathf.CeilToInt((float)dimension3 / (float)thread_per_group);

        find_kernel();
        scan_kernel = exclusive ? scan_exclusive_kernel : scan_inclusive_kernel;
        compute_shader_init();
        compute_buffer_init();
        compute_shader.Dispatch(malloc_particle_kernel, thread_per_group, 1, 1);
    }

    unsafe public void Update() /* debug_particle(); */
    {
        compute_shader.SetBuffer(clear_count_kernel, "cell_particle_count", cell_particle_count_buffer);
        compute_shader.Dispatch(clear_count_kernel, cell_dispatch_size, 1, 1);

        compute_shader.SetBuffer(compute_count_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_count_kernel, "bound", bound_buffer);
        compute_shader.SetBuffer(compute_count_kernel, "cell_particle_count", cell_particle_count_buffer);
        compute_shader.Dispatch(compute_count_kernel, dispatch_size, 1, 1);

        /* debug_compute_count();
        cpu_sequential_cell_offset();
        debug_test();
        sequential_cell_offset(); debug_cell_offset(); */

        scan(cell_particle_count_buffer, cell_offset_buffer, dimension3); /* debug_cell_offset(); */

        /* compare_gpu_cpu_cell_offset_exclusive(); */

        compute_shader.SetBuffer(count_sort_particle_index_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(count_sort_particle_index_kernel, "bound", bound_buffer);
        compute_shader.SetBuffer(count_sort_particle_index_kernel, "cell_offset", cell_offset_buffer);
        compute_shader.SetBuffer(count_sort_particle_index_kernel, "sort_particle_index", sort_particle_index_buffer);
        compute_shader.Dispatch(count_sort_particle_index_kernel, dispatch_size, 1, 1);

        /* debug_sort_particle_index(); */
        
        compute_shader.SetBuffer(compute_neighbor_count_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(compute_neighbor_count_kernel, "bound", bound_buffer);
        compute_shader.SetBuffer(compute_neighbor_count_kernel, "cell_particle_count", cell_particle_count_buffer);
        compute_shader.SetBuffer(compute_neighbor_count_kernel, "cell_offset", cell_offset_buffer);
        compute_shader.SetBuffer(compute_neighbor_count_kernel, "sort_particle_index", sort_particle_index_buffer);
        compute_shader.SetBuffer(compute_neighbor_count_kernel, "neighbor_count", neighbor_count_buffer);
        compute_shader.Dispatch(compute_neighbor_count_kernel, dispatch_size, 1, 1);

        /* debug_neighbor_count();
        sequential_neighbor_offset(); */

        scan(neighbor_count_buffer, neighbor_offset_buffer, n_particle);

        /* debug_neighbor_offset();
        compare_gpu_cpu_neighbor_offset_exclusive(); */
        
        int[] total_neighbor_count = new int[1];
        neighbor_offset_buffer.GetData(total_neighbor_count, 0, n_particle - 1, 1);
        total_neighbor = total_neighbor_count[0];
        /* Debug.Log("total_neighbor_count = " + total_neighbor_count[0]); */

        if(total_neighbor > 0)
        {
            neighbor_buffer = new ComputeBuffer(total_neighbor_count[0], sizeof(int));
            compute_shader.SetBuffer(compute_neighbor_kernel, "particles", particle_buffer);
            compute_shader.SetBuffer(compute_neighbor_kernel, "bound", bound_buffer);
            compute_shader.SetBuffer(compute_neighbor_kernel, "cell_particle_count", cell_particle_count_buffer);
            compute_shader.SetBuffer(compute_neighbor_kernel, "cell_offset", cell_offset_buffer);
            compute_shader.SetBuffer(compute_neighbor_kernel, "sort_particle_index", sort_particle_index_buffer);
            compute_shader.SetBuffer(compute_neighbor_kernel, "neighbor_count", neighbor_count_buffer);
            compute_shader.SetBuffer(compute_neighbor_kernel, "neighbor_offset", neighbor_offset_buffer);
            compute_shader.SetBuffer(compute_neighbor_kernel, "neighbors", neighbor_buffer);
            compute_shader.Dispatch(compute_neighbor_kernel, dispatch_size, 1, 1);
            /* debug_neighbor(total_neighbor_count[0]); */

            compute_shader.SetBuffer(compute_density_kernel, "particles", particle_buffer);
            compute_shader.SetBuffer(compute_density_kernel, "neighbor_offset", neighbor_offset_buffer);
            compute_shader.SetBuffer(compute_density_kernel, "neighbors", neighbor_buffer);
            compute_shader.Dispatch(compute_density_kernel, dispatch_size, 1, 1);

            compute_shader.SetBuffer(reset_acceleration_kernel, "particles", particle_buffer);
            compute_shader.Dispatch(reset_acceleration_kernel, dispatch_size, 1, 1);

            compute_shader.SetBuffer(compute_acceleration_kernel, "particles", particle_buffer);
            compute_shader.SetBuffer(compute_acceleration_kernel, "bound", bound_buffer);
            compute_shader.SetBuffer(compute_acceleration_kernel, "neighbor_offset", neighbor_offset_buffer);
            compute_shader.SetBuffer(compute_acceleration_kernel, "neighbors", neighbor_buffer);
            compute_shader.Dispatch(compute_acceleration_kernel, dispatch_size, 1, 1);
        }
        else
            Debug.Log("total neighbor is 0");

        compute_shader.SetBuffer(integrate_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(integrate_kernel, "bound", bound_buffer);
        compute_shader.Dispatch(integrate_kernel, dispatch_size, 1, 1);

        material.SetFloat(size_property, particle_size);
        material.SetBuffer(particle_buffer_property, particle_buffer); 
        Graphics.DrawMeshInstancedIndirect(particle_mesh, 0, material, new Bounds(Vector3.zero, new Vector3(1000f, 1000f, 1000f)), arg_buffer, castShadows: UnityEngine.Rendering.ShadowCastingMode.Off);
        debug_particle_neighbor();
    }

    void debug_particle()
    {
        particles = new particle[n_particle];
        particle_buffer.GetData(particles);
        for(int i = 0; i < 100; ++i)
            Debug.Log("position = " + particles[i].position
            + ", velocity = " + particles[i].velocity
            + ", acceleration = " + particles[i].acceleration
            + ", density = " + particles[i].density
            + ", pressure = " + particles[i].pressure
            + ", neighbor = " + particles[i].neighbor
            + ", cell_index = " + particles[i].cell_index
            + ", index_in_cell = " + particles[i].index_in_cell);
    }

    void debug_particle_neighbor()
    {
        particles = new particle[n_particle];
        particle_buffer.GetData(particles);
        neighbor_offset = new int[n_particle];
        neighbor_offset_buffer.GetData(neighbor_offset);
        neighbors = new int[total_neighbor];
        neighbor_buffer.GetData(neighbors);
        for(int i = 0; i < 100; ++i)
        {
            string s = "position = " + particles[i].position
            + ", velocity = " + particles[i].velocity
            + ", acceleration = " + particles[i].acceleration
            + ", density = " + particles[i].density
            + ", pressure = " + particles[i].pressure
            + ", neighbor = " + particles[i].neighbor
            + ", cell_index = " + particles[i].cell_index
            + ", index_in_cell = " + particles[i].index_in_cell;
            s += ", neighbor positions = ";
            for(int j = 0; j < particles[i].neighbor; ++j)
            {
                Vector3 neighbor_position = particles[neighbors[neighbor_offset[i] + j]].position;
                float distance = Vector3.Distance(neighbor_position, particles[i].position);
                s += particles[neighbors[neighbor_offset[i] + j]].position + ", d = " + distance + ", ";
            }
            Debug.Log(s);
        }
    }

    void debug_compute_count()
    {
        cell_particle_count = new int[dimension3];
        cell_particle_count_buffer.GetData(cell_particle_count);
        for(int i = 0; i < 1000; ++i)
            if(cell_particle_count[i] != 0)
                Debug.LogFormat("cell_particle_count[{0}] = {1}", i, cell_particle_count[i]);
    }

    void debug_test()
    {
        compute_shader.SetInt("sequential_scan_length", dimension3);
        compute_shader.SetBuffer(debug_test_kernel, "scan_in", cell_particle_count_buffer);
        compute_shader.SetBuffer(debug_test_kernel, "scan_out", cell_offset_buffer);
        compute_shader.Dispatch(debug_test_kernel, 1, 1, 1);
        cell_offset = new int[dimension3];
        cell_offset_buffer.GetData(cell_offset);
        for(int i = 0; i < 1000; ++i)
            if(cell_offset[i] != 0)
                Debug.LogFormat("cell_offset[{0}] = {1}", i, cell_offset[i]);
    }

    void cpu_sequential_cell_offset()
    {
        cell_particle_count = new int[dimension3];
        cell_particle_count_buffer.GetData(cell_particle_count);
        int[] scan_out = new int[dimension3];
        scan_out[0] = cell_particle_count[0];
        for(int i = 1; i < 1000; ++i)
        {
            scan_out[i] = scan_out[i - 1] + cell_particle_count[i];
            Debug.LogFormat("scan_out[{0}] = {1}", i - 1, scan_out[i - 1]);
        }
    }

    void compare_gpu_cpu_cell_offset()
    {
        cell_offset = new int[dimension3];
        cell_offset_buffer.GetData(cell_offset);
        cell_particle_count = new int[dimension3];
        cell_particle_count_buffer.GetData(cell_particle_count);
        int[] scan_out = new int[dimension3];
        scan_out[0] = cell_particle_count[0];
        for(int i = 1; i < 1000; ++i)
        {
            scan_out[i] = scan_out[i - 1] + cell_particle_count[i];
            Debug.LogFormat("i = {0}, cpu_scan_out = {1}, gpu_cell_offset = {2}, diff = {3}", i - 1, scan_out[i - 1], cell_offset[i - 1], scan_out[i - 1] - cell_offset[i - 1]);
        }
    }

    void compare_gpu_cpu_cell_offset_exclusive()
    {
        cell_offset = new int[dimension3];
        cell_offset_buffer.GetData(cell_offset);
        cell_particle_count = new int[dimension3];
        cell_particle_count_buffer.GetData(cell_particle_count);
        int[] scan_out = new int[dimension3];
        scan_out[0] = 0;
        for(int i = 1; i < 1000; ++i)
        {
            scan_out[i] = scan_out[i - 1] + cell_particle_count[i - 1];
            Debug.LogFormat("i = {0}, cpu_scan_out = {1}, gpu_cell_offset = {2}, diff = {3}", i - 1, scan_out[i - 1], cell_offset[i - 1], scan_out[i - 1] - cell_offset[i - 1]);
        }
    }

    void sequential_cell_offset()
    {
        compute_shader.SetInt("sequential_scan_length", dimension3);
        compute_shader.SetBuffer(sequential_scan_kernel, "scan_in", cell_particle_count_buffer);
        compute_shader.SetBuffer(sequential_scan_kernel, "scan_out", cell_offset_buffer);
        compute_shader.Dispatch(sequential_scan_kernel, 1, 1, 1);
    }

    void debug_cell_offset()
    {
        cell_offset = new int[dimension3];
        cell_offset_buffer.GetData(cell_offset);
        for(int i = 0; i < 1000; ++i)
            if(cell_offset[i] != 0)
                Debug.LogFormat("cell_offset[{0}] = {1}", i, cell_offset[i]);
    }

    void debug_sort_particle_index()
    {
        particles = new particle[n_particle];
        particle_buffer.GetData(particles);
        
        sort_particle = new int[n_particle];
        sort_particle_index_buffer.GetData(sort_particle);

        for(int i = 0; i < 100; ++i)
        {
            int id = sort_particle[i];
            Debug.Log("id = " + id
            + ", position = " + particles[id].position
            + ", velocity = " + particles[id].velocity
            + ", acceleration = " + particles[id].acceleration
            + ", density = " + particles[id].density
            + ", pressure = " + particles[id].pressure
            + ", cell_index = " + particles[id].cell_index
            + ", index_in_cell = " + particles[id].index_in_cell);
        }
    }

    void debug_neighbor_count()
    {
        neighbor_count = new int[n_particle];
        neighbor_count_buffer.GetData(neighbor_count);
        for(int i = 0; i < 1000; ++i)
            Debug.LogFormat("neighbor_count[{0}] = {1}", i, neighbor_count[i]);
    }

    void debug_neighbor_offset()
    {
        neighbor_offset = new int[n_particle];
        neighbor_offset_buffer.GetData(neighbor_offset);
        for(int i = 0; i < 1000; ++i)
            Debug.LogFormat("neighbor_offset[{0}] = {1}", i, neighbor_offset[i]);
    }

    void compare_gpu_cpu_neighbor_offset()
    {
        neighbor_offset = new int[n_particle];
        neighbor_offset_buffer.GetData(neighbor_offset);
        neighbor_count = new int[n_particle];
        neighbor_count_buffer.GetData(neighbor_count);
        int[] scan_out = new int[n_particle];
        scan_out[0] = neighbor_count[0];
        for(int i = 1; i < 1000; ++i)
        {
            scan_out[i] = scan_out[i - 1] + neighbor_count[i];
            Debug.LogFormat("i = {0}, cpu_scan_out = {1}, gpu_neighbor_offset = {2}, diff = {3}", i - 1, scan_out[i - 1], neighbor_offset[i - 1], scan_out[i - 1] - neighbor_offset[i - 1]);
        }
    }

    void compare_gpu_cpu_neighbor_offset_exclusive()
    {
        neighbor_offset = new int[n_particle];
        neighbor_offset_buffer.GetData(neighbor_offset);
        neighbor_count = new int[n_particle];
        neighbor_count_buffer.GetData(neighbor_count);
        int[] scan_out = new int[n_particle];
        scan_out[0] = 0;
        for(int i = 1; i < 1000; ++i)
        {
            scan_out[i] = scan_out[i - 1] + neighbor_count[i - 1];
            Debug.LogFormat("i = {0}, cpu_scan_out = {1}, gpu_neighbor_offset = {2}, diff = {3}", i - 1, scan_out[i - 1], neighbor_offset[i - 1], scan_out[i - 1] - neighbor_offset[i - 1]);
        }
    }

    void sequential_neighbor_offset()
    {
        compute_shader.SetInt("sequential_scan_length", n_particle);
        compute_shader.SetBuffer(sequential_scan_kernel, "scan_in", neighbor_count_buffer);
        compute_shader.SetBuffer(sequential_scan_kernel, "scan_out", neighbor_offset_buffer);
        compute_shader.Dispatch(sequential_scan_kernel, 1, 1, 1);
    }

    void debug_total_neighbor_count()
    {
        compute_shader.SetBuffer(compute_total_neighbor_count_kernel, "neighbor_offset", neighbor_offset_buffer);
        compute_shader.SetBuffer(compute_total_neighbor_count_kernel, "total_neighbor_count", total_neighbor_count_buffer);
        compute_shader.Dispatch(compute_total_neighbor_count_kernel, 1, 1, 1);
        int[] total_neighbor_count = new int[1];
        total_neighbor_count_buffer.GetData(total_neighbor_count);
        Debug.Log("total_neighbor_count = " + total_neighbor_count[0]);
    }

    void debug_neighbor(int n_neighbor)
    {
        neighbors = new int[n_neighbor];
        neighbor_buffer.GetData(neighbors);
        for(int i = 0; i < 1000; ++i)
            Debug.LogFormat("neighbor_buffer[{0}] = {1}", i, neighbors[i]);
    }

    void find_kernel()
    {
        malloc_particle_kernel = compute_shader.FindKernel("malloc_particle");
        clear_count_kernel = compute_shader.FindKernel("clear_count");
        compute_count_kernel = compute_shader.FindKernel("compute_count");
        scan_inclusive_kernel = compute_shader.FindKernel("scan_inclusive");
        scan_exclusive_kernel = compute_shader.FindKernel("scan_exclusive");
        scan_result_kernel = compute_shader.FindKernel("scan_result");
        scan_add_result_kernel = compute_shader.FindKernel("scan_add_result");
        sequential_scan_kernel = compute_shader.FindKernel("sequential_scan");
        count_sort_particle_index_kernel = compute_shader.FindKernel("count_sort_particle_index");
        compute_neighbor_count_kernel = compute_shader.FindKernel("compute_neighbor_count");
        compute_neighbor_kernel = compute_shader.FindKernel("compute_neighbor");
        debug_test_kernel = compute_shader.FindKernel("debug_test");
        compute_density_kernel = compute_shader.FindKernel("compute_density");
        reset_acceleration_kernel = compute_shader.FindKernel("reset_acceleration");
        compute_acceleration_kernel = compute_shader.FindKernel("compute_acceleration");
        integrate_kernel = compute_shader.FindKernel("integrate");
    }

    void compute_shader_init()
    {
        compute_shader.SetFloat("grid_size", grid_size);
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
        compute_shader.SetFloat("g", g);
        compute_shader.SetFloat("epsilon", Mathf.Epsilon);
        compute_shader.SetFloat("e", Mathf.Exp(1));
        compute_shader.SetFloat("pi", Mathf.PI);
        compute_shader.SetFloat("bulk_modulus", bulk_modulus);
        compute_shader.SetVector("time", Shader.GetGlobalVector("_Time"));
        compute_shader.SetInts("dimension", dimension_array);
        compute_shader.SetInt("n_particle", n_particle);
        compute_shader.SetInt("seed", UnityEngine.Random.Range(0, int.MaxValue));
    }

    unsafe void compute_buffer_init()
    {
        uint[] arg = {particle_mesh.GetIndexCount(0), (uint)n_particle, particle_mesh.GetIndexStart(0), particle_mesh.GetBaseVertex(0), 0};
        particle_buffer = new ComputeBuffer(n_particle, sizeof(particle));
        arg_buffer = new ComputeBuffer(1, arg.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        arg_buffer.SetData(arg);
        bound_buffer = new ComputeBuffer(n_bound, sizeof(float));
        bound_buffer.SetData(bound);
        cell_particle_count_buffer = new ComputeBuffer(dimension3, sizeof(int));
        cell_offset_buffer = new ComputeBuffer(dimension3, sizeof(int));
        sort_particle_index_buffer = new ComputeBuffer(n_particle, sizeof(int));
        neighbor_count_buffer = new ComputeBuffer(n_particle, sizeof(int));
        neighbor_offset_buffer = new ComputeBuffer(n_particle, sizeof(int));
        total_neighbor_count_buffer = new ComputeBuffer(1, sizeof(int));

        compute_shader.SetBuffer(malloc_particle_kernel, "particles", particle_buffer);
        compute_shader.SetBuffer(malloc_particle_kernel, "bound", bound_buffer);
    }

    void scan(ComputeBuffer in_buffer, ComputeBuffer out_buffer, int scan_length)
    {
        aux_buffer = new ComputeBuffer(scan_length, sizeof(int));
        int scan_thread_group = Mathf.CeilToInt((float)scan_length / (float)thread_per_group);
        compute_shader.SetBuffer(scan_kernel, "scan_in", in_buffer);
        compute_shader.SetBuffer(scan_kernel, "scan_out", out_buffer);
        compute_shader.Dispatch(scan_kernel, scan_thread_group, 1, 1);
        
        compute_shader.SetBuffer(scan_result_kernel, "scan_in", out_buffer);
        compute_shader.SetBuffer(scan_result_kernel, "scan_out", aux_buffer);
        compute_shader.Dispatch(scan_result_kernel, 1, 1, 1);
        
        compute_shader.SetBuffer(scan_add_result_kernel, "scan_in", aux_buffer);
        compute_shader.SetBuffer(scan_add_result_kernel, "scan_out", out_buffer);
        compute_shader.Dispatch(scan_add_result_kernel, scan_thread_group, 1, 1);
        aux_buffer.Dispose();
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
        bound_buffer.Dispose();
        aux_buffer.Dispose();
        cell_particle_count_buffer.Dispose();
        cell_offset_buffer.Dispose();
        sort_particle_index_buffer.Dispose();
        neighbor_count_buffer.Dispose();
        neighbor_offset_buffer.Dispose();
        neighbor_buffer.Dispose();
        total_neighbor_count_buffer.Dispose();
    }
}