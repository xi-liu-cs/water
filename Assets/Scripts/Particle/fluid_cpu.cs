using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;

public class fluid_cpu : MonoBehaviour
{
    [Header("particle")]
    public Mesh particle_mesh,
    fluid_mesh;
    public Material material;
    public float mass = 4f,
    viscosity_coefficient = 0.25f,
    particle_size = 8f,
    radius = 4f,
    grid_size = 16f, /* 4 * radius */
    gas_constant = 2000f,
    dt = 0.0008f,
    rest_density = 9f,
    damping = -1f;
    public Vector3 position_offset = new Vector3(-140, -230, -60),
    velocity_initial = new Vector3(0, 500, 0),
    g = new Vector3(0f, -9.81f * 2000f, 0f);
    public float[] bound = {-100, 200, -300, -100, -100, 100}; /* -200, 200, -300, -12, -150, 150 */
    float radius2,
    radius3,
    radius4,
    radius5,
    mass2;

    [Header("fluid")]
    public int n_particle = 50000, /* max = 130000 */
    max_particles_per_grid = 500,
    intial_particles_per_grid = 5;
    public Vector3Int dimension;

    public struct particle
    {
        public Vector3 position;
    }
    
    int[] neighbor_list,
    neighbor_tracker; /* neighbor_tracker[i] = number of neighors that particles[i] currently have */
    particle[] particles;
    Dictionary<int, List<int>> hash_grid = new Dictionary<int, List<int>>();

    ComputeBuffer particle_buffer,
    arg_buffer;
    int size_property = Shader.PropertyToID("size"),
    particle_buffer_property = Shader.PropertyToID("particle_buffer");

    int n = 8,
    grid_size_over_2,
    ceil_grid_size;
    float[] density,
    pressure;
    Vector3[] force,
    velocity;

    void Awake()
    {
        dimension = new Vector3Int((int)((bound[1] - bound[0]) / grid_size), (int)((bound[3] - bound[2]) / grid_size), (int)((bound[5] - bound[4]) / grid_size));
        hash_grid.Add(0, new List<int>());
        grid_size_over_2 = (int)grid_size / 2;
        ceil_grid_size = (int)Math.Ceiling(grid_size);
        radius2 = radius * radius;
        radius3 = radius2 * radius;
        radius4 = radius3 * radius;
        radius5 = radius4 * radius;
        malloc_particle();
        neighbor_hash_init();
        compute_buffer_init();
    }

    void malloc_particle()
    {
        particles = new particle[n_particle];
        density = new float[n_particle];
        pressure = new float[n_particle];
        velocity = new Vector3[n_particle];
        force = new Vector3[n_particle]; 
        /* int particle_per_dimension = Mathf.CeilToInt(Mathf.Pow(n_particle, 1f / 3f));
        while(i < n_particle)
        {
            for(int x = 0; x < particle_per_dimension; ++x)
                for(int y = 0; y < particle_per_dimension; ++y)
                    for(int z = 0; z < particle_per_dimension; ++z)
                    {
                        Vector3 pos = new Vector3(dimension - 1, dimension - 1, dimension - 1) - new Vector3(x / 2f, y / 2f, z / 2f)  - new Vector3(UnityEngine.Random.Range(0f, 0.01f), UnityEngine.Random.Range(0f, 0.01f), UnityEngine.Random.Range(0f, 0.01f));
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
        } */
        
        /* int i = 0;
        while(i < n_particle)
        {
            for(int x = (int)((int)(bound[0] / grid_size) * grid_size); x < bound[1]; x += ceil_grid_size)
                for(int y = (int)((int)(bound[2] / grid_size) * grid_size); y < bound[3]; y += ceil_grid_size)
                    for(int z = (int)((int)(bound[4] / grid_size) * grid_size); z < bound[5]; z += ceil_grid_size)
                        for(int a = 0; a < intial_particles_per_grid; ++a)
                        {
                            Vector3 pos = new Vector3(x, y, z) - new Vector3(grid_size_over_2, grid_size_over_2, grid_size_over_2) + new Vector3(UnityEngine.Random.Range(0f, grid_size_over_2), UnityEngine.Random.Range(0f, grid_size_over_2), UnityEngine.Random.Range(0f, grid_size_over_2));
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
        } */

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
                            particles[i] = new particle{position = pos};
                            density[i] = -1;
                            pressure[i] = 0;
                            force[i] = Vector3.zero;
                            velocity[i] = velocity_initial;
                            if(++i == n_particle) return;
                        }
                    }
        }
    }

    void neighbor_hash_init()
    {
        hash_grid.Clear();
        neighbor_list = new int[n_particle * max_particles_per_grid * n];
        neighbor_tracker = new int[n_particle];
        spatial_hash.grid_size = radius4;
        spatial_hash.dimension = dimension;
        spatial_hash.bound = bound;
        /* for(int x = (int)((int)(bound[0] / grid_size) * grid_size); x < bound[1]; x += ceil_grid_size)
            for(int y = (int)((int)(bound[2] / grid_size) * grid_size); y < bound[3]; y += ceil_grid_size)
                for(int z = (int)((int)(bound[4] / grid_size) * grid_size); z < bound[5]; z += ceil_grid_size)
                {
                    int hash = spatial_hash.hash(spatial_hash.get_cell(new Vector3Int(x, y, z)));
                    if(hash_grid.ContainsKey(hash)) continue;
                    hash_grid.Add(hash, new List<int>());
                } */
        for(int x = 0; x < dimension.x; ++x)
            for(int y = 0; y < dimension.y; ++y)
                for(int z = 0; z < dimension.z; ++z)
                {
                    int hash = spatial_hash.hash(new Vector3Int(x, y, z));
                    if(hash_grid.ContainsKey(hash)) continue;
                    hash_grid.Add(hash, new List<int>());
                }
    }

    unsafe void compute_buffer_init()
    {
        uint[] arg = {particle_mesh.GetIndexCount(0), (uint)n_particle, particle_mesh.GetIndexStart(0), particle_mesh.GetBaseVertex(0), 0};
        arg_buffer = new ComputeBuffer(1, arg.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        arg_buffer.SetData(arg);
        particle_buffer = new ComputeBuffer(n_particle, sizeof(particle));
        if(particles != null) particle_buffer.SetData(particles);
    }

    void Update()
    {
        compute_neighbor_list();
        compute_density_pressure();
        compute_force();
        integrate();
        particle_buffer.SetData(particles);
        material.SetFloat(size_property, particle_size);
        material.SetBuffer(particle_buffer_property, particle_buffer);
        Graphics.DrawMeshInstancedIndirect(particle_mesh, 0, material, new Bounds(Vector3.zero, new Vector3(100f, 100f, 100f)), arg_buffer, castShadows: UnityEngine.Rendering.ShadowCastingMode.On);
        
        /* print_hash_grid();
        int debug = 100;
        Debug.Log("neighbor list");
        for(int i = 0; i < debug; ++i)
            Debug.Log(neighbor_list[i]);
        Debug.Log("density");
        for(int i = 0; i < debug; ++i)
            Debug.Log(density[i]); */
    }

    void print_hash_grid()
    {
        foreach(var cell in hash_grid)
        {
            /* Debug.Log("cell = " + cell); */
            Debug.Log("count = " + cell.Value.Count);
            Debug.Log("key = " + cell.Key);
            for(int i = 0; i < cell.Value.Count; ++i)
                Debug.Log(String.Format("cell[{0}] = {1}", i, cell.Value[i]));
        }
    }

    /* origin_index = grid index of particles[particle_i].position 
    position = particles[particle_i].position 
    (origin_index.x + 0.5f) * grid_size = mid point of x position in origin_index-th grid 
    (origin_index.y + 0.5f) * grid_size = mid point of y position in origin_index-th grid
    (origin_index.z + 0.5f) * grid_size = mid point of z position in origin_index-th grid */
    int[] get_neighbor_key(Vector3Int origin_index, Vector3 position)
    { 
        Vector3Int[] neighbor_index = new Vector3Int[n];
        for(int i = 0; i < n; ++i)
            neighbor_index[i] = origin_index;
        if((origin_index.x + 0.5f) * grid_size <= position.x)
        {
            ++neighbor_index[4].x;
            ++neighbor_index[5].x;
            ++neighbor_index[6].x;
            ++neighbor_index[7].x;
        }
        else
        {
            --neighbor_index[4].x;
            --neighbor_index[5].x;
            --neighbor_index[6].x;
            --neighbor_index[7].x;
        }
        if((origin_index.y + 0.5f) * grid_size <= position.y)
        {
            ++neighbor_index[2].y;
            ++neighbor_index[3].y;
            ++neighbor_index[6].y;
            ++neighbor_index[7].y;
        }
        else
        {
            --neighbor_index[2].y;
            --neighbor_index[3].y;
            --neighbor_index[6].y;
            --neighbor_index[7].y;
        }
        if((origin_index.z + 0.5f) * grid_size <= position.z)
        {
            ++neighbor_index[1].z;
            ++neighbor_index[3].z;
            ++neighbor_index[5].z;
            ++neighbor_index[7].z;
        }
        else
        {
            --neighbor_index[1].z;
            --neighbor_index[3].z;
            --neighbor_index[5].z;
            --neighbor_index[7].z;
        }
        int[] neighbor_key = new int[n];
        for(int i = 0; i < n; ++i)
            neighbor_key[i] = spatial_hash.hash(neighbor_index[i]);
        return neighbor_key;
    }

    void compute_neighbor_list()
    {
        foreach(var cell in hash_grid) cell.Value.Clear();
        for(int i = 0; i < particles.Length; ++i)
        {
            /* Debug.Log("pos = " + particles[i].position); */
            int hash = spatial_hash.hash(spatial_hash.get_cell(particles[i].position));
            if(!hash_grid.ContainsKey(hash)) hash_grid.Add(hash, new List<int>());
            if(hash_grid[hash].Count == max_particles_per_grid) continue;
            hash_grid[hash].Add(i);
        }
        for(int particle_i = 0; particle_i < particles.Length; ++particle_i)
        {
            neighbor_tracker[particle_i] = 0;
            Vector3Int cell = spatial_hash.get_cell(particles[particle_i].position);
            int[] cells = get_neighbor_key(cell, particles[particle_i].position);
            for(int cell_i = 0; cell_i < cells.Length; ++cell_i)
            {
                if(!hash_grid.ContainsKey(cells[cell_i])) continue;
                List<int> neighbor_cell = hash_grid[cells[cell_i]];
                Debug.Log("neighbor cell count " + neighbor_cell.Count);
                for(int neighbor_i = 0; neighbor_i < neighbor_cell.Count; ++neighbor_i)
                {
                    int potential_neighbor = neighbor_cell[neighbor_i];
                    if(potential_neighbor == particle_i) continue;
                    Debug.Log(String.Format("potential_neighbor: particles[{0}].position = {1}, particle_i: particles[{2}].position = {3}, distance = {4}", potential_neighbor, particles[potential_neighbor].position, particle_i, particles[particle_i].position, (particles[potential_neighbor].position - particles[particle_i].position).sqrMagnitude));
                    /* Debug.Log("neighbor_i = " + neighbor_i);
                    Debug.Log("potential_neighbor = " + potential_neighbor);
                    Debug.Log("potential = " + particles[potential_neighbor].position);
                    Debug.Log("particles[particle_i].position = " + particles[particle_i].position);
                    Debug.Log(String.Format("potential_neighbor: particles[{0}].position = {1}, particle_i: particles[{2}].position = {3}, distance = {4}", potential_neighbor, particles[potential_neighbor].position, particle_i, particles[particle_i].position, (particles[potential_neighbor].position - particles[particle_i].position).sqrMagnitude)); */
                    if((particles[potential_neighbor].position - particles[particle_i].position).sqrMagnitude < radius2)
                        neighbor_list[particle_i * max_particles_per_grid * n + neighbor_tracker[particle_i]++] = potential_neighbor;
                }
            }
        }
    }

    float std_kernel(float distance_square)
    {
        float x = 1.0f - distance_square / radius2;
        return 315.0f / (64.0f * Mathf.PI * radius3) * x * x * x;
    }

    float spiky_kernel_first_derivative(float distance)
    {
        float x = 1.0f - distance / radius;
        return -45.0f / (Mathf.PI * radius4) * x * x;
    }

    float spiky_kernel_second_derivative(float distance)
    {
        float x = 1.0f - distance / radius;
        return 90.0f / (Mathf.PI * radius5) * x;
    }

    Vector3 spiky_kernel_gradient(float distance, Vector3 direction_from_center)
    {
        return spiky_kernel_first_derivative(distance) * direction_from_center;
    }

    void compute_density_pressure()
    {
        for(int i = 0; i < particles.Length; ++i)
        {
            Vector3 origin = particles[i].position;
            float sum = 0;
            for(int j = 0; j < neighbor_tracker[i]; ++j)
            {
                int neighbor_index = neighbor_list[i * max_particles_per_grid * n + j];
                Vector3 diff = origin - particles[neighbor_index].position;
                float distance_square = diff.sqrMagnitude;
                sum += std_kernel(distance_square);
            }
            density[i] = sum * mass + 0.000001f;
            pressure[i] = gas_constant * (density[i] - rest_density);
        }
    }

    void compute_force()
    {
        for(int i = 0; i < particles.Length; ++i)
        {
            force[i] = Vector3.zero;
            float particle_density2 = density[i] * density[i];
            for(int j = 0; j < neighbor_tracker[i]; ++j)
            {
                int neighbor_index = neighbor_list[i * max_particles_per_grid * n + j];
                float distance = (particles[i].position - particles[neighbor_index].position).magnitude;
                if(distance > 0)
                {
                    Vector3 direction = (particles[i].position - particles[neighbor_index].position) / distance;
                    force[i] -= mass2 * (pressure[i] / particle_density2 + pressure[neighbor_index] / (density[neighbor_index] * density[neighbor_index])) * spiky_kernel_gradient(distance, direction); /* compute pressure gradient force */
                    force[i] += viscosity_coefficient * mass2 * (velocity[neighbor_index] - velocity[i]) / density[neighbor_index] * spiky_kernel_second_derivative(distance);
                }
            }
            force[i] += g;
        }
    }

    void integrate()
    {
        for(int i = 0; i < particles.Length; ++i)
        {
            particle particle = particles[i];
            velocity[i] += dt * force[i] / mass;
            particle.position += dt * velocity[i];
            particles[i] = particle;
            particle = particles[i];
            Vector3 v = velocity[i];
            if(particles[i].position.x < bound[0] + float.Epsilon)
            {
                v.x *= damping;
                particle.position.x = bound[0] + float.Epsilon;
            }
            else if(particles[i].position.x > bound[1] - float.Epsilon) 
            {
                v.x *= damping;
                particle.position.x = bound[1] - float.Epsilon;
            }
            if(particles[i].position.y < bound[2] + float.Epsilon)
            {
                v.y *= damping;
                particle.position.y = bound[2] + float.Epsilon;
            }
            else if(particles[i].position.y > bound[3] - float.Epsilon) 
            {
                v.y *= damping;
                particle.position.y = bound[3] - float.Epsilon;
            }
            if(particles[i].position.z < bound[4] + float.Epsilon)
            {
                v.z *= damping;
                particle.position.z = bound[4] + float.Epsilon;
            }
            else if(particles[i].position.z > bound[5] - float.Epsilon) 
            {
                v.z *= damping;
                particle.position.z = bound[5] - float.Epsilon;
            }
            velocity[i] = v;
            particles[i] = particle;
        }
    }

    void OnDestroy()
    {
        particle_buffer.Dispose();
        arg_buffer.Dispose();
    }
}