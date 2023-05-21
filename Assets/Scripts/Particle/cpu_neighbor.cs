using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;

public class cpu_neighbor : MonoBehaviour
{
    struct particle
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
    }

    public int n_particle = 1000;
    public float[] bound = new float[]{-50, 50, -50, 50, -50, 50};
    public Vector3Int dimension;
    int dimension2,
    dimension3;
    public int[] dimension_array;
    particle[] particles;
    int[] cell_particle_count,
    cell_offset,
    sort_particle_index,
    neighbor_count,
    neighbor_offset,
    neighbors;
    int n_cell,
    total_particle_in_cell, /* 1 particle can be in multiple cells */
    total_neighbor;

    public float mass = 1f,
    viscosity_coefficient = 10f,
    particle_size = 2f,
    radius = 2f, /* h, smoothing length */
    grid_size = 2f,
    gas_constant = 2000f,
    dt = 0.01f,
    rest_density = 1f,
    damping = -0.5f,
    g = -9.81f,
    epsilon = Mathf.Epsilon,
    pi = Mathf.PI,
    bulk_modulus = 1000f;
    public float radius2,
    radius3,
    radius4,
    radius5,
    radius8,
    mass2;

    void Start()
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
        n_cell = (int)((bound[1] - bound[0]) / grid_size) + dimension.x * ((int)((bound[3] - bound[2]) / grid_size) + dimension.y * (int)((bound[5] - bound[4]) / grid_size));
        cell_particle_count = new int[n_cell + 1];
        cell_offset = new int[n_cell + 1];
        neighbor_count = new int[n_particle];
        neighbor_offset = new int[n_particle];
        malloc_particle();
    }

    void Update()
    {
        /* debug_compute_count(); */
        compute_count();
        /* Debug.Log("cell particle count"); debug_int(cell_particle_count); */
        scan(cell_particle_count, cell_offset);
        /* Debug.Log("cell offset"); debug_int(cell_offset); */
        total_particle_in_cell = cell_offset[cell_offset.Length - 2] + cell_offset[cell_offset.Length - 1];
        sort_particle_index = new int[total_particle_in_cell];
        /* debug_count_sort_particle_index(); */
        count_sort_particle_index();
        /* Debug.Log("sort_particle_index"); debug_int(sort_particle_index); */
        compute_neighbor_count();
        /* Debug.Log("neighbor_count"); debug_int(neighbor_count); */
        scan(neighbor_count, neighbor_offset);
        /* Debug.Log("neighbor_offset"); debug_int(neighbor_offset); */
        total_neighbor = neighbor_offset[n_particle - 2] + neighbor_offset[n_particle - 1];
        /* Debug.LogFormat("total neighbor = {0}", total_neighbor); */
        neighbors = new int[total_neighbor];
        compute_neighbor();
        compute_density();
        compute_acceleration();
        /* debug_integrate(); */
        integrate();
        debug_particle_neighbor();
    }

    void malloc_particle()
    {
        particles = new particle[n_particle];
        for(int i = 0; i < n_particle; ++i)
        {
            particles[i].position = new Vector3(0.2f * UnityEngine.Random.Range(bound[0], bound[1]), 0.2f * UnityEngine.Random.Range(bound[2], bound[3]), 0.2f * UnityEngine.Random.Range(bound[4], bound[5]));
            /* if(float.IsNaN(particles[i].position[0]) || float.IsNaN(particles[i].position[1]) || float.IsNaN(particles[i].position[2])) Debug.LogFormat("particles[{0}].position = {1}", i, particles[i].position); */
            particles[i].velocity = new Vector3(0, 0, 0);
            particles[i].acceleration = new Vector3(0, 0, 0);
            particles[i].mass = 1;
            particles[i].density = rest_density;
            particles[i].pressure = 0;
            particles[i].neighbor = 0;
            particles[i].cell_index = hash(get_cell(particles[i].position));
            particles[i].index_in_cell = 0;
        }
    }

    void debug_particle_neighbor()
    {
        for(int i = 0; i < n_particle; ++i)
            Debug.LogFormat("position = {0}, neighbor = {1}", particles[i].position, particles[i].neighbor);
    }

    Vector3Int get_cell(Vector3 position)
    {
        return new Vector3Int((int)((position.x - bound[0]) / grid_size), (int)((position.y - bound[2]) / grid_size), (int)((position.z - bound[4]) / grid_size));
    }

    int hash(Vector3Int cell)
    {
        return cell.x + dimension.x * (cell.y + dimension.y * cell.z);
    }

    void scan(int[] scan_in, int[] scan_out)
    {
        scan_out[0] = scan_in[0];
        for(int i = 1; i < scan_in.Length; ++i)
            scan_out[i] = scan_out[i - 1] + scan_in[i - 1];
    }

    void debug_int(int[] a)
    {
        string s = "";
        int count = 0;
        for(int i = 0; i < a.Length; ++i)
            if(a[i] != 0 && count++ < 500)
                s += a[i] + " ";
        Debug.Log(s);
    }

    void compute_count()
    {
        for(int i = 0; i < n_particle; ++i)
        {
            if(particles[i].cell_index >= cell_particle_count.Length || particles[i].cell_index < 0) Debug.LogFormat("particles[{0}].cell_index = {1}, cell_particle_count.Length = {2}, position = {3}", i, particles[i].cell_index, cell_particle_count.Length, particles[i].position);
            particles[i].index_in_cell = cell_particle_count[particles[i].cell_index]++;
        }
    }

    void debug_compute_count()
    {
        Debug.LogFormat("cell_particle_count.Length = {0}", cell_particle_count.Length);
        string s = "";
        for(int i = 0; i < n_particle; ++i)
        {
            particles[i].index_in_cell = cell_particle_count[particles[i].cell_index]++;
            s += particles[i].cell_index + " ";
        }
        Debug.Log(s);
    }

    void count_sort_particle_index()
    {
        for(int i = 0; i < n_particle; ++i)
            sort_particle_index[cell_offset[particles[i].cell_index] + particles[i].index_in_cell] = i;
    }

    void debug_count_sort_particle_index()
    {
        for(int i = 0; i < n_particle; ++i)
        {
            int j = cell_offset[particles[i].cell_index] + particles[i].index_in_cell;
            if(j >= sort_particle_index.Length || j < 0)
                Debug.LogFormat("j = {0}, sort_particle_index.Length = {1}, cell_offset[particles[i].cell_index] = {2}, particles[i].index_in_cell = {3}", j, sort_particle_index.Length, cell_offset[particles[i].cell_index], particles[i].index_in_cell);
            sort_particle_index[j] = i;
        }
    }

    void compute_neighbor_count()
    {
        for(int i = 0; i < n_particle; ++i)
        {
            Vector3Int origin_index = get_cell(particles[i].position);
            int count = 0;
            for(int x = -1; x < 2; ++x)
                for(int y = -1; y < 2; ++y)
                    for(int z = -1; z < 2; ++z)
                    {
                        Vector3Int cell = origin_index + new Vector3Int(x, y, z);
                        if(cell.x < 0 || cell.x >= dimension.x || cell.y < 0 || cell.y >= dimension.y || cell.z < 0 || cell.z >= dimension.z) continue;
                        int cell_index = hash(cell),
                        cell_count = cell_particle_count[cell_index],
                        cell_start = cell_offset[cell_index];
                        for(int j = cell_start; j < cell_start + cell_count; ++j)
                        {
                            Vector3 diff = particles[sort_particle_index[j]].position - particles[i].position;
                            float square_distance = Vector3.Dot(diff, diff);
                            if(square_distance < radius2)
                                ++count;
                        }
                    }
            particles[i].neighbor = count;
            neighbor_count[i] = count;
        }
    }

    void compute_neighbor()
    {
        for(int i = 0; i < n_particle; ++i)
        {
            Vector3Int origin_index = get_cell(particles[i].position);
            int count = 0,
            write_offset = neighbor_offset[i];
            for(int x = -1; x < 2; ++x)
                for(int y = -1; y < 2; ++y)
                    for(int z = -1; z < 2; ++z)
                    {
                        Vector3Int cell = origin_index + new Vector3Int(x, y, z);
                        if(cell.x < 0 || cell.x >= dimension.x || cell.y < 0 || cell.y >= dimension.y || cell.z < 0 || cell.z >= dimension.z) continue;
                        int cell_index = hash(cell),
                        cell_count = cell_particle_count[cell_index],
                        cell_start = cell_offset[cell_index];
                        for(int j = cell_start; j < cell_start + cell_count; ++j)
                        {
                            Vector3 diff = particles[sort_particle_index[j]].position - particles[i].position;
                            float square_distance = Vector3.Dot(diff, diff);
                            if(square_distance < radius2)
                            {
                                neighbors[write_offset + count] = sort_particle_index[j];
                                ++count;
                            }
                        }
                    }
        }
    }

    float density(int a, int b)
    {
        Vector3 diff = particles[a].position - particles[b].position;
        float diff_square = Vector3.Dot(diff, diff);
        if(diff_square > radius2) return 0;
        return ((4 * mass) / (pi * radius8)) * Mathf.Pow(radius2 - diff_square, 3);
    }

    void compute_density()
    {
        for(int i = 0; i < n_particle; ++i)
        {
            int write_offset = neighbor_offset[i];
            float sum = 0;
            particles[i].density = (4 * mass) / (pi * radius2);
            for(int j = 0; j < particles[i].neighbor; ++j)
            {
                int neighbor_index = neighbors[write_offset + j];
                if(i == neighbor_index) continue;
                particles[i].density += density(i, neighbor_index);
            }
        }
    }

    Vector3 acceleration(int a, int b)
    {
        float rho_i = particles[a].density,
        rho_j = particles[b].density;
        Vector3 diff_pos = particles[a].position - particles[b].position,
        diff_vel = particles[a].velocity - particles[b].velocity;
        float diff_square = Vector3.Dot(diff_pos, diff_pos);
        if(diff_square > radius2) return new Vector3(0, 0, 0);
        float q = Mathf.Sqrt(diff_square) / radius,
        q2 = 1 - q;
        return mass * q2 / (pi * radius4 * rho_j)
        * (15 * bulk_modulus * (rho_i + rho_j - 2 * rest_density)
        * q2 * diff_pos / q - 40 * viscosity_coefficient * diff_vel);
    }

    void compute_acceleration()
    {
        for(int i = 0; i < n_particle; ++i)
        {
            int write_offset = neighbor_offset[i];
            particles[i].acceleration = new Vector3(0, g, 0);
            for(int j = 0; j < particles[i].neighbor; ++j)
            {
                int neighbor_index = neighbors[write_offset + j];
                if(i == neighbor_index) continue;
                particles[i].acceleration += acceleration(i, neighbor_index);
            }
            /* if(float.IsNaN(particles[i].acceleration[0]) || float.IsNaN(particles[i].acceleration[1]) || float.IsNaN(particles[i].acceleration[2])) Debug.LogFormat("particles[{0}].acceleration = {1}", i, particles[i].acceleration); */
        }
    }

    void integrate()
    {
        for(int i = 0; i < n_particle; ++i)
        {
            particle p = particles[i];
            p.velocity += p.acceleration * dt;
            p.position += p.velocity * dt;
            if(p.position.x < bound[0])
            {
                p.velocity.x *= damping;
                p.position.x = bound[0] + epsilon;
            }
            else if(p.position.x > bound[1])
            {
                p.velocity.x *= damping;
                p.position.x = bound[1] - epsilon;
            }
            if(p.position.y < bound[2])
            {
                p.velocity.y *= damping;
                p.position.y = bound[2] + epsilon;
            }
            else if(p.position.y > bound[3]) 
            {
                p.velocity.y *= damping;
                p.position.y = bound[3] - epsilon;
            }
            if(p.position.z < bound[4])
            {
                p.velocity.z *= damping;
                p.position.z = bound[4] + epsilon;
            }
            else if(p.position.z > bound[5]) 
            {
                p.velocity.z *= damping;
                p.position.z = bound[5] - epsilon;
            }
            particles[i] = p;
            /* if(float.IsNaN(particles[i].position[0]) || float.IsNaN(particles[i].position[1]) || float.IsNaN(particles[i].position[2])) Debug.LogFormat("particles[{0}].position = {1}", i, particles[i].position); */
            particles[i].cell_index = hash(get_cell(particles[i].position));
        }
    }

    void debug_integrate()
    {
        for(int i = 0; i < n_particle; ++i)
        {
            particle p = particles[i];
            Vector3 old_pos = p.position;
            p.velocity += p.acceleration * dt;
            p.position += p.velocity * dt;
            if(p.position.x < bound[0])
            {
                p.velocity.x *= damping;
                p.position.x = bound[0] + epsilon;
            }
            else if(p.position.x > bound[1])
            {
                p.velocity.x *= damping;
                p.position.x = bound[1] - epsilon;
            }
            if(p.position.y < bound[2])
            {
                p.velocity.y *= damping;
                p.position.y = bound[2] + epsilon;
            }
            else if(p.position.y > bound[3]) 
            {
                p.velocity.y *= damping;
                p.position.y = bound[3] - epsilon;
            }
            if(p.position.z < bound[4])
            {
                p.velocity.z *= damping;
                p.position.z = bound[4] + epsilon;
            }
            else if(p.position.z > bound[5]) 
            {
                p.velocity.z *= damping;
                p.position.z = bound[5] - epsilon;
            }
            particles[i] = p;
            if(float.IsNaN(particles[i].position[0]) || float.IsNaN(particles[i].position[1]) || float.IsNaN(particles[i].position[2]))
            {
                Debug.LogFormat("particles[{0}].position = {1}, old_pos = {2}, velocity = {3}, acceleration = {4}",
                i, particles[i].position, old_pos, p.velocity, p.acceleration);
            }
            particles[i].cell_index = hash(get_cell(particles[i].position));
        }
    }

    void OnDrawGizmos()
    {
        Vector3 bound_center = new Vector3(0.5f * (bound[0] + bound[1]), 0.5f * (bound[2] + bound[3]), 0.5f * (bound[4] + bound[5])),
        bound_size = new Vector3(bound[1] - bound[0], bound[3] - bound[2], bound[5] - bound[4]);
        Gizmos.color = new Color(0, 0, 0, 1);
        Gizmos.DrawWireCube(bound_center, bound_size);
        for(int i = 0; i < n_particle; ++i)
        {
            Gizmos.color = new Color(0, 0, 1, 1);
            Gizmos.DrawCube(particles[i].position, new Vector3(particle_size, particle_size, particle_size));
        }
    }
}