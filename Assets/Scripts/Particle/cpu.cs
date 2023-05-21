using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;

public class cpu : MonoBehaviour
{
    struct particle
    {
        public Vector3 position,
        velocity,
        acceleration;
        public float mass,
        density,
        pressure;
    }

    public int n_particle = 1000;
    public float[] bound = {-50, 50, -50, 50, -50, 50};
    public Vector3Int dimension;
    int dimension2,
    dimension3;
    public int[] dimension_array;
    particle[] particles;
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
        malloc_particle();
    }

    void Update()
    {
        compute_density();
        compute_acceleration();
        integrate();
    }

    void malloc_particle()
    {
        particles = new particle[n_particle];
        for(int i = 0; i < n_particle; ++i)
        {
            particles[i].position = new Vector3(0.1f * UnityEngine.Random.Range(bound[0], bound[1]), 0.1f * UnityEngine.Random.Range(bound[2], bound[3]), 0.1f * UnityEngine.Random.Range(bound[4], bound[5]));
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
            float sum = 0;
            particles[i].density = (4 * mass) / (pi * radius2);
            for(int j = 0; j < n_particle; ++j)
            {
                if(i == j) continue;
                particles[i].density += density(i, j);
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
            particles[i].acceleration = new Vector3(0, g, 0);
            for(int j = 0; j < n_particle; ++j)
            {
                if(i == j) continue;
                particles[i].acceleration += acceleration(i, j);
            }
            /* Debug.LogFormat("particles[{0}].acceleration = {1}", i, particles[i].acceleration); */
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