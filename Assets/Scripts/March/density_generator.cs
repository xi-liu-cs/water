using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public abstract class density_generator : MonoBehaviour
{
    const int thread_group_size = 8;
    public ComputeShader densityShader;
    public ComputeBuffer particle_buffer,
    voxel_density_buffer,
    cube_corner_neighbor_list_buffer,
    cube_corner_neighbor_tracker_buffer;
    protected List<ComputeBuffer> buffersToRelease;
    /* public particle[] particles; */
    public mesh_generator mesh_gen;
    public fluid_gpu fluid_cs;
    public float sphere_radius = 50f;

    public int max_particles_per_cube = 10,
    n_voxel,
    clear_cube_corner_neighbor_tracker_kernel,
    compute_neighbor_list_kernel,
    compute_density_kernel;

    public struct particle
    {
        public Vector3 position;
        public Vector4 color;
    }

    public void Awake()
    {
        n_voxel = mesh_gen.n_voxel;
        find_kernel();
        voxel_density_buffer = mesh_gen.voxel_density_buffer;
        particle_buffer = fluid_cs.particle_buffer;
        cube_corner_neighbor_list_buffer = new ComputeBuffer(n_voxel * max_particles_per_cube, sizeof(int));
        cube_corner_neighbor_tracker_buffer = new ComputeBuffer(n_voxel, sizeof(int));
        /* int n = mesh_gen.n_point;
        particles = new particle[n];
        for(int i = 0; i < n; ++i)
        {
            Vector3 pos = new Vector3(UnityEngine.Random.Range(0f, 10f), UnityEngine.Random.Range(0f, 10f), UnityEngine.Random.Range(0f, 10f));
            particles[i] = new particle
            {
                position = pos,
                color = new Vector4(0.3f, 0.5f, 1f, 0.5f)
            };
        } */
    }

    public void find_kernel()
    {
        clear_cube_corner_neighbor_tracker_kernel = densityShader.FindKernel("clear_cube_corner_neighbor_tracker");
        compute_neighbor_list_kernel = densityShader.FindKernel("compute_neighbor_list");
        compute_density_kernel = densityShader.FindKernel("compute_density");
    }

    public virtual ComputeBuffer generate (ComputeBuffer point_buffer, int n_point_per_axis, float boundsSize, Vector3 worldBounds, Vector3 center, Vector3 offset, float spacing) {
        int n_point = n_point_per_axis * n_point_per_axis * n_point_per_axis;
        int numThreadsPerAxis = Mathf.CeilToInt (n_point_per_axis / (float) thread_group_size);
        /* particle_buffer.SetData(particles); */
        densityShader.SetFloat("mass", fluid_cs.mass);
        densityShader.SetFloat("radius", fluid_cs.radius);
        densityShader.SetFloat("radius2", fluid_cs.radius2);
        densityShader.SetFloat("radius3", fluid_cs.radius3);
        densityShader.SetFloat("particle_size", fluid_cs.particle_size);
        densityShader.SetFloat("boundsSize", boundsSize);
        densityShader.SetFloat("sphere_radius", sphere_radius);
        densityShader.SetFloat("pi", Mathf.PI);
        densityShader.SetFloat("max_particles_per_grid", fluid_cs.max_particles_per_grid);
        densityShader.SetInt("n_point_per_axis", n_point_per_axis);
        densityShader.SetVector ("center", new Vector4(center.x, center.y, center.z));
        densityShader.SetVector ("offset", new Vector4(offset.x, offset.y, offset.z));
        densityShader.SetFloat ("spacing", spacing);
        densityShader.SetVector("worldSize", worldBounds);
        densityShader.SetFloat("grid_size", fluid_cs.grid_size);
        densityShader.SetInts("dimension", fluid_cs.dimension_array);
        densityShader.SetInt("max_particles_per_cube", max_particles_per_cube);

        densityShader.SetBuffer(clear_cube_corner_neighbor_tracker_kernel, "cube_corner_neighbor_tracker", cube_corner_neighbor_tracker_buffer);

        densityShader.SetBuffer(compute_neighbor_list_kernel, "particles", fluid_cs.particle_buffer);
        densityShader.SetBuffer(compute_neighbor_list_kernel, "bound", fluid_cs.bound_buffer);
        densityShader.SetBuffer(compute_neighbor_list_kernel, "cube_corner_neighbor_list", cube_corner_neighbor_list_buffer);
        densityShader.SetBuffer(compute_neighbor_list_kernel, "cube_corner_neighbor_tracker", cube_corner_neighbor_tracker_buffer);

        densityShader.SetBuffer(compute_density_kernel, "particles", fluid_cs.particle_buffer);
        densityShader.SetBuffer(compute_density_kernel, "voxel_density", mesh_gen.voxel_density_buffer);
        densityShader.SetBuffer(compute_density_kernel, "points", point_buffer);
        densityShader.SetBuffer(compute_density_kernel, "neighbor_list", fluid_cs.neighbor_list_buffer);
        densityShader.SetBuffer(compute_density_kernel, "neighbor_tracker", fluid_cs.neighbor_tracker_buffer);
        densityShader.SetBuffer(compute_density_kernel, "bound", fluid_cs.bound_buffer);
        densityShader.SetBuffer(compute_density_kernel, "cube_corner_neighbor_list", cube_corner_neighbor_list_buffer);
        densityShader.SetBuffer(compute_density_kernel, "cube_corner_neighbor_tracker", cube_corner_neighbor_tracker_buffer);

        densityShader.Dispatch(clear_cube_corner_neighbor_tracker_kernel, numThreadsPerAxis, numThreadsPerAxis, numThreadsPerAxis);
        densityShader.Dispatch(compute_neighbor_list_kernel, numThreadsPerAxis, numThreadsPerAxis, numThreadsPerAxis);
        densityShader.Dispatch(compute_density_kernel, numThreadsPerAxis, numThreadsPerAxis, numThreadsPerAxis);

        if (buffersToRelease != null) {
            foreach (var b in buffersToRelease) {
                b.Release();
            }
        }
        return point_buffer;
    }
}