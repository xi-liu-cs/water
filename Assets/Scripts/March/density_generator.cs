using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public abstract class density_generator : MonoBehaviour
{
    [Header("== References ==")]
    public ComputeShader densityShader;
    public mesh_generator mesh_gen;
    public fluid_gpu fluid_cs;

    [Header("== Density Configurations ==")]
    public float sphere_radius = 50f;
    public int max_particles_per_cube = 10;

    // == Kernel Indices ==
    private int clear_cube_corner_neighbor_tracker_kernel;
    private int compute_neighbor_list_kernel;
    private int compute_density_kernel;

    const int thread_group_size = 8;
    public ComputeBuffer cube_corner_neighbor_list_buffer;
    public ComputeBuffer cube_corner_neighbor_tracker_buffer;
    protected List<ComputeBuffer> buffersToRelease;

    private bool initialized = false;

    public void Initialize() {
        InitializeKernels();
        cube_corner_neighbor_list_buffer = new ComputeBuffer(mesh_gen.numVoxels * max_particles_per_cube, sizeof(int));
        cube_corner_neighbor_tracker_buffer = new ComputeBuffer(mesh_gen.numVoxels, sizeof(int));
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
        initialized = true;
    }

    public void InitializeKernels() {
        clear_cube_corner_neighbor_tracker_kernel = densityShader.FindKernel("clear_cube_corner_neighbor_tracker");
        compute_neighbor_list_kernel = densityShader.FindKernel("compute_neighbor_list");
        compute_density_kernel = densityShader.FindKernel("compute_density");
    }

    public virtual ComputeBuffer generate (ComputeBuffer point_buffer, int n_point_per_axis, float boundsSize, Vector3 worldBounds, Vector3 center, Vector3 offset, float spacing) {
        int n_point = n_point_per_axis * n_point_per_axis * n_point_per_axis;
        int numThreadsPerAxis = Mathf.CeilToInt (n_point_per_axis / (float) thread_group_size);

        // World constants
        densityShader.SetFloat("pi", Mathf.PI);

        // Retrieved from `fluid_gpu.cs` 
        densityShader.SetFloat("gridCellSize", fluid_cs.gridCellSize);
        densityShader.SetFloat("radius", fluid_cs.smoothingRadius);
        densityShader.SetFloat("radius2", fluid_cs.radius2);
        densityShader.SetFloat("radius3", fluid_cs.radius3);
        densityShader.SetFloat("particleRenderRadius", fluid_cs.particleRenderRadius);
        densityShader.SetFloats("bounds", fluid_cs.bounds);

        densityShader.SetFloat("boundsSize", boundsSize);
        densityShader.SetFloat("sphere_radius", sphere_radius);
        densityShader.SetInt("n_point_per_axis", n_point_per_axis);
        densityShader.SetFloat ("spacing", spacing);
        densityShader.SetVector("worldSize", worldBounds);
        densityShader.SetInt("max_particles_per_cube", max_particles_per_cube);
        
        densityShader.SetBuffer(clear_cube_corner_neighbor_tracker_kernel, "cube_corner_neighbor_tracker", cube_corner_neighbor_tracker_buffer);

        densityShader.SetBuffer(compute_neighbor_list_kernel, "particles", fluid_cs.particle_buffer);
        densityShader.SetBuffer(compute_neighbor_list_kernel, "cube_corner_neighbor_list", cube_corner_neighbor_list_buffer);
        densityShader.SetBuffer(compute_neighbor_list_kernel, "cube_corner_neighbor_tracker", cube_corner_neighbor_tracker_buffer);

        densityShader.SetBuffer(compute_density_kernel, "particles", fluid_cs.particle_buffer);
        densityShader.SetBuffer(compute_density_kernel, "voxel_density", mesh_gen.voxel_density_buffer);
        densityShader.SetBuffer(compute_density_kernel, "points", point_buffer);
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