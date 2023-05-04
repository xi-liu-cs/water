using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

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

    //ComputeBuffer temp_particle_buffer, temp_buffer, temp_pos_buffer;

    public struct particle
    {
        public Vector3 position;
    }

    public void Awake()
    {
        n_voxel = mesh_gen.n_voxel;
        find_kernel();
        voxel_density_buffer = mesh_gen.voxel_density_buffer;
        particle_buffer = fluid_cs.particle_buffer;
        cube_corner_neighbor_list_buffer = new ComputeBuffer(mesh_gen.n_point * max_particles_per_cube, sizeof(int));
        cube_corner_neighbor_tracker_buffer = new ComputeBuffer(mesh_gen.n_point, sizeof(int));

        //temp_particle_buffer = new ComputeBuffer(fluid_cs.numParticles, sizeof(float)*3);
        //temp_buffer = new ComputeBuffer(fluid_cs.numParticles, sizeof(int)*3);
        //temp_pos_buffer = new ComputeBuffer(fluid_cs.numParticles, sizeof(float)*3);

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
        //int n_point = n_point_per_axis * n_point_per_axis * n_point_per_axis;
        int numThreadsPerAxis = Mathf.CeilToInt (mesh_gen.n_point_per_axis / (float) thread_group_size);

        float[] pointGridCellSizes = new float[3] {
            fluid_cs.bounds[0] / (mesh_gen.n_point_per_axis-1),
            fluid_cs.bounds[1] / (mesh_gen.n_point_per_axis-1),
            fluid_cs.bounds[2] / (mesh_gen.n_point_per_axis-1)
        };
        Debug.Log($"Calculated Grid Cell Sizes: ({pointGridCellSizes[0]},{pointGridCellSizes[1]},{pointGridCellSizes[2]})");
        int[] dimensions = new int[3] {
            mesh_gen.n_point_per_axis,
            mesh_gen.n_point_per_axis,
            mesh_gen.n_point_per_axis
        };
        Debug.Log($"Calculated Dimensions: ({dimensions[0]},{dimensions[1]},{dimensions[2]})");

        /* particle_buffer.SetData(particles); */

        float r = pointGridCellSizes[0] * 0.75f;
        
        densityShader.SetInt("numParticles", fluid_cs.numParticles);
        densityShader.SetFloat("radius", r);
        densityShader.SetFloat("radius2", r*r);
        densityShader.SetFloat("radius3", r*r*r);
        densityShader.SetFloats("bounds", fluid_cs.bounds);
        densityShader.SetFloat("max_particles_per_grid", fluid_cs.numParticlesPerGridCell);
        densityShader.SetFloats("origin",fluid_cs.origin);
        densityShader.SetFloats("pointGridCellSizes", pointGridCellSizes);
        densityShader.SetInts("dimension", dimensions);
        // densityShader.SetFloat("mass", fluid_cs.particleMass);
        densityShader.SetFloat("particle_size", fluid_cs.particleRenderSize);

        densityShader.SetFloat("pi", Mathf.PI);
        densityShader.SetInt("max_particles_per_cube", max_particles_per_cube);
        densityShader.SetInt("n_point_per_axis", mesh_gen.n_point_per_axis);

        // densityShader.SetFloat("sphere_radius", sphere_radius);
        // densityShader.SetVector ("center", new Vector4(center.x, center.y, center.z));
        densityShader.SetVector ("offset", new Vector4(offset.x, offset.y, offset.z));
        //densityShader.SetFloat ("spacing", spacing);
        //densityShader.SetVector("worldSize", worldBounds);
        // densityShader.SetFloat("grid_size", fluid_cs.gridCellSize);

        densityShader.SetBuffer(clear_cube_corner_neighbor_tracker_kernel, "cube_corner_neighbor_tracker", cube_corner_neighbor_tracker_buffer);

        densityShader.SetBuffer(compute_neighbor_list_kernel, "particles", fluid_cs.particle_buffer);
        //densityShader.SetBuffer(compute_neighbor_list_kernel, "bound", fluid_cs.bound_buffer);
        densityShader.SetBuffer(compute_neighbor_list_kernel, "cube_corner_neighbor_list", cube_corner_neighbor_list_buffer);
        densityShader.SetBuffer(compute_neighbor_list_kernel, "cube_corner_neighbor_tracker", cube_corner_neighbor_tracker_buffer);

        densityShader.SetBuffer(compute_density_kernel, "particles", fluid_cs.particle_buffer);
        densityShader.SetBuffer(compute_density_kernel, "voxel_density", mesh_gen.voxel_density_buffer);
        densityShader.SetBuffer(compute_density_kernel, "points", point_buffer);
        //densityShader.SetBuffer(compute_density_kernel, "neighbor_list", fluid_cs.neighbor_list_buffer);
        //densityShader.SetBuffer(compute_density_kernel, "neighbor_tracker", fluid_cs.neighbor_tracker_buffer);
        //densityShader.SetBuffer(compute_density_kernel, "bound", fluid_cs.bound_buffer);
        densityShader.SetBuffer(compute_density_kernel, "cube_corner_neighbor_list", cube_corner_neighbor_list_buffer);
        densityShader.SetBuffer(compute_density_kernel, "cube_corner_neighbor_tracker", cube_corner_neighbor_tracker_buffer);

        densityShader.Dispatch(clear_cube_corner_neighbor_tracker_kernel, numThreadsPerAxis, numThreadsPerAxis, numThreadsPerAxis);
        densityShader.Dispatch(compute_neighbor_list_kernel, Mathf.CeilToInt((float)fluid_cs.numParticles / 1024f), 1,1);
        densityShader.Dispatch(compute_density_kernel, numThreadsPerAxis, numThreadsPerAxis, numThreadsPerAxis);

        /*
        int[] test = new int[mesh_gen.n_point];
        cube_corner_neighbor_tracker_buffer.GetData(test);
        string top = "";
        string bottom = "";
        for(int i = 0; i < 1000; i++) {
            //if (a[i] == null) {
                top += $"{i}\t|";
                bottom += $"{test[i]}\t|";
            //}
        }
        Debug.Log($"neighbor tracker counter:\n{top}\n{bottom}");
        float3[] test_pos = new float3[fluid_cs.numParticles];
        temp_pos_buffer.GetData(test_pos);
        top = "";
        bottom = "";
        for(int i = 0; i < 1000; i++) {
            //if (a[i] == null) {
                top += $"{i}\t|";
                bottom += $"{test_pos[i]}\t|";
            //}
        }
        Debug.Log($"World Positions Per particle:\n{top}\n{bottom}");
        */

        if (buffersToRelease != null) {
            foreach (var b in buffersToRelease) {
                b.Release();
            }
        }
        return point_buffer;
    }

    void OnDestroy() {
        //temp_buffer.Dispose();
        //temp_pos_buffer.Dispose();
        cube_corner_neighbor_tracker_buffer.Dispose();
    }
}