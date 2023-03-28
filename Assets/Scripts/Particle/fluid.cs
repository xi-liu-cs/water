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
    viscosity = 2.5f,
    particle_size = 4f,
    radius = 1f,
    gas_constant = 2000f,
    dt = 0.0008f,
    rest_density = 1f,
    damping = -0.5f;
    public Vector3 g = new Vector3(0f, -9.81f, 0f);

    [Header("fluid")]
    public int n_particle = 1000,
    dimension = 10,
    maximum_particle_per_cell = 500;

    [StructLayout(LayoutKind.Sequential, Size = 28)]
    struct particle
    {
        Vector3 position;
        Vector4 color;
    }

    public ComputeShader compute_shader;

    ComputeBuffer particles,
    arg_buf,
    density,
    pressure,
    force,
    velocity;

    int malloc_particle_kernel;

    void Awake()
    {
        find_kernel();
        compute_shader_init();
        compute_buffer_init();
    }

    void Update()
    {
        compute_shader.Dispatch(malloc_particle_kernel, 1, 1, 1);
        Graphics.DrawMeshInstancedIndirect(particle_mesh, 0, material, new Bounds(Vector3.zero, new Vector3(100f, 100f, 100f)), arg_buf);
    }

    void find_kernel()
    {
        malloc_particle_kernel = compute_shader.FindKernel("malloc_particle");
    }

    void compute_shader_init()
    {
        compute_shader.SetInt("n_particle", n_particle);
        compute_shader.SetInt("dimension", dimension);
    }

    unsafe void compute_buffer_init()
    {
        uint[] arg =
        {
            particle_mesh.GetIndexCount(0),
            (uint)n_particle,
            particle_mesh.GetIndexStart(0),
            particle_mesh.GetBaseVertex(0),
            0
        };
        particles = new ComputeBuffer(n_particle, sizeof(particle));
        arg_buf = new ComputeBuffer(1, arg.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        arg_buf.SetData(arg);
        compute_shader.SetBuffer(malloc_particle_kernel, "particles", particles);
    }
}