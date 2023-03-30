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
    particle_size = 2f,
    radius = 1f,
    gas_constant = 2000f,
    dt = 0.0008f,
    rest_density = 1f,
    damping = -0.5f;
    public Vector3 offset,
    g = new Vector3(0f, -9.81f, 0f);

    [Header("fluid")]
    public int n_particle = 50000,
    dimension = 100,
    maximum_particle_per_cell = 500;

    [StructLayout(LayoutKind.Sequential, Size = 28)]
    public struct particle
    {
        public Vector3 position;
        public Vector4 color;
    }

    particle[] particles;

    int size_property = Shader.PropertyToID("size"),
    particle_buffer_property = Shader.PropertyToID("particle_buffer");

    public ComputeShader compute_shader;

    ComputeBuffer particle_buffer,
    arg_buffer;

    void Awake()
    {
        malloc_particle();
        find_kernel();
        compute_shader_init();
        compute_buffer_init();
    }

    void Update()
    {
        material.SetFloat(size_property, particle_size);
        material.SetBuffer(particle_buffer_property, particle_buffer);
        Graphics.DrawMeshInstancedIndirect(particle_mesh, 0, material, new Bounds(Vector3.zero, new Vector3(100f, 100f, 100f)), arg_buffer, castShadows: UnityEngine.Rendering.ShadowCastingMode.Off);
    }

    void malloc_particle()
    {
        particles = new particle[n_particle];
        int particle_per_dimension = Mathf.CeilToInt(Mathf.Pow(n_particle, 1f / 3f)),
        i = 0;
        while(i < n_particle)
        {
            for(int x = 0; x < particle_per_dimension; ++x)
                for(int y = 0; y < particle_per_dimension; ++y)
                    for(int z = 0; z < particle_per_dimension; ++z)
                    {
                        float r = Random.Range(0f, 5f);
                        Vector3 pos = new Vector3(1.5f * r * x, r * y, r * z) + offset;
                        particles[i] = new particle
                        {
                            position = pos,
                            color = new Vector4(0.3f, 0.5f, 1f, 0.5f)
                        };
                        if(++i == n_particle) return;
                    }
        }
    }

    void find_kernel()
    {

    }

    void compute_shader_init()
    {

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
        particle_buffer = new ComputeBuffer(n_particle, sizeof(particle));
        particle_buffer.SetData(particles);
        arg_buffer = new ComputeBuffer(1, arg.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        arg_buffer.SetData(arg);
    }

    void OnDestroy()
    {
        free();
    }

    void free()
    {
        particle_buffer.Dispose();
    }
}