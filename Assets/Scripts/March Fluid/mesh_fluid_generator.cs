using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class mesh_fluid_generator : MonoBehaviour
{
    [Header("== Particle Configurations ==")]
    public fluid_gpu fluid_cs;

    const int thread_group_size = 8;
    /* public density_generator density_gen; */
    public ComputeShader shader,
    march_shader;
    public Material mat;
    public float isolevel;
    public float boundsSize = 1;
    public Vector3 offset = Vector3.zero;
    private int n_point_per_axis;
    Mesh fluid;
    public ComputeBuffer triangle_buffer;
    public ComputeBuffer point_buffer;
    public ComputeBuffer triangle_count_buffer;
    public ComputeBuffer voxel_density_buffer;
    public int n_point;
    public int march_kernel;

    

    void Awake() {
        // We initialize our fluid_cs script. This ensures that we have everything set for our particles
        fluid_cs.Initialize();

        n_point_per_axis = fluid_cs.n_point_per_axis;
        voxel_density_buffer = fluid_cs.density_buffer;
        CreateBuffers();
        fluid = GetComponent<MeshFilter>().mesh;
        /* density_gen.Awake(); */
        find_kernel();
    }

    void Update()
    {
        fluid_cs.UpdateParticles();
        UpdateChunkMesh(fluid);
    }

    void find_kernel()
    {
        march_kernel = march_shader.FindKernel("march");
    }

    unsafe void CreateBuffers()
    {
        n_point = n_point_per_axis * n_point_per_axis * n_point_per_axis;
        int n_voxel_per_axis = n_point_per_axis - 1,
        n_voxel = n_voxel_per_axis * n_voxel_per_axis * n_voxel_per_axis,
        maxTriangleCount = n_voxel * 5;
        triangle_buffer = new ComputeBuffer(maxTriangleCount, sizeof(tri), ComputeBufferType.Append);
        triangle_count_buffer = new ComputeBuffer(1, sizeof (int), ComputeBufferType.Raw);
        march_shader.SetBuffer(march_kernel, "triangles", triangle_buffer);
        march_shader.SetInt("n_point_per_axis", n_point_per_axis);
        march_shader.SetFloat("isolevel", isolevel);
        march_shader.SetBuffer(march_kernel, "voxel_density", voxel_density_buffer);
        march_shader.SetBuffer(march_kernel, "particles", fluid_cs.particle_buffer);
    }

    public void UpdateChunkMesh(Mesh mesh)
    {
        int n_voxel_per_axis = n_point_per_axis - 1;
        int numThreadsPerAxis = Mathf.CeilToInt (n_voxel_per_axis / (float) thread_group_size);
        float pointSpacing = boundsSize / (n_point_per_axis - 1);
        Vector3 coord = Vector3.zero;
        Vector3 center = - new Vector3(boundsSize, boundsSize, boundsSize) / 2 + coord * boundsSize + Vector3.one * boundsSize / 2;
        Vector3 worldBounds = new Vector3(boundsSize, boundsSize, boundsSize);
        /* density_gen.generate(point_buffer, n_point_per_axis, boundsSize, worldBounds, center, offset, pointSpacing); */
        
        /* float[] a = new float[100];
        voxel_density_buffer.GetData(a);
        Debug.Log("voxel");
        for(int i = 0; i < 100; ++i) Debug.Log(a[i]); */

        triangle_buffer.SetCounterValue (0);
        shader.Dispatch (march_kernel, numThreadsPerAxis, numThreadsPerAxis, numThreadsPerAxis);
        ComputeBuffer.CopyCount(triangle_buffer, triangle_count_buffer, 0);
        int[] triCountArray = { 0 };
        triangle_count_buffer.GetData (triCountArray);
        int numTris = triCountArray[0];
        tri[] tris = new tri[numTris];
        triangle_buffer.GetData (tris, 0, 0, numTris);
        /* Debug.Log("triangle");
        for(int i = 0; i < numTris; ++i)
        {
            Debug.Log(tris[i].a);
            Debug.Log(tris[i].b);
            Debug.Log(tris[i].c);
        } */
        Debug.Log("triangle count = " + numTris);
        for(int i = 0; i < numTris; ++i)
        {
            Debug.Log(tris[i].a);
            Debug.Log(tris[i].b);
            Debug.Log(tris[i].c);
        }
        if(numTris == 0)
        {
            float[] a = new float[100];
            voxel_density_buffer.GetData(a);
            Debug.Log("voxel");
            for(int i = 0; i < 100; ++i) Debug.Log(a[i]);
        }

        mesh.Clear();
        var vertices = new Vector3[numTris * 3];
        var meshTriangles = new int[numTris * 3];

        for (int i = 0; i < numTris; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                meshTriangles[i * 3 + j] = i * 3 + j; /* assign index */
                vertices[i * 3 + j] = tris[i][j];
            }
        }
        mesh.vertices = vertices;
        mesh.triangles = meshTriangles;
        mesh.RecalculateNormals ();
    }

    void OnDestroy()
    {
        voxel_density_buffer.Release();
        if(triangle_buffer != null)
        {
            triangle_buffer.Release();
            triangle_count_buffer.Release();
        }
    }

    struct tri
    {
        public Vector3 a;
        public Vector3 b;
        public Vector3 c;

        public Vector3 this [int i]
        {
            get {
                switch (i) {
                    case 0:
                        return a;
                    case 1:
                        return b;
                    default:
                        return c;
                }
            }
        }
    }
}