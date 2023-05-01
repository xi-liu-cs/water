using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(MeshRenderer))]
public class mesh_generator : MonoBehaviour
{
    struct Triangle {
        public Vector3 a;
        public Vector3 b;
        public Vector3 c;

        public Vector3 this [int i] {
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

    [Header("== Particle Configurations ==")]
    [Tooltip("Reference to the script that handles all particles")]
    public fluid_gpu fluid_cs;

    [Header("== Triangle Configurations ==")]
    [Tooltip("How many points will we render across the cube grid, per axis?")]
    public int n_point_per_axis = 50;
    // read-only: the calculated number of triangles outputted by the marching cubes shader
    private int numTriangles = 0;
    // read-only: the triangles outputted by the marching cubes shader
    private Triangle[] triangles = new Triangle[0];

    [Header("Debug Configurations")]
    [Tooltip("When active, prints out details about the triangles outputted by teh marching cubes shader")]
    public bool verbose_triangles = false;

    const int thread_group_size = 8;
    public density_generator density_gen;
    public ComputeShader shader;
    public Material material;
    public float isolevel;
    public float boundsSize = 1;
    public Vector3 offset = Vector3.zero;
    Mesh fluid;
    MeshRenderer fluid_mesh_renderer;
    public ComputeBuffer triangle_buffer;
    public ComputeBuffer point_buffer;
    public ComputeBuffer triangle_count_buffer;
    public ComputeBuffer voxel_density_buffer;
    public int n_point;
    public int n_voxel_per_axis;
    public int n_voxel;

    private int size_property;
    private int particle_buffer_property; 

    private void Awake() {
        // Move to center of the world
        transform.position = Vector3.zero;

        // Initialize Particles
        fluid_cs.Initialize();

        // Initialize shader properties
        size_property = Shader.PropertyToID("size");
        particle_buffer_property = Shader.PropertyToID("particle_buffer");

        // Dunno what these do -_-
        fluid_mesh_renderer = gameObject.GetComponent<MeshRenderer>();
        fluid = GetComponent<MeshFilter>().mesh;
        CreateBuffers();
        density_gen.Awake();
    }
    

    void Update() {
        // Tell our particle manager to update the particles
        fluid_cs.UpdateParticles();
        
        UpdateChunkMesh();
        UpdateMesh(fluid);
    }

    unsafe void CreateBuffers()
    {
        n_point = n_point_per_axis * n_point_per_axis * n_point_per_axis;
        n_voxel_per_axis = n_point_per_axis - 1;
        n_voxel = n_voxel_per_axis * n_voxel_per_axis * n_voxel_per_axis;
        int maxTriangleCount = n_voxel * 5;
        triangle_buffer = new ComputeBuffer(maxTriangleCount, sizeof(Triangle), ComputeBufferType.Append);
        triangle_count_buffer = new ComputeBuffer(1, sizeof (int), ComputeBufferType.Raw);
        voxel_density_buffer = new ComputeBuffer(n_point, sizeof(float));
        point_buffer = new ComputeBuffer(n_point, 3 * sizeof(float));
        shader.SetBuffer(0, "triangles", triangle_buffer);
        shader.SetInt("n_point_per_axis", n_point_per_axis);
        shader.SetFloat("isolevel", isolevel);
        shader.SetBuffer(0, "voxel_density", voxel_density_buffer);
        shader.SetBuffer(0, "particles", fluid_cs.particle_buffer);
        shader.SetBuffer(0, "points", point_buffer);
        material.SetFloat(size_property, fluid_cs.particleRenderRadius);
        material.SetBuffer(particle_buffer_property, fluid_cs.particle_buffer);
        fluid_mesh_renderer.material = material;
    }

    unsafe public void UpdateChunkMesh()
    {
        int n_voxel_per_axis = n_point_per_axis - 1;
        int numThreadsPerAxis = Mathf.CeilToInt (n_voxel_per_axis / (float) thread_group_size);
        float pointSpacing = boundsSize / (n_point_per_axis - 1);
        Vector3 coord = Vector3.zero;
        Vector3 center = - new Vector3(boundsSize, boundsSize, boundsSize) / 2 + coord * boundsSize + Vector3.one * boundsSize / 2;
        Vector3 worldBounds = new Vector3(boundsSize, boundsSize, boundsSize);
        density_gen.generate(point_buffer, n_point_per_axis, boundsSize, worldBounds, center, offset, pointSpacing);
        
        /* float[] a = new float[100];
        voxel_density_buffer.GetData(a);
        Debug.Log("voxel");
        for(int i = 0; i < 100; ++i) Debug.Log(a[i]); */

        triangle_buffer.SetCounterValue(0);
        shader.Dispatch (0, numThreadsPerAxis, numThreadsPerAxis, numThreadsPerAxis);
        ComputeBuffer.CopyCount (triangle_buffer, triangle_count_buffer, 0);
        int[] triCountArray = { 0 };
        triangle_count_buffer.GetData (triCountArray);
        numTriangles = triCountArray[0];
        triangles = new Triangle[numTriangles];
        triangle_buffer.GetData(triangles, 0, 0, numTriangles);
        if (verbose_triangles) {
            Debug.Log($"Triangle count: {numTriangles}");
            for(int i = 0; i < numTriangles; i++) {
                Debug.Log($"a: {triangles[i].a} | b: {triangles[i].b} | c: {triangles[i].c}");
            }
        }
    }

    private void UpdateMesh(Mesh mesh) {
        mesh.Clear();
        var vertices = new Vector3[numTriangles * 3];
        var meshTriangles = new int[numTriangles * 3];

        for(int i = 0; i < numTriangles; i++) {
            for(int j = 0; j < 3; j++) {
                meshTriangles[i * 3 + j] = i * 3 + j; /* assign index */
                vertices[i * 3 + j] = triangles[i][j];
            }
        }
        mesh.vertices = vertices;
        mesh.triangles = meshTriangles;
        /* Color[] colors = new Color[vertices.Length];
        for(int i = 0; i < vertices.Length; ++i)
            colors[i] = Color.Lerp(Color.red, Color.green, vertices[i].y);
        mesh.colors = colors; */
        mesh.RecalculateNormals();
        /* Graphics.DrawMeshInstancedIndirect(mesh, 0, material, new Bounds(Vector3.zero, new Vector3(1000f, 1000f, 1000f)), fluid_cs.arg_buffer, castShadows: UnityEngine.Rendering.ShadowCastingMode.Off); */
    }

    void OnDestroy()
    {
        //particle_buffer.Release();
        voxel_density_buffer.Release();
        if(triangle_buffer != null) {
            triangle_buffer.Release();
            triangle_count_buffer.Release();
        }
    }
}