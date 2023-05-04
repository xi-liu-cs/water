using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class mesh_generator : MonoBehaviour
{
    const int thread_group_size = 8;
    public density_generator density_gen;
    public ComputeShader shader;
    public fluid_gpu fluid_cs;
    public Material material;
    public float isolevel;
    public float boundsSize = 1;
    public Vector3 offset = Vector3.zero;
    public int n_point_per_axis = 84;
    Mesh fluid;
    MeshFilter fluid_mesh_filter;
    MeshRenderer fluid_mesh_renderer;
    public ComputeBuffer triangle_buffer,
    point_buffer,
    triangle_count_buffer,
    voxel_density_buffer,
    particle_buffer;
    public int n_point,
    n_particle,
    n_voxel_per_axis,
    n_voxel;

    int size_property = Shader.PropertyToID("size"),
    particle_buffer_property = Shader.PropertyToID("particle_buffer");

    public Color gridColor = Color.yellow;
    
    public struct particle
    {
        public Vector3 position;
    }

    void OnDrawGizmos() {
        /*
        if(Application.isPlaying) {

            Vector3 pointGridCellSizes = new Vector3(
                fluid_cs.bounds[0] / (n_point_per_axis-1),
                fluid_cs.bounds[1] / (n_point_per_axis-1),
                fluid_cs.bounds[2] / (n_point_per_axis-1)
            );
            Debug.Log($"Point grid cell sizes: ({pointGridCellSizes[0]},{pointGridCellSizes[1]},{pointGridCellSizes[2]})");
            Vector3 halfMidpoint = new Vector3(pointGridCellSizes.x/2f, pointGridCellSizes.y/2f, pointGridCellSizes.z/2f);
            int[] temp_numbers = new int[n_point];
            density_gen.cube_corner_neighbor_tracker_buffer.GetData(temp_numbers);
            Gizmos.color = gridColor;
            for(int x = 0; x < n_point_per_axis; x++) {
                for(int y = 0; y < n_point_per_axis; y++) {
                    for(int z = 0; z < n_point_per_axis; z++) {
                        int proj_index = GetProjectedGridIndexFromXYZ(x,y,z);
                        float num = Mathf.Clamp(Mathf.FloorToInt(temp_numbers[proj_index]/1f),0f,1f);
                        //Debug.Log($"Gizmos: ({x},{y},{z}) => {proj_index} => {density}");
                        Gizmos.color = new Color(1f,1f,0f,num);
                        Vector3 pos = new Vector3(
                            fluid_cs.originVector3.x - (fluid_cs.bounds[0]/2f) + (x * pointGridCellSizes.x),
                            fluid_cs.originVector3.y - (fluid_cs.bounds[1]/2f) + (y * pointGridCellSizes.y),
                            fluid_cs.originVector3.z - (fluid_cs.bounds[2]/2f) + (z * pointGridCellSizes.z)
                        );
                        Gizmos.DrawSphere(pos, 1f);
                    }
                }
            }
            

            particle[] debug_particles = new particle[n_particle];
            particle_buffer.GetData(debug_particles);
            for(int i = 0; i < debug_particles.Length; ++i)
                Gizmos.DrawSphere(debug_particles[i].position, 1f);
        }
        */
        Gizmos.color = Color.white;
        Gizmos.DrawWireCube(fluid_cs.originVector3, new Vector3(fluid_cs.bounds[0], fluid_cs.bounds[1], fluid_cs.bounds[2]));
        Gizmos.color = Color.red;
        Gizmos.DrawSphere(fluid_cs.originVector3,0.1f);
    }

    int GetProjectedGridIndexFromXYZ(int x, int y, int z) {
        return x + (fluid_cs.n_point_per_axis * y) + (fluid_cs.n_point_per_axis * fluid_cs.n_point_per_axis * z);
    }

    void Awake()
    {
        fluid_cs.Initialize();
        particle_buffer = fluid_cs.particle_buffer;
        n_point_per_axis = fluid_cs.n_point_per_axis;
        n_particle = fluid_cs.numParticles;
        gameObject.transform.position = new Vector3(0, 0, 0);
        fluid_mesh_filter = gameObject.GetComponent<MeshFilter>();
        if(fluid_mesh_filter == null)
        {
            gameObject.AddComponent<MeshFilter>();
            fluid_mesh_filter = gameObject.GetComponent<MeshFilter>();
        }
        fluid_mesh_renderer = gameObject.GetComponent<MeshRenderer>();
        if(fluid_mesh_renderer == null)
        {
            gameObject.AddComponent<MeshRenderer>();
            fluid_mesh_renderer = gameObject.GetComponent<MeshRenderer>();
        }
        fluid = GetComponent<MeshFilter>().mesh;
        CreateBuffers();
        density_gen.Awake();
    }
    

    void Update()
    {
        fluid_cs.UpdateParticles();
        UpdateChunkMesh(fluid);
    }

    unsafe void CreateBuffers()
    {
        n_point = n_point_per_axis * n_point_per_axis * n_point_per_axis;
        Debug.Log($"N Points: {n_point}");
        n_voxel_per_axis = n_point_per_axis - 1;
        n_voxel = n_voxel_per_axis * n_voxel_per_axis * n_voxel_per_axis;
        int maxTriangleCount = n_voxel * 5;
        triangle_buffer = new ComputeBuffer(maxTriangleCount, sizeof(tri), ComputeBufferType.Append);
        triangle_count_buffer = new ComputeBuffer(1, sizeof (int), ComputeBufferType.Raw);
        voxel_density_buffer = new ComputeBuffer(n_point, sizeof(float));
        point_buffer = new ComputeBuffer(n_point, 3 * sizeof(float));
        shader.SetBuffer(0, "triangles", triangle_buffer);
        shader.SetInt("n_point_per_axis", n_point_per_axis);
        shader.SetFloat("isolevel", isolevel);
        shader.SetBuffer(0, "voxel_density", voxel_density_buffer);
        shader.SetBuffer(0, "particles", fluid_cs.particle_buffer);
        shader.SetBuffer(0, "points", point_buffer);
        material.SetFloat(size_property, fluid_cs.particleRenderSize);
        material.SetBuffer(particle_buffer_property, fluid_cs.particle_buffer);
        fluid_mesh_renderer.material = material;
    }

    unsafe public void UpdateChunkMesh(Mesh mesh)
    {
        int n_voxel_per_axis = n_point_per_axis - 1;
        int numThreadsPerAxis = Mathf.CeilToInt (n_voxel_per_axis / (float) thread_group_size);
        float pointSpacing = boundsSize / (n_point_per_axis - 1);
        Vector3 coord = Vector3.zero;
        Vector3 center = - new Vector3(boundsSize, boundsSize, boundsSize) / 2 + coord * boundsSize + Vector3.one * boundsSize / 2;
        Vector3 worldBounds = new Vector3(boundsSize, boundsSize, boundsSize);
        density_gen.generate(point_buffer, n_point_per_axis, boundsSize, worldBounds, center, offset, pointSpacing);

        /*
        float[] a = new float[n_point];
        voxel_density_buffer.GetData(a);
        string top = "";
        string bottom = "";
        for(int i = 0; i < n_point; i++) {
            if (a[i] > 0) {
                top += $"{i}\t|";
                bottom += $"{a[i]}\t|";
            }
        }
        Debug.Log($"Voxels:\n{top}\n{bottom}");
        */

        triangle_buffer.SetCounterValue (0);
        shader.Dispatch (0, numThreadsPerAxis, numThreadsPerAxis, numThreadsPerAxis);
        ComputeBuffer.CopyCount (triangle_buffer, triangle_count_buffer, 0);
        int[] triCountArray = { 0 };
        triangle_count_buffer.GetData (triCountArray);
        int numTris = triCountArray[0];
        tri[] tris = new tri[numTris];
        triangle_buffer.GetData(tris, 0, 0, numTris);
        // Debug.Log("triangle count = " + numTris);
        // for(int i = 0; i < numTris; ++i)
        // {
        //     Debug.Log(tris[i].a);
        //     Debug.Log(tris[i].b);
        //     Debug.Log(tris[i].c);
        // }

        mesh.Clear();
        var vertices = new Vector3[numTris * 3];
        var meshTriangles = new int[numTris * 3];

        for(int i = 0; i < numTris; ++i)
        {
            for(int j = 0; j < 3; ++j)
            {
                meshTriangles[i * 3 + j] = i * 3 + j; // assign index
                vertices[i * 3 + j] = tris[i][j];
            }
        }
        mesh.vertices = vertices;
        mesh.triangles = meshTriangles;
        // Color[] colors = new Color[vertices.Length];
        //for(int i = 0; i < vertices.Length; ++i)
        //    colors[i] = Color.Lerp(Color.red, Color.green, vertices[i].y);
        // mesh.colors = colors; 
        mesh.RecalculateNormals();
        // Graphics.DrawMeshInstancedIndirect(mesh, 0, material, new Bounds(Vector3.zero, new Vector3(1000f, 1000f, 1000f)), fluid_cs.arg_buffer, castShadows: UnityEngine.Rendering.ShadowCastingMode.Off);
        
    }

    void OnDestroy()
    {
        particle_buffer.Release();
        voxel_density_buffer.Release();
        if(triangle_buffer != null)
        {
            triangle_buffer.Release();
            triangle_count_buffer.Release();
            particle_buffer.Release();
        }
    }

    struct tri
    {
        public Vector3 a;
        public Vector3 b;
        public Vector3 c;

        public Vector3 this [int i]
        {
            get
            {
                switch (i)
                {
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