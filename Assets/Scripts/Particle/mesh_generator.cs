using System.Collections;
using System.Collections.Generic;
using UnityEngine;

class mesh_generator : fluid
{
    int thread_group_size = 8;
    public density_generator density_gen;
    public bool fixed_map_size;
    public Vector3Int n_chunk = Vector3Int.one;
    public Transform viewer;
    public float view_distance = 30;
    public bool auto_update_in_editor = true,
    auto_update_in_game = true;
    public ComputeShader shader;
    public Material mat;
    public bool generate_collider;

    [Header("voxel")]
    public float isolevel,
    bound_size = 1;
    public Vector3 offset = Vector3.zero;

    [Range(2, 100)]
    public int n_point_per_axis = 30;

    GameObject chunk_holder;
    string chunk_holder_name = "chunk holder";
    List<chunk> chunks;
    Dictionary<Vector3Int, chunk> exist_chunk;
    Queue<chunk> recycle_chunk;

    ComputeBuffer triangle_buffer,
    point_buffer,
    triangle_count_buffer;

    bool setting_updated;

    int march_kernel;

    struct triangle
    {
        public Vector3 vertex_a,
        vertex_b,
        vertex_c;
        public Vector3 this [int i]
        {
            get
            {
                switch(i)
                {
                    case 0: return vertex_a;
                    case 1: return vertex_b;
                    default: return vertex_c;
                }
            }
        }
    };

    void Awake()
    {
        if(Application.isPlaying && !fixed_map_size)
        {
            malloc_chunk();
            chunk[] old_chunk = FindObjectsOfType<chunk>();
            for(int i = old_chunk.Length - 1; i >= 0; --i)
                Destroy(old_chunk[i].gameObject);
        }    
    }

    void Update()
    {
        if(Application.isPlaying && !fixed_map_size)
            run();
        if(setting_updated)
        {
            update_mesh();
            setting_updated = false;
        }    
    }

    void find_kernel_mesh_generator()
    {
        march_kernel = shader.FindKernel("march");
    }

    public void run()
    {
        create_buffer();
        if(fixed_map_size)
        {
            init_chunk();
            update_all_chunk();
        }
        else
        {
            if(Application.isPlaying)
                init_visible_chunk();
        }
        if(!Application.isPlaying)
            release_buffer();
    }

    public void update_mesh()
    {
        if((Application.isPlaying && auto_update_in_game) || (!Application.isPlaying && auto_update_in_editor))
            run();
    }

    void malloc_chunk()
    {
        recycle_chunk = new Queue<chunk>();
        chunks = new List<chunk>();
        exist_chunk = new Dictionary<Vector3Int, chunk>();
    }

    void init_visible_chunk()
    {
        if(chunks == null) return;
        create_chunk_holder();
        Vector3 viewer_position = viewer.position,
        position_div_size = viewer_position / bound_size;
        Vector3Int viewer_coordinate = new Vector3Int(Mathf.RoundToInt(position_div_size.x), Mathf.RoundToInt(position_div_size.y), Mathf.RoundToInt(position_div_size.z));
        int max_chunk_in_view = Mathf.CeilToInt(view_distance / bound_size);
        float view_distance_square = view_distance * view_distance;
        for(int i = chunks.Count - 1; i >= 0; --i)
        {
            chunk chunk_i = chunks[i];
            Vector3 center = center_from_coordinate(chunk_i.coordinate),
            viewer_offset = viewer_position - center,
            o = new Vector3(Mathf.Abs(viewer_offset.x), Mathf.Abs(viewer_offset.y), Mathf.Abs(viewer_offset.z));
            float distance_square = new Vector3(Mathf.Max(o.x, 0), Mathf.Max(o.y, 0), Mathf.Max(o.z, 0)).sqrMagnitude;
            if(distance_square > view_distance_square)
            {
                exist_chunk.Remove(chunk_i.coordinate);
                recycle_chunk.Enqueue(chunk_i);
                chunks.RemoveAt(i);
            }
        }
        for(int x = -max_chunk_in_view; x <= max_chunk_in_view; ++x)
        {
            for(int y = -max_chunk_in_view; y <= max_chunk_in_view; ++y)
            {
                for(int z = -max_chunk_in_view; z <= max_chunk_in_view; ++z)
                {
                    Vector3Int coordinate = new Vector3Int(x, y, z) + viewer_coordinate;
                    if(exist_chunk.ContainsKey(coordinate)) continue;
                    Vector3 center = center_from_coordinate(coordinate),
                    viewer_offset = viewer_position - center,
                    o = new Vector3(Mathf.Abs(viewer_offset.x), Mathf.Abs(viewer_offset.y), Mathf.Abs(viewer_offset.z));
                    float distance_square = new Vector3(Mathf.Max(o.x, 0), Mathf.Max(o.y, 0), Mathf.Max(o.z, 0)).sqrMagnitude;

                    if(distance_square <= view_distance_square)
                    {
                        Bounds bounds = new Bounds(center_from_coordinate(coordinate), bound_size * Vector3.one);
                        if(is_visible_from(bounds, Camera.main))
                        {
                            if(recycle_chunk.Count > 0)
                            {
                                chunk a = recycle_chunk.Dequeue();
                                a.coordinate = coordinate;
                                exist_chunk.Add(coordinate, a);
                                chunks.Add(a);
                                update_chunk_mesh(a);
                            }
                            else
                            {
                                chunk a = create_chunk(coordinate);
                                a.coordinate = coordinate;
                                a.setup(mat, generate_collider);
                                exist_chunk.Add(coordinate, a);
                                chunks.Add(a);
                                update_chunk_mesh(a);
                            }
                        }
                    }
                }
            }
        }
    }

    public bool is_visible_from(Bounds bounds, Camera camera)
    {
        Plane[] planes = GeometryUtility.CalculateFrustumPlanes(camera);
        return GeometryUtility.TestPlanesAABB(planes, bounds);
    }

    public void update_chunk_mesh(chunk a)
    {
        int n_voxel_per_axis = n_point_per_axis - 1,
        n_thread_per_axis = Mathf.CeilToInt(n_voxel_per_axis / (float)thread_group_size);
        float point_spacing = bound_size / (n_point_per_axis - 1);
        Vector3Int coordinate = a.coordinate;
        Vector3 center = center_from_coordinate(coordinate);
        Vector3 world_bound = bound_size * new Vector3(n_chunk.x, n_chunk.y, n_chunk.z);
        density_gen.generate(point_buffer, n_point_per_axis, bound_size, world_bound, center, offset, point_spacing);
        triangle_buffer.SetCounterValue(0);
        shader.SetBuffer(march_kernel, "points", point_buffer);
        shader.SetBuffer(march_kernel, "triangles", triangle_buffer);
        shader.SetInt("n_point_per_axis", n_point_per_axis);
        shader.SetFloat("isolevel", isolevel);
        shader.Dispatch(march_kernel, n_thread_per_axis, n_thread_per_axis, n_thread_per_axis);
        
        ComputeBuffer.CopyCount(triangle_buffer, triangle_count_buffer, 0);
        int[] triangle_count_array = {0};
        triangle_count_buffer.GetData(triangle_count_array);
        int n_triangle = triangle_count_array[0];
        
        triangle[] triangles = new triangle[n_triangle];
        triangle_buffer.GetData(triangles, 0, 0, n_triangle);

        Mesh mesh = a.mesh;
        mesh.Clear();
        Vector3[] vertices = new Vector3[3 * n_triangle];
        int[] mesh_triangles = new int[3 * n_triangle];
        for(int i = 0; i < n_triangle; ++i)
        {
            for(int j = 0; j < n_triangle; ++j)
            {
                mesh_triangles[i * 3 + j] = i * 3 + j;
                vertices[i * 3 + j] = triangles[i][j];
            }
        }
        mesh.vertices = vertices;
        mesh.triangles = mesh_triangles;
        mesh.RecalculateNormals();
    }

    public void update_all_chunk()
    {
        for(int i = 0; i < chunks.Count; ++i)
            update_chunk_mesh(chunks[i]);
    }

    void OnDestroy()
    {
        if(Application.isPlaying)
            release_buffer();
    }

    unsafe void create_buffer()
    {
        int n_point = n_point_per_axis * n_point_per_axis * n_point_per_axis,
        n_voxel_per_axis = n_point_per_axis - 1,
        n_voxel = n_voxel_per_axis * n_voxel_per_axis * n_voxel_per_axis,
        max_triangle_count = 5 * n_voxel;
        if(!Application.isPlaying || (point_buffer == null || n_point != point_buffer.count))
        {
            if(Application.isPlaying)
                release_buffer();
            compute_buffer_init();
            triangle_buffer = new ComputeBuffer(max_triangle_count, sizeof(triangle), ComputeBufferType.Append);
            point_buffer = new ComputeBuffer(n_point, 4 * sizeof(float));
            triangle_count_buffer = new ComputeBuffer(1, sizeof(int), ComputeBufferType.Raw);
        }
    }

    void release_buffer()
    {
        if(triangle_buffer != null)
        { 
            triangle_buffer.Release();
            point_buffer.Release();
            triangle_count_buffer.Release();
        }
    }

    Vector3 center_from_coordinate(Vector3Int coordinate)
    {
        if(fixed_map_size)
        {
            Vector3 total_bound = (Vector3)n_chunk * bound_size;
            return -0.25f * total_bound + (Vector3)coordinate * bound_size + Vector3.one * bound_size; /* return -0.5 * total_bound + coordinate * bound_size + Vector3.one * 0.5 * bound_size; */
        }
        return bound_size * new Vector3(coordinate.x, coordinate.y, coordinate.z);
    }

    void create_chunk_holder()
    {
        if(chunk_holder == null)
        {
            GameObject find = GameObject.Find(chunk_holder_name);
            if(find)
                chunk_holder = find;
            else
                chunk_holder = new GameObject(chunk_holder_name);
        }
    }

    void init_chunk()
    {
        create_chunk_holder();
        chunks = new List<chunk>();
        List<chunk> old_chunk = new List<chunk>(FindObjectsOfType<chunk>());
        for(int x = 0; x < n_chunk.x; ++x)
        {
            for(int y = 0; y < n_chunk.y; ++y)
            {
                for(int z = 0; z < n_chunk.z; ++z)
                {
                    Vector3Int coordinate = new Vector3Int(x, y, z);
                    bool chunk_exist = false;
                    for(int i = 0; i < old_chunk.Count; ++i)
                    {
                        if(old_chunk[i].coordinate == coordinate)
                        {
                            chunks.Add(old_chunk[i]);
                            old_chunk.RemoveAt(i);
                            chunk_exist = true;
                            break;
                        }
                    }
                    if(!chunk_exist)
                    {
                        chunk new_chunk = create_chunk(coordinate);
                        chunks.Add(new_chunk);
                    }
                    chunks[chunks.Count - 1].setup(mat, generate_collider);
                }
            }
        }        
    }

    chunk create_chunk(Vector3Int coordinate)
    {
        GameObject a = new GameObject($"chunk ({coordinate.x}, {coordinate.y}, {coordinate.z})");
        a.transform.parent = chunk_holder.transform;
        chunk new_chunk = a.AddComponent<chunk>();
        new_chunk.coordinate = coordinate;
        return new_chunk;
    }

    void OnValidate()
    {
        setting_updated = true;    
    }
}