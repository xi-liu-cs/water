using System.Collections;
using System.Collections.Generic;
using UnityEngine;

class density_generator : MonoBehaviour
{
    public ComputeShader density_shader;
    int thread_group_size = 8;
    public List<ComputeBuffer> buffer_release;
    public int density_kernel;

    void Awake()
    {
        find_kernel();
    }

    void find_kernel()
    {
        density_kernel = density_shader.FindKernel("density");
    }

    void OnValidate()
    {
        if(FindObjectOfType<mesh_generator>())
            FindObjectOfType<mesh_generator>().update_mesh();    
    }

    public virtual ComputeBuffer generate(ComputeBuffer point_buffer, int n_point_per_axis, float bound_size, Vector3 world_bound, Vector3 center, Vector3 offset, float spacing)
    {
        int n_point = n_point_per_axis * n_point_per_axis * n_point_per_axis,
        n_thread_per_axis = Mathf.CeilToInt(n_point_per_axis / (float)thread_group_size);
        density_shader.SetBuffer(density_kernel, "points", point_buffer);
        density_shader.SetInt("n_point_per_axis", n_point_per_axis);
        density_shader.SetFloat("bound_size", bound_size);
        density_shader.SetVector("center", new Vector4(center.x, center.y, center.z));
        density_shader.SetVector("offset", new Vector4(offset.x, offset.y, offset.z));
        density_shader.SetFloat("spacing", spacing);
        density_shader.SetVector("world_size", world_bound);
        density_shader.Dispatch(density_kernel, n_thread_per_axis, n_thread_per_axis, n_thread_per_axis);
        if(buffer_release != null)
            for(int i = 0; i < buffer_release.Count; ++i)
                buffer_release[i].Release();
        return point_buffer;
    }
}