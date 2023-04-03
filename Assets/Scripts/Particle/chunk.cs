using System.Collections;
using System.Collections.Generic;
using UnityEngine;

class chunk : MonoBehaviour
{
    public Vector3Int coordinate;
    public Mesh mesh;
    MeshFilter mesh_filter;
    MeshRenderer mesh_renderer;
    MeshCollider mesh_collider;
    bool generate_collider;

    public void destroy_or_disable()
    {
        if(Application.isPlaying)
        {
            mesh.Clear();
            gameObject.SetActive(false);
        }
        else
            DestroyImmediate(gameObject, false);
    }

    public void setup(Material mat, bool generate_collider)
    {
        this.generate_collider = generate_collider;
        mesh_filter = GetComponent<MeshFilter>();
        mesh_renderer = GetComponent<MeshRenderer>();
        mesh_collider = GetComponent<MeshCollider>();
        if(!mesh_filter)
            mesh_filter = gameObject.AddComponent<MeshFilter>();
        if(!mesh_renderer)
            mesh_renderer = gameObject.AddComponent<MeshRenderer>();
        if(!mesh_collider && generate_collider)
            mesh_collider = gameObject.AddComponent<MeshCollider>();
        if(!mesh_collider && !generate_collider)
            DestroyImmediate(mesh_collider);
        mesh = mesh_filter.sharedMesh;
        if(!mesh)
        {
            mesh = new Mesh();
            mesh.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
            mesh_filter.sharedMesh = mesh;
        }
        if(generate_collider)
        {
            if(!mesh_collider.sharedMesh)
                mesh_collider.sharedMesh = mesh;
            mesh_collider.enabled = false;
            mesh_collider.enabled = true;
        }
        mesh_renderer.material = mat;
    }
}