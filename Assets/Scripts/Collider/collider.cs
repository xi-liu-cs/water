using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class collider : MonoBehaviour
{
    void Start()
    {
        foreach(Transform childObject in transform)
        {
            MeshFilter filter = childObject.gameObject.GetComponent<MeshFilter>();
            if(!filter)
                filter = childObject.gameObject.AddComponent<MeshFilter>();
            Mesh mesh = filter.mesh;
            if(mesh != null)
            {                      
                MeshCollider meshCollider = childObject.gameObject.AddComponent<MeshCollider>();
                meshCollider.sharedMesh = mesh;
            }
        }    
    }
}
