using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class water_start : MonoBehaviour
{
    Rigidbody rigid_body;
    public GameObject waterLevel;

    void Start()
    {
        rigid_body = GetComponent<Rigidbody>();            
    }

    void Update()
    {
        if(transform.position.y > waterLevel.transform.position.y) /* when vehicle is above water */
            rigid_body.useGravity = true;
        if(transform.position.y < waterLevel.transform.position.y) /* when vehicle is under water */
            rigid_body.useGravity = false;
    }
}
