using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class input : MonoBehaviour
{
    public float speed = 1f,
    sensitivity = 1f;
    public bool move = true;

    void Update()
    {
        if(Input.GetMouseButtonDown(1)) move = !move; /* right click, press esc to exit controlling the camera */
        if(!move) return;
        if(Input.GetKey(KeyCode.LeftArrow) || Input.GetKey(KeyCode.A))
        {
            transform.Translate(new Vector3(-speed * Time.deltaTime, 0, 0)); /* negative since the scene is reflected, otherwise 'a' and 'd' would be switched */
        }
        if(Input.GetKey(KeyCode.RightArrow) || Input.GetKey(KeyCode.D))
        {
            transform.Translate(new Vector3(speed * Time.deltaTime, 0, 0));
        }
        if(Input.GetKey(KeyCode.UpArrow) || Input.GetKey(KeyCode.W))
        {
            transform.Translate(new Vector3(0, speed * Time.deltaTime, 0));
        }
        if(Input.GetKey(KeyCode.DownArrow) || Input.GetKey(KeyCode.S))
        {
            transform.Translate(new Vector3(0, -speed * Time.deltaTime, 0));
        }
        if(Input.GetKey(KeyCode.Alpha1) || Input.GetKey(KeyCode.Keypad1))
        {
            transform.Translate(new Vector3(0, 0, speed * Time.deltaTime));
        }
        if(Input.GetKey(KeyCode.Alpha2) || Input.GetKey(KeyCode.Keypad2))
        {
            transform.Translate(new Vector3(0, 0, -speed * Time.deltaTime));
        }

        float factor = sensitivity / 10f;
        Transform c = Camera.main.transform;
        c.Rotate(0, Input.GetAxis("Mouse X") * sensitivity, 0);
        c.Rotate(-Input.GetAxis("Mouse Y") * sensitivity, 0, 0);
        c.Rotate(0, 0, -Input.GetAxis("QandE") * 90 * Time.deltaTime * factor);
        if(Input.GetMouseButtonDown(0))
            Cursor.lockState = CursorLockMode.Locked;
    }
}