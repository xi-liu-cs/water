using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[ExecuteInEditMode]
public class PlaneObstacle : MonoBehaviour
{

    public Vector3 size;
    public float radius;
    public Transform target;

    public bool isIntersecting = false;

    // Update is called once per frame
    void Update() {
        size = 10f * new Vector3(
            transform.lossyScale.z, 
            transform.lossyScale.x, 
            transform.lossyScale.y / 5
        );
        if (target.hasChanged || transform.hasChanged) {
            Debug.Log("CHANGING");
            Vector3 targetPosition = target.position - transform.position;
            float d = Vector3.Dot(targetPosition, transform.up);

            Vector2 projection = new Vector2(
                Vector3.Dot(targetPosition, transform.forward),
                Vector3.Dot(targetPosition, transform.right)
            );

            Vector2 start = new Vector2(size.x, size.y)/(-2f) - new Vector2(radius, radius);
            Vector2 end = new Vector2(size.x, size.y)/(2f) + new Vector2(radius, radius);

            isIntersecting = (d <= 0 && d >= -size.z && projection.x >= start.x && projection.y >= start.y && projection.x <= end.x && projection.y <= end.y);

            target.hasChanged = false;
            transform.hasChanged = false;
        }
    }
}
