using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[ExecuteInEditMode]
public class PlaneObstacle : MonoBehaviour
{

    public Vector3 size;
    public float radius;
    [ReadOnly] public Vector3 forward, up, right;
    public Transform normalTarget, particleTarget;
    [ReadOnly] public Vector3 normalVector, targetVector, projectionPoint;
    [ReadOnly] public float dotBetweenParticleAndNormal;

    public bool isIntersecting = false;

    void OnDrawGizmos() {
        Gizmos.color = Color.blue;
        Gizmos.DrawLine(transform.position, transform.position + normalVector);

        Gizmos.color = Color.yellow;
        Gizmos.DrawLine(transform.position, transform.position + targetVector);

        Gizmos.color = Color.red;
        Gizmos.DrawSphere(projectionPoint, 0.1f);
    }

    // Update is called once per frame
    void Update() {
        // https://forum.unity.com/threads/projection-of-point-on-plane.855958/
        targetVector = particleTarget.position - transform.position;
        normalVector = normalTarget.position - transform.position;
        dotBetweenParticleAndNormal = Vector3.Dot(targetVector, normalVector);
        projectionPoint = ClosestPointOnPlane(targetVector, normalVector, particleTarget.position);
        /*
        size = new Vector3(
            transform.lossyScale.x, 
            transform.lossyScale.y, 
            transform.lossyScale.z
        );
        if (target.hasChanged || transform.hasChanged) {
            
            up = transform.up;
            forward = -1f * transform.forward;
            right = transform.right;

            Vector3 targetPosition = target.position - transform.position;
            float d = Vector3.Dot(targetPosition, forward);

            Vector2 projection = new Vector2(
                Vector3.Dot(targetPosition, up),
                Vector3.Dot(targetPosition, right)
            );

            Vector2 start = new Vector2(size.x, size.y)/(-2f) - new Vector2(radius, radius);
            Vector2 end = new Vector2(size.x, size.y)/(2f) + new Vector2(radius, radius);

            isIntersecting = (d <= 0 && d >= -size.z && projection.x >= start.x && projection.y >= start.y && projection.x <= end.x && projection.y <= end.y);

            target.hasChanged = false;
            transform.hasChanged = false;
        }
        */
    }

    public Vector3 ClosestPointOnPlane(Vector3 planeOffset, Vector3 planeNormal, Vector3 point) {
        return point + DistanceFromPlane(planeOffset, planeNormal, point) * planeNormal;
    }

    public float DistanceFromPlane(Vector3 planeOffset, Vector3 planeNormal, Vector3 point) {
        return Vector3.Dot(planeOffset - point, planeNormal);
    }
}
