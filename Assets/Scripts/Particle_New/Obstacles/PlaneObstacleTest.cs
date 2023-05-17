using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[ExecuteInEditMode]
public class PlaneObstacleTest : MonoBehaviour
{
    public Transform[] vertices = new Transform[3];
    public Vector3 normalVector;
    public Vector3 size;
    public float radius;
    public Transform particleTarget;
    [ReadOnly] public Vector3 centroid, targetVector, projectionPoint;
    [ReadOnly] public float dotBetweenParticleAndNormal;

    public bool isIntersecting = false;

    void OnDrawGizmos() {
        Gizmos.color = Color.white;
        Gizmos.DrawSphere(centroid, 0.05f);

        Gizmos.color = Color.green;
        Gizmos.DrawSphere(vertices[0].position, 0.05f);
        Gizmos.DrawSphere(vertices[1].position, 0.05f);
        Gizmos.DrawSphere(vertices[2].position, 0.05f);

        Gizmos.color = Color.blue;
        Gizmos.DrawLine(centroid, centroid + normalVector);

        Gizmos.color = Color.yellow;
        Gizmos.DrawLine(centroid, centroid + targetVector);

        Gizmos.color = (isIntersecting) ? Color.red : Color.black;
        Gizmos.DrawSphere(projectionPoint, 0.1f);
    }

    // Update is called once per frame
    void Update() {
        centroid = (vertices[0].position + vertices[1].position + vertices[2].position)/3f;
        normalVector = (vertices[0].forward + vertices[1].forward + vertices[2].forward)/3f;
        normalVector = normalVector.normalized;

        // https://forum.unity.com/threads/projection-of-point-on-plane.855958/
        targetVector = (particleTarget.position - centroid).normalized;
        dotBetweenParticleAndNormal = Vector3.Dot(targetVector, normalVector);
        projectionPoint = ClosestPointOnPlane(centroid, normalVector, particleTarget.position);
        isIntersecting = dotBetweenParticleAndNormal <= 0f 
            && ObstacleHelper.PointInTriangle(
                projectionPoint + ((projectionPoint - particleTarget.position).normalized * radius), 
                vertices[0].position, 
                vertices[1].position, 
                vertices[2].position
            );
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
