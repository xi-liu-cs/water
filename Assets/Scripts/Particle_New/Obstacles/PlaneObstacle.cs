using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(MeshFilter))]
[RequireComponent(typeof(MeshRenderer))]
[ExecuteInEditMode]
public class PlaneObstacle : MonoBehaviour
{

    public struct ObsPlane {
        public Vector3[] vertices;
        public Vector3 centroid;
        public Vector3 normalVector;
    }

    public enum ObstacleType {
        Static,
        OnTransformChange
    }

    [ReadOnly, SerializeField] private MeshFilter meshFilter;
    [ReadOnly, SerializeField] private  MeshRenderer meshRenderer;
    private Vector4 rotationV4;
    public ObstacleType obstacleType = ObstacleType.Static;
    public ObsPlane[] obstaclePlanes;

    public Transform particle = null;
    private bool isIntersecting = false;

    void OnDrawGizmos() {
        if (obstaclePlanes.Length == 0) return;
        foreach(ObsPlane p in obstaclePlanes) {
            Gizmos.color = Color.green;

            Gizmos.DrawSphere(LocalPointToWorldPoint(transform.position, rotationV4, transform.localScale, p.vertices[0]),0.01f);
            Gizmos.DrawSphere(LocalPointToWorldPoint(transform.position, rotationV4, transform.localScale, p.vertices[1]),0.01f);
            Gizmos.DrawSphere(LocalPointToWorldPoint(transform.position, rotationV4, transform.localScale, p.vertices[2]),0.01f);
            Gizmos.color = (isIntersecting) ? Color.yellow : Color.red;
            Gizmos.DrawSphere(LocalPointToWorldPoint(transform.position, rotationV4, transform.localScale, p.centroid),0.01f);
            Gizmos.color = Color.blue;
            Gizmos.DrawLine(LocalPointToWorldPoint(transform.position, rotationV4, transform.localScale, p.centroid), LocalPointToWorldPoint(transform.position, rotationV4, transform.localScale, p.centroid + p.normalVector));
        }
    }

    
    private void Awake() {
        meshFilter = GetComponent<MeshFilter>();
        meshRenderer = GetComponent<MeshRenderer>();
        GeneratePlanes();
    }

    private void Update() {
        rotationV4 =  new Vector4(transform.rotation.x, transform.rotation.y, transform.rotation.z, transform.rotation.w);
        if (obstacleType != ObstacleType.Static && transform.hasChanged) {
            GeneratePlanes();
            transform.hasChanged = false;
        }

        if (particle == null) return;

        // First, find closest centroid to current target
        ObsPlane closest = obstaclePlanes[FindClosestPlane(particle.position)];

        // Secondly, convert centroid to world space
        Vector3 worldCentroid = LocalPointToWorldPoint(transform.position, rotationV4, transform.localScale, closest.centroid);

        // Secondly, Calculate dot product between normal vector and vector b/w centroid and target
        Vector3 targetVector = (particle.position - worldCentroid).normalized;
        float dotBetweenParticleAndNormal = Vector3.Dot(targetVector, LocalVectorToWorldVector(transform.position, rotationV4, transform.localScale, closest.normalVector));
        
        // Thirdly, calculate projection onto plane
        Vector3 projectionPoint = ClosestPointOnPlane(worldCentroid, LocalVectorToWorldVector(transform.position, rotationV4, transform.localScale, closest.normalVector), particle.position);
        
        // Finally, check if intersecting or not
        isIntersecting = dotBetweenParticleAndNormal <= 0f 
            && ObstacleHelper.PointInTriangle(
                projectionPoint, 
                transform.TransformPoint(closest.vertices[0]), 
                transform.TransformPoint(closest.vertices[1]), 
                transform.TransformPoint(closest.vertices[2])
        );
    }

    private void GeneratePlanes() {
        Mesh mesh = meshFilter.sharedMesh;
        Vector3 centroid, normalVector;
        ObsPlane newPlane;
        
        obstaclePlanes = new ObsPlane[mesh.triangles.Length / 3];
        for(int t = 0; t < mesh.triangles.Length; t+=3) {
            centroid = (mesh.vertices[mesh.triangles[t]] + mesh.vertices[mesh.triangles[t+1]] + mesh.vertices[mesh.triangles[t+2]])/3f;
            normalVector = (mesh.normals[mesh.triangles[t]] + mesh.normals[mesh.triangles[t+1]] + mesh.normals[mesh.triangles[t+2]])/3f;
            newPlane = new ObsPlane();
            newPlane.vertices = new Vector3[3]{
                mesh.vertices[mesh.triangles[t]],
                mesh.vertices[mesh.triangles[t+1]],
                mesh.vertices[mesh.triangles[t+2]]
            };
            newPlane.centroid = centroid;
            newPlane.normalVector = normalVector;
            obstaclePlanes[t/3] = newPlane;
        }
    }

    private int FindClosestPlane(Vector3 t) {
        int closest = 0;
        float distance, closestDistance = Mathf.Infinity;
        for(int i = 0; i < obstaclePlanes.Length; i++) {
            distance = Vector3.Distance(t, transform.TransformPoint(obstaclePlanes[i].centroid));
            if (distance < closestDistance) {
                closestDistance = distance;
                closest = i;
            }
        }
        return closest;
    }

    // https://forum.unity.com/threads/projection-of-point-on-plane.855958/
    public Vector3 ClosestPointOnPlane(Vector3 planeOffset, Vector3 planeNormal, Vector3 point) {
        return point + DistanceFromPlane(planeOffset, planeNormal, point) * planeNormal;
    }
    public float DistanceFromPlane(Vector3 planeOffset, Vector3 planeNormal, Vector3 point) {
        return Vector3.Dot(planeOffset - point, planeNormal);
    }

    // https://forum.unity.com/threads/whats-the-math-behind-transform-transformpoint.107401/
    public static Vector3 LocalPointToWorldPoint(Vector3 pos, Vector4 rot, Vector3 scale, Vector3 localPoint) {
        Vector3 s = new Vector3(localPoint.x * scale.x, localPoint.y * scale.y, localPoint.z * scale.z);
        return RotMultVec3(rot, s) + pos;
    }
    // https://answers.unity.com/questions/372371/multiply-quaternion-by-vector3-how-is-done.html
    public static Vector3 RotMultVec3(Vector4 quat, Vector3 vec){
        float num = quat.x * 2f;
        float num2 = quat.y * 2f;
        float num3 = quat.z * 2f;
        float num4 = quat.x * num;
        float num5 = quat.y * num2;
        float num6 = quat.z * num3;
        float num7 = quat.x * num2;
        float num8 = quat.x * num3;
        float num9 = quat.y * num3;
        float num10 = quat.w * num;
        float num11 = quat.w * num2;
        float num12 = quat.w * num3;
        Vector3 result;
        result.x = (1f - (num5 + num6)) * vec.x + (num7 - num12) * vec.y + (num8 + num11) * vec.z;
        result.y = (num7 + num12) * vec.x + (1f - (num4 + num6)) * vec.y + (num9 - num10) * vec.z;
        result.z = (num8 - num11) * vec.x + (num9 + num10) * vec.y + (1f - (num4 + num5)) * vec.z;
        return result;
    }
    public static Vector3 LocalVectorToWorldVector(Vector3 pos, Vector4 rot, Vector3 scale, Vector3 v) {
        Vector3 start = LocalPointToWorldPoint(pos, rot, scale, Vector3.zero);
        Vector3 end = LocalPointToWorldPoint(pos, rot, scale, v);
        return end - start;
    }
    // https://github.com/HardlyDifficult/Tutorials/blob/master/Quaternions.md#361-quaternioninverse
    public static Quaternion InverseQuaternion(Quaternion rotation) {
        // Split the Quaternion component
        Vector3 vector = new Vector3(rotation.x, rotation.y, rotation.z);
        float scalar = rotation.w;
        // Calculate inverse
        vector = -vector;
        // Store results
        Quaternion inverseRotation = new Quaternion(vector.x, vector.y, vector.z, scalar);
        return inverseRotation;
    }
}
