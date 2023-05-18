using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(MeshFilter))]
[RequireComponent(typeof(MeshRenderer))]
[ExecuteInEditMode]
public class PlaneObstacle : MonoBehaviour
{
    [System.Serializable]
    public struct Obs {
        public Vector3 position;
        public Vector4 rotation;
        public Vector3 scale;
    }
    [System.Serializable]
    public struct ObsPlane {
        public Vector3 vertex1, vertex2, vertex3;
        public Vector3 centroid;
        public Vector3 normalVector;
    }

    public enum ObstacleType {
        Static,
        OnTransformChange
    }

    private MeshFilter meshFilter;
    private  MeshRenderer meshRenderer;

    public Obs obs;
    public ObstacleType obstacleType = ObstacleType.Static;
    public List<ObsPlane> obstaclePlanes = new List<ObsPlane>();
    public float dotBetweenParticleAndNormal = 0f;

    public Transform particle = null;
    public bool isIntersecting = false;
    public Vector3 projectionPoint = Vector3.zero;
    public Vector3 worldCentroid = Vector3.zero;

    [SerializeField] private bool printPlaneGizmos = false;

    void OnDrawGizmos() {
        if (obstaclePlanes.Count == 0) return;
        Gizmos.color = Color.black;
        Gizmos.DrawSphere(worldCentroid, 0.025f);
        Gizmos.color = Color.yellow;
        Gizmos.DrawSphere(projectionPoint, 0.01f);
        if (!printPlaneGizmos) return;
        int max = Mathf.Min(obstaclePlanes.Count,1000);
        for(int i = 0; i < max; i++) {
            ObsPlane p = obstaclePlanes[i];
            Gizmos.color = (isIntersecting) ? Color.yellow : Color.red;
            Gizmos.DrawSphere(LocalPointToWorldPoint(obs, p.centroid),0.01f);
            Gizmos.color = Color.blue;
            Gizmos.DrawLine(LocalPointToWorldPoint(obs, p.centroid), LocalPointToWorldPoint(obs, p.centroid + p.normalVector));
        }
    }

    
    private void Awake() {
        meshFilter = GetComponent<MeshFilter>();
        meshRenderer = GetComponent<MeshRenderer>();
        obs = new Obs();
        UpdateObs();
        GeneratePlanes();
    }

    private void Update() {
        // Update `obs` if our transform has changed in some way
        if (transform.hasChanged) {
            UpdateObs();
            if (obstacleType != ObstacleType.Static) GeneratePlanes();
            transform.hasChanged = false;
        }

        // Exit early if we don't have a particle to begin with, or if we don't have any obstacle planes to consider
        if (particle == null || obstaclePlanes.Count == 0) return;

        // First, find closest centroid to current target
        ObsPlane closest = obstaclePlanes[FindClosestPlane(particle.position)];

        // Secondly, convert centroid to world space
        worldCentroid = LocalPointToWorldPoint(obs, closest.centroid);

        // Secondly, Calculate dot product between normal vector and vector b/w centroid and target
        Vector3 targetVector = (particle.position - worldCentroid).normalized;
        dotBetweenParticleAndNormal = Vector3.Dot(targetVector, LocalVectorToWorldVector(obs, closest.normalVector));
        
        // Thirdly, calculate projection onto plane
        projectionPoint = ClosestPointOnPlane(
            worldCentroid, 
            LocalVectorToWorldVector(obs, closest.normalVector), 
            particle.position
        );
        
        float d = DistanceFromPlane(worldCentroid, LocalVectorToWorldVector(obs, closest.normalVector), particle.position);

        // Finally, check if intersecting or not
        isIntersecting = dotBetweenParticleAndNormal <= 0f 
            && ObstacleHelper.PointInTriangle(
                projectionPoint, 
                LocalPointToWorldPoint(obs, closest.vertex1), 
                LocalPointToWorldPoint(obs, closest.vertex2), 
                LocalPointToWorldPoint(obs, closest.vertex3)
        );
    }

    private void UpdateObs() {
        obs.position = transform.position;
        obs.rotation = new Vector4(transform.rotation.x, transform.rotation.y, transform.rotation.z, transform.rotation.w);
        obs.scale = transform.lossyScale;
    }

    public bool GeneratePlanes() {
        Mesh mesh = meshFilter.sharedMesh;
        Vector3 centroid, normalVector;
        ObsPlane newPlane;
        
        obstaclePlanes = new List<ObsPlane>();
        for(int t = 0; t < mesh.triangles.Length; t+=3) {
            centroid = (mesh.vertices[mesh.triangles[t]] + mesh.vertices[mesh.triangles[t+1]] + mesh.vertices[mesh.triangles[t+2]])/3f;
            normalVector = (mesh.normals[mesh.triangles[t]] + mesh.normals[mesh.triangles[t+1]] + mesh.normals[mesh.triangles[t+2]])/3f;
            newPlane = new ObsPlane();
            newPlane.vertex1 = mesh.vertices[mesh.triangles[t]];
            newPlane.vertex2 = mesh.vertices[mesh.triangles[t+1]];
            newPlane.vertex3 = mesh.vertices[mesh.triangles[t+2]];
            newPlane.centroid = centroid;
            newPlane.normalVector = normalVector;
            obstaclePlanes.Add(newPlane);
        }

        return true;
    }

    private int FindClosestPlane(Vector3 t) {
        int closest = 0;
        float distance, closestDistance = Mathf.Infinity;
        for(int i = 0; i < obstaclePlanes.Count; i++) {
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
        return point + DistanceFromPlane(planeOffset, planeNormal.normalized, point) * planeNormal.normalized;
    }
    public float DistanceFromPlane(Vector3 planeOffset, Vector3 planeNormal, Vector3 point) {
        return Vector3.Dot(planeOffset - point, planeNormal);
    }

    // https://forum.unity.com/threads/whats-the-math-behind-transform-transformpoint.107401/
    public static Vector3 LocalPointToWorldPoint(Vector3 pos, Vector4 rot, Vector3 scale, Vector3 localPoint) {
        Vector3 s = new Vector3(localPoint.x * scale.x, localPoint.y * scale.y, localPoint.z * scale.z);
        return RotMultVec3(rot, s) + pos;
    }
    public static Vector3 LocalPointToWorldPoint(Obs obstacle, Vector3 localPoint) {
        Vector3 s = new Vector3(localPoint.x * obstacle.scale.x, localPoint.y * obstacle.scale.y, localPoint.z * obstacle.scale.z);
        return RotMultVec3(obstacle.rotation, s) + obstacle.position;
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
    public static Vector3 LocalVectorToWorldVector(Obs obstacle, Vector3 v) {
        Vector3 start = LocalPointToWorldPoint(obstacle, Vector3.zero);
        Vector3 end = LocalPointToWorldPoint(obstacle, v);
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
