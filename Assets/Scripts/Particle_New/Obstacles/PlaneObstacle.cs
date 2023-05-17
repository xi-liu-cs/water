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
    public ObstacleType obstacleType = ObstacleType.Static;
    public ObsPlane[] obstaclePlanes;

    public Transform particle = null;
    private bool isIntersecting = false;

    void OnDrawGizmos() {
        if (obstaclePlanes.Length == 0) return;
        foreach(ObsPlane p in obstaclePlanes) {
            Gizmos.color = Color.green;
            Gizmos.DrawSphere(transform.TransformPoint(p.vertices[0]),0.01f);
            Gizmos.DrawSphere(transform.TransformPoint(p.vertices[1]),0.01f);
            Gizmos.DrawSphere(transform.TransformPoint(p.vertices[2]),0.01f);
            Gizmos.color = (isIntersecting) ? Color.yellow : Color.red;
            Gizmos.DrawSphere(transform.TransformPoint(p.centroid),0.01f);
            Gizmos.color = Color.blue;
            Gizmos.DrawLine(transform.TransformPoint(p.centroid), transform.TransformPoint(p.centroid + p.normalVector));
        }
    }

    
    private void Awake() {
        meshFilter = GetComponent<MeshFilter>();
        meshRenderer = GetComponent<MeshRenderer>();
        GeneratePlanes();
    }

    private void Update() {
        if (obstacleType != ObstacleType.Static && transform.hasChanged) {
            GeneratePlanes();
            transform.hasChanged = false;
        }

        if (particle == null) return;

        // First, find closest centroid to current target
        ObsPlane closest = obstaclePlanes[FindClosestPlane(particle.position)];

        // Secondly, Calculate dot product between normal vector and vector b/w centroid and target
        Vector3 targetVector = (particle.position - transform.TransformPoint(closest.centroid)).normalized;
        float dotBetweenParticleAndNormal = Vector3.Dot(targetVector, transform.TransformVector(closest.normalVector));
        
        // Thirdly, calculate projection onto plane
        Vector3 projectionPoint = ClosestPointOnPlane(transform.TransformPoint(closest.centroid), transform.TransformVector(closest.normalVector), particle.position);
        
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

    public Vector3 ClosestPointOnPlane(Vector3 planeOffset, Vector3 planeNormal, Vector3 point) {
        return point + DistanceFromPlane(planeOffset, planeNormal, point) * planeNormal;
    }

    public float DistanceFromPlane(Vector3 planeOffset, Vector3 planeNormal, Vector3 point) {
        return Vector3.Dot(planeOffset - point, planeNormal);
    }
}
