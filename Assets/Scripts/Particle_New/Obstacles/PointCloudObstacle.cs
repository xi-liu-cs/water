using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using System.Linq;

[RequireComponent(typeof(MeshRenderer))]
[RequireComponent(typeof(MeshFilter))]
[ExecuteInEditMode]
public class PointCloudObstacle : MonoBehaviour
{
    [System.Serializable]
    public class ObstaclePlane {
        public Vector3[] vertices = new Vector3[3];
        public Plane plane;
        public Vector3 a1, a2, a1Closest;
        public int a1Res, a2Res;
        public ObstaclePlane(Vector3[] vertices, Plane plane, Vector3 a1, Vector3 a2, Vector3 a1Closest) {
            this.vertices = vertices;
            this.plane = plane;
            this.a1 = a1;
            this.a2 = a2;
            this.a1Closest = a1Closest;
        }
        public float A1_NumCells(float res) => Mathf.Ceil(a1.magnitude / (res));
        public Vector3 A1_Direction => a1.normalized;
        public float A2_NumCells(float res) => Mathf.Ceil(a2.magnitude / (res));
        public Vector3 A2_Direction => a2.normalized;
    }

    [System.Serializable]
    public struct Point {
        public int vertexIndex;
        public Vector3 worldPosition, localPosition;
        public List<int> neighbors;
        public Point(int vertexIndex, Vector3 worldPosition, Vector3 localPosition) {
            this.vertexIndex = vertexIndex;
            this.worldPosition = worldPosition;
            this.localPosition = localPosition;
            this.neighbors = new List<int>();
        }
    }

    public enum GenerateFreq {
        Always,
        OnChange,
        OnceOnly
    }
    public enum FilterSettings {
        FilterClosePoints,
        NoFiltering
    }

    public GenerateFreq generationFrequency = GenerateFreq.Always;
    public FilterSettings filterSettings = FilterSettings.FilterClosePoints;
    public MeshFilter meshFilter;
    public MeshRenderer meshRenderer;
    public List<ObstaclePlane> planesFromMesh = new List<ObstaclePlane>();
    [ReadOnly, SerializeField]
    private int numPlanes;
    [ReadOnly, SerializeField]
    Vector3 newBounds, newBoundsMin, newBoundsMax;
    public List<Vector3> points = new List<Vector3>();
    public List<Point> pPoints = new List<Point>();
    public float3[] pointsF3 = new float3[0];

    public float pointCloudResolution = 1f;
    public int numX, numY, numZ;

    
    void OnDrawGizmos() {
        MeshRenderer r = GetComponent<MeshRenderer>();
        if (r == null) return;
        var bounds = r.bounds;
        Gizmos.matrix = Matrix4x4.identity;
        Gizmos.color = Color.blue;
        Gizmos.DrawWireCube(bounds.center, bounds.extents * 2);
        Gizmos.color = Color.yellow;
        Gizmos.DrawSphere(transform.position - newBounds/2f, 0.5f);
        Gizmos.DrawSphere(transform.position + newBounds/2f, 0.5f);
        if (points.Count == 0) return;
        Gizmos.color = new Vector4(1f,0f,0f,1f);
        foreach(Vector3 point in points) {
            Gizmos.DrawSphere(transform.TransformPoint(point), pointCloudResolution);
        }
    }
    

    private void Awake() {
        meshFilter = GetComponent<MeshFilter>();
        meshRenderer = GetComponent<MeshRenderer>();
    }

    private void Update() {
        switch(generationFrequency) {
            case GenerateFreq.Always:
                //GeneratePlanes();
                //GenerateGrid();
                GeneratePoints();
                break;
            case GenerateFreq.OnChange:
                if (!transform.hasChanged) return;
                //GeneratePlanes();
                //GenerateGrid();
                GeneratePoints();
                transform.hasChanged = false;
                break;
            case GenerateFreq.OnceOnly:
                break;
        }
    }

    private void GeneratePoints() {
        // Initializations
        points = new List<Vector3>();
        List<float3> tempPointsF3 = new List<float3>();

        Mesh mesh = meshFilter.sharedMesh;
        int[] triangles = mesh.triangles;
        Vector3[] vertices = mesh.vertices;

        // world positions for each triangle vertex, as well as centroid
        Vector3 wv1, wv2, wv3, centroid;
        // vector b/w wv1 and wv2, and wv3-closestPointOn_wv1-wv2
        Vector3 wv1wv2, wv3_closest, wv3_wv1wv2;
        // How many points can we fit on each axis?
        int wv1wv2_num, wv3_wv1wv2_num;
        // When we're adding points, we need reference points!
        Vector3 tempPos;
        // Plane data
        Plane plane;

        for(int t = 0; t < triangles.Length; t += 3) {
            // Calculate world positions for each triangle vertex
            wv1 = transform.TransformPoint(mesh.vertices[triangles[t]]);
            wv2 = transform.TransformPoint(mesh.vertices[triangles[t+1]]);
            wv3 = transform.TransformPoint(mesh.vertices[triangles[t+2]]);
            // Calculate centroid of three points
            centroid = (wv1 + wv2 + wv3)/3f;
            // Calculate world-space vector b/w wv1 and wv2
            wv1wv2 = wv2 - wv1;
            // Calculate closest point to wv3 from `wv1_wv2`
            wv3_closest = ObstacleHelper.FindNearestPointOnLine(wv1,wv2,wv3);
            // Calculate vector from wv3 to wv3_closest
            wv3_wv1wv2 = wv3 - wv3_closest;

            // Calculate how many points we are allowed on each axis
            wv1wv2_num = Mathf.CeilToInt(wv1wv2.magnitude / pointCloudResolution);
            wv3_wv1wv2_num = Mathf.CeilToInt(wv3_wv1wv2.magnitude / pointCloudResolution);

            // Get the plane that encompasses wv1,wv2,wv3
            plane = new Plane(wv1, wv2, wv3);

            // Add centroid
            points.Add(transform.InverseTransformPoint(centroid));

            // Add points by iterating across rectangle
            for(int x = 0; x < wv1wv2_num; x++) {
                float wv3_wv1wv2_oddOffset = (x % 2 == 1) ? pointCloudResolution : 0f;
                for(int y = 0; y < wv3_wv1wv2_num; y++) {
                    tempPos = wv1 + (wv1wv2.normalized * (x*pointCloudResolution + pointCloudResolution)) + (wv3_wv1wv2.normalized * (y*pointCloudResolution + pointCloudResolution));
                    tempPos = plane.ClosestPointOnPlane(tempPos);
                    if(ObstacleHelper.PointInTriangle(tempPos, wv1,wv2,wv3)) {
                        // Convert to local space
                        points.Add(transform.InverseTransformPoint(tempPos));
                    }
                 }
            }
        }
    }
    
    private void GeneratePlanes() {
        planesFromMesh = new List<ObstaclePlane>();
        Mesh mesh = meshFilter.sharedMesh;
        int[] triangles = mesh.triangles;
        Vector3[] worldPoints;
        Plane plane;
        Vector3 p1, p2, p3, a1, a2, a1Closest;
        for(int t = 0; t < triangles.Length; t += 3) {
            // Get world space coordinates of each triangle vertex at this frame
            worldPoints = new Vector3[3] {
                transform.TransformPoint(mesh.vertices[triangles[t]]),
                transform.TransformPoint(mesh.vertices[triangles[t+1]]),
                transform.TransformPoint(mesh.vertices[triangles[t+2]])
            };
            // Calculate vector b/w p1 and p2
            a1 = worldPoints[1] - worldPoints[0];
            // Find the closest point between p3 and the vector p1p2
            a1Closest = ObstacleHelper.FindNearestPointOnLine(worldPoints[0],worldPoints[1],worldPoints[2]);
            a2 = worldPoints[2] - a1Closest;
            // Get the plane from each world point
            plane = new Plane(worldPoints[0], worldPoints[1], worldPoints[2]);
            // Save these results
            planesFromMesh.Add(new ObstaclePlane(worldPoints, plane, a1, a2, a1Closest));
        }
        numPlanes = planesFromMesh.Count;
    }

    private void GenerateGrid() {
        Vector3 testPos, closestPointOnPlane;
        points = new List<Vector3>();
        List<float3> tempPointsF3 = new List<float3>();
        
        for(int i = 0; i < planesFromMesh.Count; i++) {
            
            // Place at least one point on the center of the mesh
            closestPointOnPlane = planesFromMesh[i].vertices[0] 
                + planesFromMesh[i].a1.normalized * planesFromMesh[i].a1.magnitude * 0.5f
                + planesFromMesh[i].a2.normalized * planesFromMesh[i].a2.magnitude * 0.5f;
            points.Add(transform.InverseTransformPoint(closestPointOnPlane));
            tempPointsF3.Add(new(closestPointOnPlane.x, closestPointOnPlane.y, closestPointOnPlane.z));

            for(int x = 0; x < planesFromMesh[i].A1_NumCells(pointCloudResolution); x++) {
                for(int y = 0; y < planesFromMesh[i].A2_NumCells(pointCloudResolution); y++) {
                    testPos = planesFromMesh[i].vertices[0] 
                        + (planesFromMesh[i].a1.normalized * (x*pointCloudResolution + pointCloudResolution)) 
                        + (planesFromMesh[i].a2.normalized * (y*pointCloudResolution + pointCloudResolution));
                    closestPointOnPlane = planesFromMesh[i].plane.ClosestPointOnPlane(testPos);
                    if(ObstacleHelper.PointInTriangle(closestPointOnPlane,planesFromMesh[i].vertices)) {
                        // Convert to local space
                        points.Add(transform.InverseTransformPoint(closestPointOnPlane));
                        tempPointsF3.Add(new(closestPointOnPlane.x, closestPointOnPlane.y, closestPointOnPlane.z));
                    }
                 }
            }
        }

        pointsF3 = tempPointsF3.ToArray();
        

        foreach(ObstaclePlane obstaclePlane in planesFromMesh) {
            
            for(int x = 0; x < obstaclePlane.A1_NumCells(pointCloudResolution*2f); x++) {
                for(int y = 0; y < obstaclePlane.A2_NumCells(pointCloudResolution*2f); y++) {
                    testPos = obstaclePlane.vertices[0] 
                        + (obstaclePlane.A1_Direction * (x*pointCloudResolution*2f + pointCloudResolution)) 
                        + (obstaclePlane.A2_Direction * (y*pointCloudResolution*2f + pointCloudResolution));
                    closestPointOnPlane = obstaclePlane.plane.ClosestPointOnPlane(testPos);
                    
                    if(ObstacleHelper.PointInTriangle(closestPointOnPlane,obstaclePlane.vertices)) {
                        points.Add(transform.InverseTransformPoint(closestPointOnPlane));
                    }
                    
                    points.Add(transform.InverseTransformPoint(testPos));
                }
            }
            
            //points.Add(transform.InverseTransformPoint(obstaclePlane.vertices[0]));
           // points.Add(transform.InverseTransformPoint(obstaclePlane.vertices[1]));
            //points.Add(transform.InverseTransformPoint(obstaclePlane.vertices[0]+obstaclePlane.a1.normalized*0.1f));
            //points.Add(transform.InverseTransformPoint(obstaclePlane.a1Closest));
            //points.Add(transform.InverseTransformPoint(obstaclePlane.vertices[0]+obstaclePlane.a2.normalized*0.1f));
        }
    }
}

public class ObstacleHelper {
    public static bool SameSide(Vector3 p1, Vector3 p2, Vector3 a, Vector3 b) {
        Vector3 cp1 = Vector3.Cross(b-a, p1-a);
        Vector3 cp2 = Vector3.Cross(b-a, p2-a);
        return Vector3.Dot(cp1, cp2) >= 0;
    }

    public static bool PointInTriangle(Vector3 p, Vector3 a, Vector3 b, Vector3 c) {
        return SameSide(p,a, b,c) && SameSide(p,b, a,c) && SameSide(p,c, a,b);
    }
    public static bool PointInTriangle(Vector3 p, Vector3[] triangle) {
        if (triangle.Length != 3) return false;
        return SameSide(p,triangle[0], triangle[1],triangle[2]) && SameSide(p,triangle[1], triangle[0],triangle[2]) && SameSide(p,triangle[2], triangle[0],triangle[1]);
    }

    public static Vector3 FindNearestPointOnLine(Vector3 origin, Vector3 end, Vector3 point) {
        //Get heading
        Vector3 heading = (end - origin);
        float magnitudeMax = heading.magnitude;
        heading.Normalize();

        //Do projection from the point but clamp it
        Vector3 lhs = point - origin;
        float dotP = Vector3.Dot(lhs, heading);
        dotP = Mathf.Clamp(dotP, 0f, magnitudeMax);
        return origin + heading * dotP;
    }
}