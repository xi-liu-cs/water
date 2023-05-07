using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

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
    public float3[] pointsF3 = new float3[0];

    public float pointCloudResolution = 1f;
    public int numX, numY, numZ;

    /*
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
            Gizmos.DrawSphere(transform.TransformPoint(point), 0.1f);
        }
    }
    */

    private void Awake() {
        meshFilter = GetComponent<MeshFilter>();
        meshRenderer = GetComponent<MeshRenderer>();
    }

    private void Update() {
        if (meshFilter == null) meshFilter = GetComponent<MeshFilter>();
        if (meshRenderer == null) meshRenderer = GetComponent<MeshRenderer>();
        switch(generationFrequency) {
            case GenerateFreq.Always:
                GeneratePlanes();
                GenerateGrid();
                break;
            case GenerateFreq.OnChange:
                if (!transform.hasChanged) return;
                GeneratePlanes();
                GenerateGrid();
                transform.hasChanged = false;
                break;
            case GenerateFreq.OnceOnly:
                break;
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
            // Calculate magnitude and direction of p1p2
            //a1_mag = Mathf.Ceil(p1p2.magnitude / (pointCloudResolution*2f)) * (pointCloudResolution*2f);
            //a1_dir = p1p2.normalized;
            // Find the closest point between p3 and the vector p1p2
            a1Closest = ObstacleHelper.FindNearestPointOnLine(worldPoints[0],worldPoints[1],worldPoints[2]);
            a2 = worldPoints[2] - a1Closest;
            // Calcualte magnitude and direction of p3_p1p2
            //a2_mag = Mathf.Ceil(p3_p1p2.magnitude / (pointCloudResolution*2f)) * (pointCloudResolution*2f);
            //a2_dir = p3_p1p2.normalized;
            // Get the plane from each world point
            plane = new Plane(worldPoints[0], worldPoints[1], worldPoints[2]);
            // Save these results
            planesFromMesh.Add(new ObstaclePlane(worldPoints, plane, a1, a2, a1Closest));
        }
        /*
        for(int t = 0; t < triangles.Length; t+=3) {
            // triangles[t] = first index, trianges[t+1] = second index, triangles[t+2] = 3rd index
            worldPoints = new Vector3[3] {
                transform.TransformPoint(mesh.vertices[triangles[t]]),
                transform.TransformPoint(mesh.vertices[triangles[t+1]]),
                transform.TransformPoint(mesh.vertices[triangles[t+2]])
            };
            Plane plane = new Plane(worldPoints[0], worldPoints[1], worldPoints[2]);
            planesFromMesh.Add(new ObstaclePlane(plane, worldPoints));
        }
        */
        numPlanes = planesFromMesh.Count;
    }

    private void GenerateGrid() {
        Vector3 testPos, closestPointOnPlane;
        points = new List<Vector3>();
        
        for(int i = 0; i < planesFromMesh.Count; i++) {
            /*
            points.Add(transform.InverseTransformPoint(planesFromMesh[i].vertices[0]));
            //points.Add(transform.InverseTransformPoint(planesFromMesh[i].vertices[1]));
            points.Add(transform.InverseTransformPoint(planesFromMesh[i].vertices[2]));
            //points.Add(transform.InverseTransformPoint(planesFromMesh[i].a1Closest));
            points.Add(transform.InverseTransformPoint(planesFromMesh[i].vertices[0]+planesFromMesh[i].a1.normalized*0.1f));
            points.Add(transform.InverseTransformPoint(planesFromMesh[i].vertices[0]+planesFromMesh[i].a2.normalized*0.1f));
            */
            // Place at least one point on the center of the mesh
            points.Add(transform.InverseTransformPoint(
                planesFromMesh[i].vertices[0] 
                + planesFromMesh[i].a1.normalized * planesFromMesh[i].a1.magnitude * 0.05f
                + planesFromMesh[i].a2.normalized * planesFromMesh[i].a2.magnitude * 0.05f 
            ));
            for(int x = 0; x < planesFromMesh[i].A1_NumCells(pointCloudResolution); x++) {
                for(int y = 0; y < planesFromMesh[i].A2_NumCells(pointCloudResolution); y++) {
                    testPos = planesFromMesh[i].vertices[0] 
                        + (planesFromMesh[i].a1.normalized * (x*pointCloudResolution + pointCloudResolution)) 
                        + (planesFromMesh[i].a2.normalized * (y*pointCloudResolution + pointCloudResolution));
                    closestPointOnPlane = planesFromMesh[i].plane.ClosestPointOnPlane(testPos);
                    if(ObstacleHelper.PointInTriangle(closestPointOnPlane,planesFromMesh[i].vertices)) {
                        points.Add(transform.InverseTransformPoint(closestPointOnPlane));
                    }
                 }
            }
        }
        
        pointsF3 = new float3[points.Count];
        for(int i = 0; i < points.Count; i++) {
            Vector3 worldPoint = transform.TransformPoint(points[i]);
            pointsF3[i] = new(worldPoint.x, worldPoint.y, worldPoint.z);
        }
        

        //foreach(ObstaclePlane obstaclePlane in planesFromMesh) {
            /*
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
            */
            //points.Add(transform.InverseTransformPoint(obstaclePlane.vertices[0]));
           // points.Add(transform.InverseTransformPoint(obstaclePlane.vertices[1]));
            //points.Add(transform.InverseTransformPoint(obstaclePlane.vertices[0]+obstaclePlane.a1.normalized*0.1f));
            //points.Add(transform.InverseTransformPoint(obstaclePlane.a1Closest));
            //points.Add(transform.InverseTransformPoint(obstaclePlane.vertices[0]+obstaclePlane.a2.normalized*0.1f));
        //}

        /*
        numX = Mathf.CeilToInt((meshRenderer.bounds.max.x - meshRenderer.bounds.min.x) / (pointCloudResolution*2f)) + 10;
        numY = Mathf.CeilToInt((meshRenderer.bounds.max.y - meshRenderer.bounds.min.y) / (pointCloudResolution*2f)) + 10;
        numZ = Mathf.CeilToInt((meshRenderer.bounds.max.z - meshRenderer.bounds.min.z) / (pointCloudResolution*2f)) + 10;
        Vector3Int numGridCellsPerAxis = new Vector3Int(numX, numY, numZ);

        newBounds = new Vector3(
            (float)numX * pointCloudResolution,
            (float)numY * pointCloudResolution,
            (float)numZ * pointCloudResolution
        );

        Vector3 halfBounds = newBounds / 2f;
        newBoundsMin = transform.position - halfBounds;
        newBoundsMax = transform.position + halfBounds;
        Vector3 boundPos, closestPointOnPlane;
        float posX, posY, posZ;
        int projectedIndex;

        points = new List<Vector3>();
        bool[] pointsOccupied = new bool[numX * numY * numZ];

        for(int x = 0; x <= numX; x++) {
            posX = newBoundsMin.x + (x * pointCloudResolution*2f) + pointCloudResolution;
            for(int y = 1; y < numY; y++) {
                posY = newBoundsMin.y + (y * pointCloudResolution*2f) + pointCloudResolution;
                boundPos = new Vector3(posX, posY, newBoundsMin.z);
                foreach(ObstaclePlane obsPlane in planesFromMesh) {
                    closestPointOnPlane = obsPlane.plane.ClosestPointOnPlane(boundPos);
                    if (ObstacleHelper.PointInTriangle(closestPointOnPlane,obsPlane.points)) {
                        if (filterSettings == FilterSettings.FilterClosePoints) {
                            projectedIndex = Grid3DHelpers.GetProjectedGridIndexFromGivenPosition(numGridCellsPerAxis,  transform.position, pointCloudResolution, closestPointOnPlane);
                            if(!pointsOccupied[projectedIndex]) {
                                points.Add(transform.InverseTransformPoint(closestPointOnPlane));
                                pointsOccupied[projectedIndex] = true;
                            }
                        } else {
                            points.Add(transform.InverseTransformPoint(closestPointOnPlane));
                        }
                    }
                }
                boundPos = new Vector3(posX, posY, newBoundsMax.z);
                foreach(ObstaclePlane obsPlane in planesFromMesh) {
                    closestPointOnPlane = obsPlane.plane.ClosestPointOnPlane(boundPos);
                    if (ObstacleHelper.PointInTriangle(closestPointOnPlane,obsPlane.points)) {
                        if (filterSettings == FilterSettings.FilterClosePoints) {
                            projectedIndex = Grid3DHelpers.GetProjectedGridIndexFromGivenPosition(numGridCellsPerAxis,  transform.position, pointCloudResolution, closestPointOnPlane);
                            if(!pointsOccupied[projectedIndex]) {
                                points.Add(transform.InverseTransformPoint(closestPointOnPlane));
                                pointsOccupied[projectedIndex] = true;
                            }
                        } else {
                            points.Add(transform.InverseTransformPoint(closestPointOnPlane));
                        }
                    }
                }
            }
            for(int z = 0; z <= numZ; z++) {
                posZ = newBoundsMin.z + (z * pointCloudResolution*2f) + pointCloudResolution;
                boundPos = new Vector3(posX, newBoundsMin.y, posZ);
                foreach(ObstaclePlane obsPlane in planesFromMesh) {
                    closestPointOnPlane = obsPlane.plane.ClosestPointOnPlane(boundPos);
                    if (ObstacleHelper.PointInTriangle(closestPointOnPlane,obsPlane.points)) {
                        if (filterSettings == FilterSettings.FilterClosePoints) {
                            projectedIndex = Grid3DHelpers.GetProjectedGridIndexFromGivenPosition(numGridCellsPerAxis,  transform.position, pointCloudResolution, closestPointOnPlane);
                            if(!pointsOccupied[projectedIndex]) {
                                points.Add(transform.InverseTransformPoint(closestPointOnPlane));
                                pointsOccupied[projectedIndex] = true;
                            }
                        } else {
                            points.Add(transform.InverseTransformPoint(closestPointOnPlane));
                        }
                    }
                }
                boundPos = new Vector3(posX, newBoundsMax.y, posZ);
                foreach(ObstaclePlane obsPlane in planesFromMesh) {
                    closestPointOnPlane = obsPlane.plane.ClosestPointOnPlane(boundPos);
                    if (ObstacleHelper.PointInTriangle(closestPointOnPlane,obsPlane.points)) {
                        if (filterSettings == FilterSettings.FilterClosePoints) {
                            projectedIndex = Grid3DHelpers.GetProjectedGridIndexFromGivenPosition(numGridCellsPerAxis,  transform.position, pointCloudResolution, closestPointOnPlane);
                            if(!pointsOccupied[projectedIndex]) {
                                points.Add(transform.InverseTransformPoint(closestPointOnPlane));
                                pointsOccupied[projectedIndex] = true;
                            }
                        } else {
                            points.Add(transform.InverseTransformPoint(closestPointOnPlane));
                        }
                    }
                }
            }
        }
        for(int z = 0; z <= numZ; z++) {
            posZ = newBoundsMin.z + (z * pointCloudResolution*2f) + pointCloudResolution;
            for(int y = 0; y <= numY; y++) {
                posY = newBoundsMin.y + (y * pointCloudResolution*2f) + pointCloudResolution;
                boundPos = new Vector3(newBoundsMin.x, posY, posZ);
                foreach(ObstaclePlane obsPlane in planesFromMesh) {
                    closestPointOnPlane = obsPlane.plane.ClosestPointOnPlane(boundPos);
                    if (ObstacleHelper.PointInTriangle(closestPointOnPlane, obsPlane.points)) {
                        if (filterSettings == FilterSettings.FilterClosePoints) {
                            projectedIndex = Grid3DHelpers.GetProjectedGridIndexFromGivenPosition(numGridCellsPerAxis,  transform.position, pointCloudResolution, closestPointOnPlane);
                            if(!pointsOccupied[projectedIndex]) {
                                points.Add(transform.InverseTransformPoint(closestPointOnPlane));
                                pointsOccupied[projectedIndex] = true;
                            }
                        } else {
                            points.Add(transform.InverseTransformPoint(closestPointOnPlane));
                        }
                    }
                }
                boundPos = new Vector3(newBoundsMax.x, posY, posZ);
                foreach(ObstaclePlane obsPlane in planesFromMesh) {
                    closestPointOnPlane = obsPlane.plane.ClosestPointOnPlane(boundPos);
                    if (ObstacleHelper.PointInTriangle(closestPointOnPlane, obsPlane.points)) {
                        if (filterSettings == FilterSettings.FilterClosePoints) {
                            projectedIndex = Grid3DHelpers.GetProjectedGridIndexFromGivenPosition(numGridCellsPerAxis,  transform.position, pointCloudResolution, closestPointOnPlane);
                            if(!pointsOccupied[projectedIndex]) {
                                points.Add(transform.InverseTransformPoint(closestPointOnPlane));
                                pointsOccupied[projectedIndex] = true;
                            }
                        } else {
                            points.Add(transform.InverseTransformPoint(closestPointOnPlane));
                        }
                    }
                }
            }
        }
        */
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