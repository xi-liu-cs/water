using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(MeshRenderer))]
[RequireComponent(typeof(MeshFilter))]
[ExecuteInEditMode]
public class PointCloudObstacle : MonoBehaviour
{
    [System.Serializable]
    public class ObstaclePlane {
        public Vector3[] points = new Vector3[3];
        public Plane plane;
        public ObstaclePlane(Plane plane, Vector3[] points) {
            this.plane = plane;
            this.points = points;
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
            Gizmos.DrawSphere(transform.TransformPoint(point), 0.1f);
        }
    }

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
        numPlanes = planesFromMesh.Count;
    }

    private void GenerateGrid() {
        numX = Mathf.CeilToInt((meshRenderer.bounds.max.x - meshRenderer.bounds.min.x) / pointCloudResolution) + 10;
        numY = Mathf.CeilToInt((meshRenderer.bounds.max.y - meshRenderer.bounds.min.y) / pointCloudResolution) + 10;
        numZ = Mathf.CeilToInt((meshRenderer.bounds.max.z - meshRenderer.bounds.min.z) / pointCloudResolution) + 10;
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
            posX = newBoundsMin.x + (x * pointCloudResolution);
            for(int y = 1; y < numY; y++) {
                posY = newBoundsMin.y + (y * pointCloudResolution);
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
                posZ = newBoundsMin.z + (z * pointCloudResolution);
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
            posZ = newBoundsMin.z + (z * pointCloudResolution);
            for(int y = 0; y <= numY; y++) {
                posY = newBoundsMin.y + (y * pointCloudResolution);
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
}