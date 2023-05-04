using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TestGrid : MonoBehaviour
{

    [Header("GRID PROPERTIES")]
    public Vector3 bounds = new Vector3(10f,10f,10f);
    public float gridCellSize = 0.5f;
    public Vector3Int numBufferCellsPerAxis = new Vector3Int(10,10,10);
    
    private Vector3 gridCellSize3D;
    private float gridCellDiagonalHalf;
    private Vector3Int numCellsPerAxis;
    private int numGridCells;
    private Vector3[] obstacleGrid;

    public List<Transform> particles = new List<Transform>();
    public Vector3[] particlePositions = new Vector3[0];
    
    public Color cellColor = Color.red;

    [System.Serializable]
    public class Obstacle {
        public MeshFilter obstacle;
        public bool active;

        public Vector3Int currentCenterGridIndices;
        public Quaternion currentRotation;
        public Vector3 currentScale;
        public Vector3[] vertexGridWorldPositions;
    }
    public Obstacle[] obstacles = new Obstacle[0];
    public MeshFilter testMesh = null;

    public bool verbose = false;

    void OnDrawGizmos() {

        Vector3Int dimensions = new Vector3Int(
            Mathf.CeilToInt(bounds.x/gridCellSize)+numBufferCellsPerAxis.x,
            Mathf.CeilToInt(bounds.y/gridCellSize)+numBufferCellsPerAxis.y,
            Mathf.CeilToInt(bounds.z/gridCellSize)+numBufferCellsPerAxis.z
        );

        Vector3 minBound = Boid3DHelpers.GetGridCellWorldPositionFromGivenPosition(dimensions, transform.position, gridCellSize, new Vector3(-bounds.x/2f, -bounds.y/2f, -bounds.z/2f));
        Vector3 maxBound = Boid3DHelpers.GetGridCellWorldPositionFromGivenPosition(dimensions, transform.position, gridCellSize, new Vector3(bounds.x/2f, bounds.y/2f, bounds.z/2f));
        Vector3 gridBounds = new Vector3(
            Mathf.Abs(maxBound.x - minBound.x),
            Mathf.Abs(maxBound.y - minBound.y),
            Mathf.Abs(maxBound.z - minBound.z)
        );
        Gizmos.color = Color.white;
        Gizmos.DrawWireCube(transform.position, gridBounds);

        Gizmos.color = Color.yellow;
        Gizmos.DrawWireCube(transform.position, new Vector3(
            gridCellSize * dimensions.x,
            gridCellSize * dimensions.y,
            gridCellSize * dimensions.z
        ));

        if (!Application.isPlaying) return;

        Gizmos.color = new Vector4(1f,1f,0f,0.5f);
        for(int o = 0; o < obstacles.Length; o++) { 
            Obstacle obs = obstacles[o];
            //Vector3Int currentCenterGridIndices = Boid3DHelpers.GetGridXYZIndices(numCellsPerAxis, transform.position, gridCellSize, obs.obstacle.transform.position);
            //for(int j = 0; j < obs.relativeGridIndices.Length; j++) {
            //    Gizmos.DrawSphere(
            //        Boid3DHelpers.GetGridCellWorldPositionFromXYZIndices(
            //            numCellsPerAxis, 
            //            transform.position, 
            //            gridCellSize, 
            //            obs.relativeGridIndices[j] + currentCenterGridIndices
            //        ),
            //        0.1f
            //    );
            //}
            for(int t = 0; t < obs.obstacle.sharedMesh.triangles.Length; t+=3) {
                Plane plane = new Plane(
                    obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t]],
                    obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+1]],
                    obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+2]]
                );
                for(int i = 0; i <= 2; i++) {
                    for (int j = 0; j <= 2; j++) {
                        if(i==j) continue;
                        //Gizmos.color = Color.white;
                        //Gizmos.DrawLine(
                        //    obs.obstacle.transform.position + obs.obstacle.sharedMesh.vertices[obs.obstacle.sharedMesh.triangles[t+i]], 
                        //    obs.obstacle.transform.position + obs.obstacle.sharedMesh.vertices[obs.obstacle.sharedMesh.triangles[t+j]]
                        //);
                        float xDiff = obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+j]].x - obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+i]].x;
                        float yDiff = obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+j]].y - obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+i]].y;
                        float zDiff = obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+j]].z - obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+i]].z;
                        float xSign = Mathf.Sign(xDiff);
                        float ySign = Mathf.Sign(yDiff);
                        float zSign = Mathf.Sign(zDiff);
                        float xDist = Mathf.Abs(xDiff);
                        float yDist = Mathf.Abs(yDiff);
                        float zDist = Mathf.Abs(zDiff);
                        int numXSteps = (int)(xDist / gridCellSize);
                        int numYSteps = (int)(yDist / gridCellSize);
                        int numZSteps = (int)(zDist / gridCellSize);

                        for(int x = 0; x <= numXSteps ; x++) {
                            for(int y = 0; y <= numYSteps; y++) {
                                for(int z = 0; z <= numZSteps; z++) {
                                    Vector3 tempPos = new Vector3(
                                        obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+i]].x + (x*gridCellSize*xSign),
                                        obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+i]].y + (y*gridCellSize*ySign),
                                        obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+i]].z + (z*gridCellSize*zSign)
                                    );
                                    Vector3 closestPointToLine = FindNearestPointOnLine(
                                        obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+i]], 
                                        obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+j]], 
                                        tempPos
                                    );
                                    Vector3 closestPointToPlane = plane.ClosestPointOnPlane(tempPos);
                                    float distToLine = Vector3.Distance(tempPos, closestPointToLine);
                                    float distToPlane = plane.GetDistanceToPoint(tempPos);
                                    //Gizmos.color = Color.white;
                                    //Gizmos.DrawSphere(closestPointToPlane, 0.05f);
                                    //Gizmos.color = Color.yellow;
                                    //Gizmos.DrawSphere(closestPointToLine, 0.05f);
                                    if (
                                        distToLine <= gridCellDiagonalHalf
                                        || (
                                            PointInTriangle(
                                                closestPointToPlane, 
                                                obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t]], 
                                                obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+1]], 
                                                obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+2]]
                                            ) 
                                            && Mathf.Abs(distToPlane) <= gridCellDiagonalHalf
                                        )
                                    ) {
                                        Gizmos.color = new Vector4(0f,0f,1f,0.1f);
                                        Gizmos.DrawCube(tempPos, gridCellSize3D); 
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        /*
        if (testMesh == null) return;

        Mesh mesh = testMesh.sharedMesh;
        particlePositions = new Vector3[mesh.vertices.Length];
        Gizmos.color = cellColor;
        for(int i = 0; i < mesh.vertices.Length; i++) {
            Vector3 gridPos = Boid3DHelpers.GetGridCellWorldPositionFromGivenPosition(
                dimensions, 
                transform.position, 
                bounds, 
                cellSize, 
                testMesh.transform.position + mesh.vertices[i]
            );
            particlePositions[i] = gridPos;
            Gizmos.DrawCube(particlePositions[i], cellSize3D);
        }
        
        
        for(int t = 0; t < mesh.triangles.Length; t += 3) {
            Plane plane = new Plane(
                particlePositions[mesh.triangles[t]], 
                particlePositions[mesh.triangles[t+1]], 
                particlePositions[mesh.triangles[t+2]]
            );
            for(int i = 0; i <= 2; i++) {
                for (int j = 0; j <= 2; j++) {
                    if(i==j) continue;
                    Gizmos.color = Color.white;
                    Gizmos.DrawLine(particlePositions[mesh.triangles[t+i]], particlePositions[mesh.triangles[t+j]]);
                    if (!Application.isPlaying) continue;

                    float xDiff = particlePositions[mesh.triangles[t+j]].x - particlePositions[mesh.triangles[t+i]].x;
                    float yDiff = particlePositions[mesh.triangles[t+j]].y - particlePositions[mesh.triangles[t+i]].y;
                    float zDiff = particlePositions[mesh.triangles[t+j]].z - particlePositions[mesh.triangles[t+i]].z;
                    float xSign = Mathf.Sign(xDiff);
                    float ySign = Mathf.Sign(yDiff);
                    float zSign = Mathf.Sign(zDiff);
                    float xDist = Mathf.Abs(xDiff);
                    float yDist = Mathf.Abs(yDiff);
                    float zDist = Mathf.Abs(zDiff);
                    int numXSteps = (int)(xDist / cellSize);
                    int numYSteps = (int)(yDist / cellSize);
                    int numZSteps = (int)(zDist / cellSize);
                    for(int x = 0; x <= numXSteps ; x++) {
                        for(int y = 0; y <= numYSteps; y++) {
                            for(int z = 0; z <= numZSteps; z++) {
                                Vector3 tempPos = new Vector3(
                                    particlePositions[mesh.triangles[t+i]].x + (x*cellSize*xSign),
                                    particlePositions[mesh.triangles[t+i]].y + (y*cellSize*ySign),
                                    particlePositions[mesh.triangles[t+i]].z + (z*cellSize*zSign)
                                );
                                Vector3 closestPointToLine = FindNearestPointOnLine(particlePositions[mesh.triangles[t+i]], particlePositions[mesh.triangles[t+j]], tempPos);
                                Vector3 closestPointToPlane = plane.ClosestPointOnPlane(tempPos);
                                float distToLine = Vector3.Distance(tempPos, closestPointToLine);
                                float distToPlane = plane.GetDistanceToPoint(tempPos);
                                //Gizmos.color = Color.white;
                                //Gizmos.DrawSphere(closestPointToPlane, 0.05f);
                                //Gizmos.color = Color.yellow;
                                //Gizmos.DrawSphere(closestPointToLine, 0.05f);
                                if (
                                    distToLine <= cellDiagonalHalf
                                    || (
                                        PointInTriangle(closestPointToPlane, particlePositions[mesh.triangles[t]], particlePositions[mesh.triangles[t+1]], particlePositions[mesh.triangles[2]]) 
                                        && Mathf.Abs(distToPlane) <= cellDiagonalHalf
                                    )
                                ) {
                                    Gizmos.color = new Vector4(0f,0f,1f,0.1f);
                                    Gizmos.DrawCube(tempPos, cellSize3D); 
                                }
                            }
                        }
                    }
                }
            }
        }
        */
    }

    private void Start() {
        InitializeGrid();
        for(int i = 0; i < obstacles.Length; i++) {
            PreprocessObstacle(ref obstacles[i]);
        }
    }

    private void InitializeGrid() {
        // Determine cell properties
        gridCellSize3D = new Vector3(gridCellSize, gridCellSize, gridCellSize);
        gridCellDiagonalHalf = Mathf.Pow(3f,0.5f) * gridCellSize * 0.5f;
        // We determine grid properties.
        // 1) we calculate how many cells we can get
        numCellsPerAxis = new Vector3Int(
            Mathf.CeilToInt(bounds.x/gridCellSize)+numBufferCellsPerAxis.x,
            Mathf.CeilToInt(bounds.y/gridCellSize)+numBufferCellsPerAxis.y,
            Mathf.CeilToInt(bounds.z/gridCellSize)+numBufferCellsPerAxis.z
        );
        numGridCells = numCellsPerAxis.x * numCellsPerAxis.y * numCellsPerAxis.z;

        if (verbose) {
            Debug.Log($"Number of Cells Per Axis: {numCellsPerAxis}");
            Debug.Log($"Total number of cells: {numGridCells}");
        }
    }

    private void PreprocessObstacle(ref Obstacle obs) {
        Mesh mesh = obs.obstacle.mesh;
        obs.vertexGridWorldPositions = new Vector3[mesh.vertices.Length];
        obs.currentCenterGridIndices = Boid3DHelpers.GetGridXYZIndices(numCellsPerAxis, transform.position, gridCellSize, obs.obstacle.transform.position);
        obs.currentRotation = obs.obstacle.transform.rotation;
        obs.currentScale = obs.obstacle.transform.localScale;
        
        for(int i = 0; i < mesh.vertices.Length; i++) {
            // Get the world space position of the current vertex. By default, vertices are local to the obstacle.
            Vector3 vertexCurrentWorldPos = obs.obstacle.transform.TransformPoint(mesh.vertices[i]);
            // We pre-populate the vertex's grid cell world position as well.
            obs.vertexGridWorldPositions[i] = Boid3DHelpers.GetGridCellWorldPositionFromGivenPosition(
                numCellsPerAxis,
                transform.position,
                gridCellSize,
                vertexCurrentWorldPos
            );
        }
    }

    private void Update() {
        for(int i = 0; i < obstacles.Length; i++) {
            Obstacle obs = obstacles[i];
            UpdateObstacle(ref obs);

        }
    }

    private void UpdateObstacle(ref Obstacle obs) {
        
        Vector3Int newCenterGridIndices = Boid3DHelpers.GetGridXYZIndices(numCellsPerAxis, transform.position, gridCellSize, obs.obstacle.transform.position);
        Quaternion newRotation = obs.obstacle.transform.rotation;
        Vector3 newScale = obs.obstacle.transform.localScale;

        if (
            newCenterGridIndices != obs.currentCenterGridIndices
            || newRotation != obs.currentRotation 
            || newScale != obs.currentScale
        ) {
            Debug.Log("UPDATING!");
            // The position has updated. We need to re-update the vertex grid world positions
            Mesh mesh = obs.obstacle.mesh;
            for(int i = 0; i < obs.vertexGridWorldPositions.Length; i++) {
                Vector3 vertexCurrentWorldPos = obs.obstacle.transform.TransformPoint(mesh.vertices[i]);
                obs.vertexGridWorldPositions[i] = Boid3DHelpers.GetGridCellWorldPositionFromGivenPosition(
                    numCellsPerAxis,
                    transform.position,
                    gridCellSize,
                    vertexCurrentWorldPos
                );
            }
            /*
            for(int t = 0; t < obs.obstacle.sharedMesh.triangles.Length; t+=3) {
                Plane plane = new Plane(
                    obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t]],
                    obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+1]],
                    obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+2]]
                );
                for(int i = 0; i <= 2; i++) {
                    for (int j = 0; j <= 2; j++) {
                        if(i==j) continue;
                        Vector3 avgNormal = (mesh.normals[])/3f;
                        float xDiff = obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+j]].x - obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+i]].x;
                        float yDiff = obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+j]].y - obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+i]].y;
                        float zDiff = obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+j]].z - obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+i]].z;
                        float xSign = Mathf.Sign(xDiff);
                        float ySign = Mathf.Sign(yDiff);
                        float zSign = Mathf.Sign(zDiff);
                        float xDist = Mathf.Abs(xDiff);
                        float yDist = Mathf.Abs(yDiff);
                        float zDist = Mathf.Abs(zDiff);
                        int numXSteps = (int)(xDist / gridCellSize);
                        int numYSteps = (int)(yDist / gridCellSize);
                        int numZSteps = (int)(zDist / gridCellSize);
                        for(int x = 0; x <= numXSteps ; x++) {
                            for(int y = 0; y <= numYSteps; y++) {
                                for(int z = 0; z <= numZSteps; z++) {
                                    Vector3 tempPos = new Vector3(
                                        obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+i]].x + (x*gridCellSize*xSign),
                                        obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+i]].y + (y*gridCellSize*ySign),
                                        obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+i]].z + (z*gridCellSize*zSign)
                                    );
                                    Vector3 closestPointToLine = FindNearestPointOnLine(
                                        obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+i]], 
                                        obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+j]], 
                                        tempPos
                                    );
                                    Vector3 closestPointToPlane = plane.ClosestPointOnPlane(tempPos);
                                    float distToLine = Vector3.Distance(tempPos, closestPointToLine);
                                    float distToPlane = plane.GetDistanceToPoint(tempPos);
                                    if (
                                        distToLine <= gridCellDiagonalHalf
                                        || (
                                            PointInTriangle(
                                                closestPointToPlane, 
                                                obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t]], 
                                                obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+1]], 
                                                obs.vertexGridWorldPositions[obs.obstacle.sharedMesh.triangles[t+2]]
                                            ) 
                                            && Mathf.Abs(distToPlane) <= gridCellDiagonalHalf
                                        )
                                    ) {
                                        //Gizmos.color = new Vector4(0f,0f,1f,0.1f);
                                        //Gizmos.DrawCube(tempPos, gridCellSize3D); 
                                    }
                                }
                            }
                        }
                    }
                }
            }
            */
        }
        
        obs.currentCenterGridIndices = newCenterGridIndices;
        obs.currentRotation = newRotation;
        obs.currentScale = newScale;
    }

    public Vector3 FindNearestPointOnLine(Vector3 origin, Vector3 end, Vector3 point) {
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

    public bool SameSide(Vector3 p1, Vector3 p2, Vector3 a, Vector3 b) {
        Vector3 cp1 = Vector3.Cross(b-a, p1-a);
        Vector3 cp2 = Vector3.Cross(b-a, p2-a);
        return Vector3.Dot(cp1, cp2) >= 0;
    }

    public bool PointInTriangle(Vector3 p, Vector3 a, Vector3 b, Vector3 c) {
        return SameSide(p,a, b,c) && SameSide(p,b, a,c) && SameSide(p,c, a,b);
    }
}
