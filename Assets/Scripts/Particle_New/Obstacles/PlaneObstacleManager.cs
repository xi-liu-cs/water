using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using System.Linq;

public class PlaneObstacleManager : Grid3D
{   
    // We store a global reference to this specific script via the `current` method.
    public static PlaneObstacleManager current;

    public List<PlaneObstacle> obstacles = new List<PlaneObstacle>();
    public List<float3> vertices = new List<float3>();
    public List<PlaneObstacle.ObsPlane> obstaclePlanes = new List<PlaneObstacle.ObsPlane>();
    public List<int> vertexTriangleMap = new List<int>();
    public List<int3> vertexTriangleCount = new List<int3>();

    public List<int> vertexBuckets;
    public int2[] vertexBucketCounts;
    public int maxVerticesPerBucket;

    public List<Transform> debugParticles = new List<Transform>();
    public bool awaitInitialization = false;

    public List<Vector3> closestVertices = new List<Vector3>();
    public List<Vector3> gizmoNormals = new List<Vector3>();

    public ComputeBuffer obstacleBuffer;
    public ComputeBuffer verticesBuffer;
    public ComputeBuffer planesBuffer;
    public ComputeBuffer vertexTriangleMapBuffer;
    public ComputeBuffer vertexTriangleCountBuffer;
    public ComputeBuffer vertexBucketsBuffer;
    public ComputeBuffer vertexBucketCountsBuffer;
    // Purely for debugging
    public ComputeBuffer debugParticlesBuffer;
    public ComputeBuffer particleIsIntersectingBuffer;
    public ComputeBuffer testWorldPosBuffer;
    // Kernels
    public int dampenParticlesByObstaclesKernel;
    public int testParticleLocalToWorldKernel;
    private int _BLOCK_SIZE = 1024;

    public ComputeShader obstacleShader;

    public bool drawGizmos = false;

    void OnDrawGizmos() {
        if (!drawGizmos) return;
        
        DrawGridAxes();
        DrawGridBounds();
        DrawGridInnerBounds();
        
        /*
        if (closestVertices.Count == 0) return;
        Gizmos.color = Color.blue;
        for(int i = 0; i < closestVertices.Count; i++) {
            Vector3 closestVertex = closestVertices[i];
            Vector3 closestNormal = gizmoNormals[i];
            Gizmos.DrawSphere(closestVertex, 0.1f);
            Gizmos.DrawRay(closestVertex, closestNormal);
        }
        */
        if(Application.isPlaying) {
            Gizmos.color = Color.blue;
            float3[] tempTestWorldPos = new float3[obstaclePlanes.Count];
            testWorldPosBuffer.GetData(tempTestWorldPos);

            int obsPlCount = 0;
            for(int i = 0; i < tempTestWorldPos.Length; i++) {
                if (tempTestWorldPos[i][0] != -100f && tempTestWorldPos[i][1] != -100f && tempTestWorldPos[i][2] != -100f) {
                    Gizmos.DrawSphere(
                        new Vector3(tempTestWorldPos[i][0], tempTestWorldPos[i][1], tempTestWorldPos[i][2]), 
                        0.1f
                    );
                    obsPlCount += 1;
                }
            }
            Debug.Log($"ObsPlCount: {obsPlCount}");
        }
    }

    public void PreprocessPlanes() {
        vertices = new List<float3>();
        vertexTriangleMap = new List<int>();
        obstaclePlanes = new List<PlaneObstacle.ObsPlane>();
        vertexTriangleCount = new List<int3>();
        for(int i = 0; i < obstacles.Count; i++) {
            PlaneObstacle obstacle = obstacles[i];
            obstacle.Initialize(i, vertices.Count);
            vertices.AddRange(obstacle.vertices);
            
            for(int j = 0; j < obstacle.vertexTriangleCount.Length; j++) {
                vertexTriangleCount.Add(new(
                    obstacle.vertexTriangleCount[j][0],
                    obstacle.vertexTriangleCount[j][1]+vertexTriangleMap.Count, 
                    obstacle.vertexTriangleCount[j][2])
                );
            }
            // We add these last because we need to update vertexTriangleCount first
            for(int j = 0; j < obstacle.vertexTriangleMap.Count; j++) {
                vertexTriangleMap.Add(obstacle.vertexTriangleMap[j] + obstaclePlanes.Count);
            }
            //vertexTriangleMap.AddRange(obstacle.vertexTriangleMap);
            obstaclePlanes.AddRange(obstacle.obstaclePlanes);
        }
        ReorderVertices();
    }

    public void ReorderVertices() {
        // Given a list of vertices, let's re-sort them so that they're in a prescribed order.
        // The order is defined by their projected 1D index. The ordering within these buckets doesn't matter.
        
        // Firstly, let's create our buckets. We need to determine how many cells there exactly are, then go from there.
        List<int>[] tempVertexBuckets = new List<int>[numGridCells];
        // Now, we iterate through all vertices and assign them to a bucket.
        for(int i = 0; i < vertices.Count; i++) {
            PlaneObstacle obstacle = obstacles[vertexTriangleCount[i][0]];
            Vector3 worldPos = LocalPointToWorldPoint(obstacle.obs, vertices[i]);
            // HOWEVER: if our point is outside our grid bounds, we won't consider them.
            if (
                worldPos.x < origin.x - bounds.x/2f || worldPos.x > origin.x + bounds.x/2f ||
                worldPos.y < origin.y - bounds.y/2f || worldPos.y > origin.y + bounds.y/2f ||
                worldPos.z < origin.z - bounds.z/2f || worldPos.z > origin.z + bounds.z/2f
            ) continue;
            int projectedIndex = Grid3DHelpers.GetProjectedGridIndexFromGivenPosition(
                numGridCellsPerAxis,
                origin,
                gridCellSize,
                worldPos
            );
            if (tempVertexBuckets[projectedIndex] == null) tempVertexBuckets[projectedIndex] = new List<int>();
            tempVertexBuckets[projectedIndex].Add(i);
        }

        for(int i = 0; i < numGridCells; i++) {
            if (tempVertexBuckets[i] == null) tempVertexBuckets[i] = new List<int>();
        }

        // Now we can condense the buckets into a single list!
        maxVerticesPerBucket = new List<List<int>>(tempVertexBuckets).Max(x => x.Count);
        vertexBuckets = new List<int>();
        vertexBucketCounts = new int2[numGridCells];
        for(int i = 0; i < numGridCells; i++) {
            if (tempVertexBuckets[i] == null) {
                vertexBucketCounts[i] = new(vertexBuckets.Count, 0);
                // there's nothing to add to the vertex buckets list...
            } else {
                vertexBucketCounts[i] = new(vertexBuckets.Count, tempVertexBuckets[i].Count);
                vertexBuckets.AddRange(tempVertexBuckets[i]);
            }
        }

        Debug.Log("Vertices Reordered!");
    }

    public void DebugCheckIfIntersecting() {
        if (debugParticles.Count == 0) {
            Debug.Log("Cannot check anything if the debug particle isn't set!");
            return;
        }

        // Debugging. Remove later
        closestVertices = new List<Vector3>();
        gizmoNormals = new List<Vector3>();

        // Initialize some variables
        PlaneObstacle.ObsPlane p;
        Vector3 worldCentroid;
        Vector3 projectionPoint;
        Vector3 targetVector;
        float dotBetweenParticleAndNormal, distanceToPlane;

        // 1) Determine all vertices close to the particle's current position
        int3[] cVertices = GetClosestVertexIndex(debugParticles[0].position);
        if (cVertices.Length == 0) {
            Debug.Log("We couldn't find any vertices close to our debug particle's position...");
            return;
        }
        
        // We need to iterate through all vertices in this case
        int isIntersecting = 0;
        PlaneObstacle obstacle;
        for(int i = 0; i < cVertices.Length; i++) {
            // If the returned value has a 1st int == -1, then we know it's a dud and all future ones are duds as well
            if (cVertices[i][0] == -1) break;
            // Get the obstacle associated with this vertex
            obstacle = obstacles[cVertices[i][0]];
            // Now iterate through all planes associated with this vertex
            for(int c = cVertices[i][1]; c < cVertices[i][1] + cVertices[i][2]; c++) {
                p = obstaclePlanes[vertexTriangleMap[c]];
                // For each plane `p`, we:
                // 1) convert its centroid to world position
                worldCentroid = LocalPointToWorldPoint(obstacle.obs, p.centroid);
                closestVertices.Add(worldCentroid);
                gizmoNormals.Add(LocalVectorToWorldVector(obstacle.obs, p.normalVector));
                // 2) Calculate dot product between normal vector and vector b/w centroid and particle
                targetVector = (debugParticles[0].position - worldCentroid).normalized;
                dotBetweenParticleAndNormal = Vector3.Dot(targetVector, LocalVectorToWorldVector(obstacle.obs, p.normalVector));
                // 3) Calculate the projection of the debug particle onto the plane
                projectionPoint = ClosestPointOnPlane(
                    worldCentroid, 
                    LocalVectorToWorldVector(obstacle.obs, p.normalVector), 
                    debugParticles[0].position
                );
                //closestVertices.Add(projectionPoint);
                // 4) Calculate distance from the debug particle to its projection
                distanceToPlane = DistanceFromPlane(worldCentroid, LocalVectorToWorldVector(obstacle.obs, p.normalVector), debugParticles[0].position);
                // 5) Finally check the intersection
                if(dotBetweenParticleAndNormal <= 0f && ObstacleHelper.PointInTriangle(
                        projectionPoint, 
                        LocalPointToWorldPoint(obstacle.obs, p.vertex1), 
                        LocalPointToWorldPoint(obstacle.obs, p.vertex2), 
                        LocalPointToWorldPoint(obstacle.obs, p.vertex3)
                )) isIntersecting += 1;
            }
        }

        // Let us know if we're intersecting!
            if (isIntersecting > 0)  Debug.Log("The particle is intersecting with one or more obstacles!");
            else                 Debug.Log("The particle isn't intersecting with an obstacle...");
    } 
    
    private void Awake() {
        current = this;
        if (!awaitInitialization) Initialize();
    }

    public void Initialize() {
        InitializeKernels();
        InitializeShaderVariables();
        InitializeBuffers();
    }

    private void InitializeKernels() {
        dampenParticlesByObstaclesKernel = obstacleShader.FindKernel("DampenParticlesByObstacles");
        testParticleLocalToWorldKernel = obstacleShader.FindKernel("TestParticleLocalToWorld");
    }

    private void InitializeShaderVariables() {
        obstacleShader.SetInt("numParticles",debugParticles.Count);
        obstacleShader.SetInt("numBoids",0);
        obstacleShader.SetInt("numVertices",vertices.Count);
        obstacleShader.SetInt("numPlanes", obstaclePlanes.Count);
        UpdateShaderVariables();
    }

    private void InitializeBuffers() {
        // Initialize obstacles buffer
        PlaneObstacle.Obs[] obstaclesArray = new PlaneObstacle.Obs[obstacles.Count];
        for(int i = 0; i < obstacles.Count; i++) obstaclesArray[i] = obstacles[i].obs;
        obstacleBuffer = new ComputeBuffer(obstacles.Count, sizeof(int)*4 + sizeof(float)*10);
        obstacleBuffer.SetData(obstaclesArray);

        // Initialize vertices buffer
        verticesBuffer = new ComputeBuffer(vertices.Count, sizeof(float)*3);
        verticesBuffer.SetData(vertices.ToArray());

        // Initialize planes buffer
        planesBuffer = new ComputeBuffer(obstaclePlanes.Count, sizeof(int) + sizeof(float)*15);
        planesBuffer.SetData(obstaclePlanes.ToArray());

        // Initialize Vertex-Triangle-Map
        vertexTriangleMapBuffer = new ComputeBuffer(vertexTriangleMap.Count, sizeof(int));
        vertexTriangleMapBuffer.SetData(vertexTriangleMap.ToArray());

        // Initialize Vertex-Triangle-Count
        vertexTriangleCountBuffer = new ComputeBuffer(vertexTriangleCount.Count, sizeof(int)*3);
        vertexTriangleCountBuffer.SetData(vertexTriangleCount.ToArray());

        // Initialize vertex buckets
        vertexBucketsBuffer = new ComputeBuffer(vertexBuckets.Count, sizeof(int));
        vertexBucketsBuffer.SetData(vertexBuckets.ToArray());

        // Initialize vertex bucket counts
        vertexBucketCountsBuffer = new ComputeBuffer(vertexBucketCounts.Length, sizeof(int)*2);
        vertexBucketCountsBuffer.SetData(vertexBucketCounts);

        // Purely for debugging
        particleIsIntersectingBuffer = new ComputeBuffer(debugParticles.Count, sizeof(int));
        debugParticlesBuffer = new ComputeBuffer(debugParticles.Count, sizeof(float)*3);
        float3[] debugParticlesF = new float3[debugParticles.Count];
        for(int i = 0; i < debugParticles.Count; i++) {
            debugParticlesF[i] = new(
                debugParticles[i].position.x,
                debugParticles[i].position.y,
                debugParticles[i].position.z
            );
        }
        debugParticlesBuffer.SetData(debugParticlesF);

        testWorldPosBuffer = new ComputeBuffer(obstaclePlanes.Count, sizeof(float)*3);

        // Now to set to kernel
        obstacleShader.SetBuffer(dampenParticlesByObstaclesKernel, "particles", debugParticlesBuffer);
        obstacleShader.SetBuffer(dampenParticlesByObstaclesKernel, "obstacles", obstacleBuffer);
        obstacleShader.SetBuffer(dampenParticlesByObstaclesKernel, "obstaclePlanes", planesBuffer);
        obstacleShader.SetBuffer(dampenParticlesByObstaclesKernel, "vertexBucketCounts", vertexBucketCountsBuffer);
        obstacleShader.SetBuffer(dampenParticlesByObstaclesKernel, "vertexBuckets", vertexBucketsBuffer);
        obstacleShader.SetBuffer(dampenParticlesByObstaclesKernel, "vertexTriangleCount", vertexTriangleCountBuffer);
        obstacleShader.SetBuffer(dampenParticlesByObstaclesKernel, "vertexTriangleMap", vertexTriangleMapBuffer);
        //obstacleShader.SetBuffer(dampenParticlesByObstaclesKernel, "particleIsIntersecting", particleIsIntersectingBuffer);

        obstacleShader.SetBuffer(dampenParticlesByObstaclesKernel, "testWorldPos", testWorldPosBuffer);

        obstacleShader.SetBuffer(testParticleLocalToWorldKernel, "obstacles", obstacleBuffer);
        obstacleShader.SetBuffer(testParticleLocalToWorldKernel, "obstacleVertices", verticesBuffer);
        obstacleShader.SetBuffer(testParticleLocalToWorldKernel, "obstaclePlanes", planesBuffer);
        obstacleShader.SetBuffer(testParticleLocalToWorldKernel, "vertexBucketCounts", vertexBucketCountsBuffer);
        obstacleShader.SetBuffer(testParticleLocalToWorldKernel, "vertexBuckets", vertexBucketsBuffer);
        obstacleShader.SetBuffer(testParticleLocalToWorldKernel, "vertexTriangleCount", vertexTriangleCountBuffer);
        obstacleShader.SetBuffer(testParticleLocalToWorldKernel, "vertexTriangleMap", vertexTriangleMapBuffer);
        obstacleShader.SetBuffer(testParticleLocalToWorldKernel, "testWorldPos", testWorldPosBuffer);

    }

    private void Update() {
        if (!awaitInitialization) CustomUpdate();
    }

    private void UpdateShaderVariables() {
        // Grid-related changes
        obstacleShader.SetFloats("origin", originF);
        obstacleShader.SetFloats("numCellsPerAxis", numGridCellsPerAxisF);
        obstacleShader.SetFloat("gridCellSize", gridCellSize);
        obstacleShader.SetFloats("bounds", boundsF);

        // Vertex-related
        obstacleShader.SetInt("maxVerticesPerCell", maxVerticesPerBucket);
    }

    private void UpdateBuffers() {
        // The only things we need to update are the obstacles buffer and, if any obstacles were changed, that we re-order the vertices
        PlaneObstacle.Obs[] obstaclesArray = new PlaneObstacle.Obs[obstacles.Count];
        bool needToReorder = false;
        for(int i = 0; i < obstacles.Count; i++) {
            obstaclesArray[i] = obstacles[i].obs;
            needToReorder |= obstacles[i].obs.hasChanged == 1;
        }
        obstacleBuffer.SetData(obstaclesArray);
        
        if (needToReorder) {
            Debug.Log("We need to reorder...");
            ReorderVertices();

            // Vertex-related
            obstacleShader.SetInt("maxVerticesPerCell", maxVerticesPerBucket);

            // Update vertex buckets
            vertexBucketsBuffer.Release();
            vertexBucketsBuffer = new ComputeBuffer(vertexBuckets.Count, sizeof(int));
            vertexBucketsBuffer.SetData(vertexBuckets.ToArray());

            // Update vertex bucket counts
            vertexBucketCountsBuffer.Release();
            vertexBucketCountsBuffer = new ComputeBuffer(vertexBucketCounts.Length, sizeof(int)*2);
            vertexBucketCountsBuffer.SetData(vertexBucketCounts);

            obstacleShader.SetBuffer(dampenParticlesByObstaclesKernel, "vertexBucketCounts", vertexBucketCountsBuffer);
            obstacleShader.SetBuffer(dampenParticlesByObstaclesKernel, "vertexBuckets", vertexBucketsBuffer);

            obstacleShader.SetBuffer(testParticleLocalToWorldKernel, "vertexBucketCounts", vertexBucketCountsBuffer);
            obstacleShader.SetBuffer(testParticleLocalToWorldKernel, "vertexBuckets", vertexBucketsBuffer);
        }

        if (!awaitInitialization) {
            // We're doing this with our own particles
            float3[] debugParticlesF = new float3[debugParticles.Count];
            for(int i = 0; i < debugParticles.Count; i++) {
                debugParticlesF[i] = new(
                    debugParticles[i].position.x,
                    debugParticles[i].position.y,
                    debugParticles[i].position.z
                );
            }
            debugParticlesBuffer.SetData(debugParticlesF);
        }
    }

    public void CustomUpdate() {
        // We have to update the compute shader
        UpdateShaderVariables();
        UpdateBuffers();

        obstacleShader.Dispatch(testParticleLocalToWorldKernel, Mathf.CeilToInt((float)obstaclePlanes.Count / (float)_BLOCK_SIZE),1,1);

        // Now, we can perform a dispatch
        int n = Mathf.CeilToInt((float)debugParticles.Count / (float)_BLOCK_SIZE);
        obstacleShader.Dispatch(dampenParticlesByObstaclesKernel, n,1,1);

        /*
        // To check our status, we grab the data from the `particleIsIntsersecting` buffer
        int[] areParticlesIntersecting = new int[debugParticles.Count];
        particleIsIntersectingBuffer.GetData(areParticlesIntersecting);
        for(int i = 0; i < debugParticles.Count; i++) {
            Debug.Log($"Particle {i+1} is... {areParticlesIntersecting[i]}");
        }
        */

       
        
        /*
        int2[] tempVertexBucketCounts = new int2[5];
        vertexBucketCountsBuffer.GetData(tempVertexBucketCounts);
        Debug.Log(tempVertexBucketCounts[0]);
        */
        //DebugCheckIfIntersecting();
    }

    private void OnDestroy() {
        ClearBuffers();
    }

    public void ClearBuffers() {
        obstacleBuffer.Release();
        verticesBuffer.Release();
        planesBuffer.Release();
        vertexTriangleMapBuffer.Release();
        vertexTriangleCountBuffer.Release();
        vertexBucketsBuffer.Release();
        vertexBucketCountsBuffer.Release();
        debugParticlesBuffer.Release();
    }

    private int GetClosestVertexIndex(PlaneObstacle obstacle, Vector3 pos) {
        // Initialize some variables
        int closestIndex = -1;
        float closestDistance = Mathf.Infinity;
        float dist;
        Vector3 worldVertex;
        for(int i = obstacle.obs.verticesStart; i < obstacle.obs.verticesStart + obstacle.obs.verticesCount; i++) {
            // Calculate world position of this vertex
            worldVertex = LocalPointToWorldPoint(obstacle.obs, vertices[i]);
            // Get distance
            dist = Vector3.Distance(pos, worldVertex);
            // If closer, replace
            if (dist < closestDistance) {
                closestIndex = i;
                closestDistance = dist;
            }
        }
        return closestIndex;
    }

    private int3[] GetClosestVertexIndex(Vector3 pos) {
        float closestDistance = Mathf.Infinity;
        //closestObstacle = obstacles[vertexTriangleCount[0][0]];
        
        float dist;
        Vector3 worldVertex;
        PlaneObstacle obstacle;

        // return early if our position is outside the bounds of the grid
        if (
            pos.x < origin.x - bounds.x/2f || pos.x > origin.x + bounds.x/2f ||
            pos.y < origin.y - bounds.y/2f || pos.y > origin.y + bounds.y/2f ||
            pos.z < origin.z - bounds.z/2f || pos.z > origin.z + bounds.z/2f 
        ) return new int3[0];

        // Get the projected position of our particle
        Vector3Int xyz = Grid3DHelpers.GetGridXYZIndices(numGridCellsPerAxis, origin, gridCellSize, pos);
        int minX = Mathf.Max(0,xyz.x - 1);
        int maxX = Mathf.Min(numGridCellsPerAxis.x-1, xyz.x+1);
        int minY = Mathf.Max(0, xyz.y - 1);
        int maxY = Mathf.Min(numGridCellsPerAxis.y-1, xyz.y+1);
        int minZ = Mathf.Max(0,xyz.z - 1);
        int maxZ = Mathf.Min(numGridCellsPerAxis.z-1, xyz.z+1);

        int maxPossibleVertices = (maxX - minX) * (maxY-minY) * (maxZ-minZ) * maxVerticesPerBucket;
        int3[] closestVertices = new int3[maxPossibleVertices];
        int closestVerticesCount = 0;

        int projectedIndex;
        for(int x = minX; x <= maxX; x++) {
            for(int y = minY; y <= maxY; y++) {
                for(int z = minZ; z <= maxZ; z++) {
                    projectedIndex = Grid3DHelpers.GetProjectedGridIndexFromXYZ(numGridCellsPerAxis,x,y,z);
                    if (vertexBucketCounts[projectedIndex][1] == 0) continue;
                    // We'll iterate through 
                    for(int i = vertexBucketCounts[projectedIndex][0]; i < vertexBucketCounts[projectedIndex][0] + vertexBucketCounts[projectedIndex][1]; i++) {
                        // Get the vertex index from vertexbuckets
                        int index = vertexBuckets[i];
                        // Add this vertex triangle count to `closestvertices` and increment count
                        closestVertices[closestVerticesCount] = vertexTriangleCount[index];
                        closestVerticesCount += 1;
                        /*
                        // Get the obstacle associated with this vertex
                        obstacle = obstacles[vertexTriangleCount[index][0]];
                        // Calculate world position of this vertex
                        worldVertex = LocalPointToWorldPoint(obstacle.obs, vertices[index]);
                        // Get distance
                        dist = Vector3.Distance(pos, worldVertex);
                        // If closer, replace
                        if (dist < closestDistance) {
                            closestIndex = index;
                            closestDistance = dist;
                            closestObstacle = obstacle;
                        }
                        */
                    }
                }
            }
        }

        /*
        int projectedIndex = Grid3DHelpers.GetProjectedGridIndexFromGivenPosition(numGridCellsPerAxis, origin, gridCellSize, pos);

        // Return early if we don't have any vertices either
        if (vertexBucketCounts[projectedIndex][1] == 0) return -1;

        // We'll iterate through 
        for(int i = vertexBucketCounts[projectedIndex][0]; i < vertexBucketCounts[projectedIndex][0] + vertexBucketCounts[projectedIndex][1]; i++) {
            // Get the vertex index from vertexbuckets
            int index = vertexBuckets[i];
            // Get the obstacle associated with this vertex
            obstacle = obstacles[vertexTriangleCount[index][0]];
            // Calculate world position of this vertex
            worldVertex = LocalPointToWorldPoint(obstacle.obs, vertices[index]);
            // Get distance
            dist = Vector3.Distance(pos, worldVertex);
            // If closer, replace
            if (dist < closestDistance) {
                closestIndex = index;
                closestDistance = dist;
                closestObstacle = obstacle;
            }
        }
        */

        // Make sure to fill remainders with (-1,-1,-1) to ensure that we don't accidently do something stupid
        for(int j = closestVerticesCount; j < maxPossibleVertices; j++) {
            closestVertices[j] = new(-1,-1,-1);
        }

        // Finally return
        return closestVertices;

        //return closestIndex;
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
    public static float3 LocalPointToWorldPoint(PlaneObstacle.Obs obstacle, Vector3 localPoint) {
        Vector3 s = new Vector3(localPoint.x * obstacle.scale[0], localPoint.y * obstacle.scale[1], localPoint.z * obstacle.scale[2]);
        Vector3 rotvec = RotMultVec3(obstacle.rotation, s);
        return new(rotvec.x + obstacle.position[0], rotvec.y + obstacle.position[1], rotvec.z + obstacle.position[2]);
    }
    public static Vector3 LocalPointToWorldPoint(PlaneObstacle.Obs obstacle, float3 localPoint) {
        Vector3 s = new Vector3(localPoint[0] * obstacle.scale[0], localPoint[1] * obstacle.scale[1], localPoint[2] * obstacle.scale[2]);
        Vector3 rotvec = RotMultVec3(obstacle.rotation, s);
        return new Vector3(rotvec.x + obstacle.position[0], rotvec.y + obstacle.position[1], rotvec.z + obstacle.position[2]);
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
    public static Vector3 RotMultVec3(float4 quat, Vector3 vec){
        float num = quat[0] * 2f;
        float num2 = quat[1] * 2f;
        float num3 = quat[2] * 2f;
        float num4 = quat[0] * num;
        float num5 = quat[1] * num2;
        float num6 = quat[2] * num3;
        float num7 = quat[0] * num2;
        float num8 = quat[0] * num3;
        float num9 = quat[1] * num3;
        float num10 = quat[3] * num;
        float num11 = quat[3] * num2;
        float num12 = quat[3] * num3;
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
    public static Vector3 LocalVectorToWorldVector(PlaneObstacle.Obs obstacle, Vector3 v) {
        Vector3 start = LocalPointToWorldPoint(obstacle, Vector3.zero);
        Vector3 end = LocalPointToWorldPoint(obstacle, v);
        return end - start;
    }
    public static Vector3 LocalVectorToWorldVector(PlaneObstacle.Obs obstacle, float3 v) {
        Vector3 start = LocalPointToWorldPoint(obstacle, Vector3.zero);
        Vector3 end = LocalPointToWorldPoint(obstacle,new Vector3(v[0],v[1],v[2]));
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


    /*
    

    // We store all obstacles inside of this list. We'll reference this list many times.
    // Note: we store both static and dynamic obstacles in this script, so we need to filter through them during each update.
    public List<PlaneObstacle> obstacles = new List<PlaneObstacle>();
    



    
    private List<PlaneObstacle.ObsPlane> obstaclePlanes;
    private PlaneObstacle.Obs[] obs;
    private int[] obsPlaneMap;

    // These are compute buffers that'll be passed onto the particle GPU
    private ComputeBuffer obstaclesBuffer;
    private ComputeBuffer obstaclePlanesBuffer;
    private ComputeBuffer obstaclePlanesMapBuffer;

    // Stores the # of obstacle planes in total
    [ReadOnly, SerializeField] private int numPlanes;

    public bool waitForInitialization = false;

    private void Awake() {
        current = this;
        if (!waitForInitialization) Initialize();
    }

    public void Initialize() {

        UpdateBuffers();
    }
    
    public void UpdateBuffers() {
        // Update obstacles buffer
        obs = new PlaneObstacle.Obs[obstacles.Count];
        obstaclePlanes = new List<PlaneObstacle.ObsPlane>();

        int numPlanes = 0;
        for(int i = 0; i < obstacles.Count; i++) {
            numPlanes += obstacles[i].obstaclePlanes.Count;
        }
        obsPlaneMap = new int[numPlanes];
        
        int start = 0;
        for(int i = 0; i < obstacles.Count; i++) {
            obs[i] = obstacles[i].obs;
            obstaclePlanes.AddRange(obstacles[i].obstaclePlanes);
            for(int j = start; j < obstaclePlanes.Count; j++) {
                obsPlaneMap[j] = i;
                start = j;
            }
            
        }
        obstaclesBuffer = new ComputeBuffer(obs.Length, sizeof(float)*10);
        obstaclesBuffer.SetData(obs);
        obstaclePlanesBuffer = new ComputeBuffer(obstaclePlanes.Count, sizeof(float)*15);
        obstaclePlanesBuffer.SetData(obstaclePlanes.ToArray());
        obstaclePlanesMapBuffer = new ComputeBuffer(obsPlaneMap.Length, sizeof(int));
        obstaclePlanesMapBuffer.SetData(obsPlaneMap);

        numPlanes = obstaclePlanes.Count;
    }

    void OnDestroy() {
        DisposeBuffers();
    }
    public void DisposeBuffers() {
        obstaclesBuffer.Dispose();
        obstaclePlanesBuffer.Dispose();
        obstaclePlanesMapBuffer.Dispose();
    }
    */
}
