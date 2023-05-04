using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Random = UnityEngine.Random;
using System.Runtime.InteropServices;

public class BoidManager_GPU : MonoBehaviour
{
    public struct Boid {
        public float px;
        public float py;
        public float pz;

        public float vx;
        public float vy;
        public float vz;

        public Vector3 position;
        public Vector3Int hashPosition;
        public int projectedHashPosition;
        public Vector3Int oldHashPosition;
        public int oldProjectedHashPosition;
        public Vector3 velocity;

        public int minX;
        public int maxX;
        public int minY;
        public int maxY;
        public int minZ;
        public int maxZ;

        public Vector3 separationInfluence;
        public Vector3 alignmentInfluence;
        public Vector3 cohesionInfluence;
        //public Vector3 targetInfluence;
        public Vector3 edgeInfluence;

        public int numNeighbors;
        public int numCloseNeighbors;

        public int targetIndex;
        public int numUpdates;
    };

    struct GridCell {
        public Vector3 separation;
        public Vector3 alignment;
        public Vector3 cohesion;
        public Vector3 position;
        public int n;
    };

    // Private list of structs that represent targets. Sent into a custom targets buffer
    private struct Target {
        public Vector3 position;
        public int numBoids;
        public int isActive;
    }

    [System.Serializable]
    public class Obstacle {
        public MeshFilter obstacle;
        public bool isStatic;
    }
    
    [Header("= CORE SETTINGS =")]
    [Tooltip("Total # of boids in the simulation")] 
        public int numBoids = 32;
    [Tooltip("The minimum limits of the world space boundaries of the simulation. Stored in XYZ format")] 
        public float[] minBounds = { -100f,-100f,-100f };
    [Tooltip("The maximum limits of the world space boundaries of the simulation. Stored in XYZ format")] 
        public float[] maxBounds = { 100f,100f,100f };
    [Tooltip("The inner XYZ grid dimension where boids are restricted to stay in")]
        public float[] innerMinBounds = { 2f, 2f, 2f };
    [Tooltip("The inner XYZ grid dimension where boids are restricted to stay in")]
        public float[] innerMaxBounds = { 8f, 8f, 8f };

    // The world space size of the simulation
    private float[] dimensions;
    [Tooltip("The # of cells along the X, Y, and Z axes each for the BOID grid")]
        public int[] boidGridDiscretizations = { 10,10,10 };
    [Tooltip("The # of cells along the X, Y, and Z axis each for the OBSTACLE grid")]
        public int[] obstacleGridDiscretizations = {20,20,20};
    // Tracking the size of each grid cell based on discretizations
    private float[] boidGridCellSize;
    private float[] obstacleGridSize;
    [Tooltip("The Compute Shader for Boid Behavior")]
        public ComputeShader boidShader;
    
    [SerializeField] private float updateTimeDelay = 0.1f;
    [Header("= BOID CONFIGURATIONS =")]
    [Tooltip("The visual range in XYZ grid dimensions that the boid is allowed to see")]
        public int visualRange = 5;
    [Tooltip("The range in XYZ grid dimensions that the boid considers to be too close")]
        public int closeRange = 0;
    [Tooltip("Minimum speed a boid can move")]
        public float minSpeed = 0.1f;
    [Tooltip("Maximum speed a boid can move")]
        public float maxSpeed = 2f;

    [Header("= WEIGHT CONFIGURATIONS =")]
    [Tooltip("Weight for cohesion between boids")] 
        public float cohesionFactor = 0.0005f;
    [Tooltip("Weight for alignment between boids")] 
        public float alignmentFactor = 0.05f;
    [Tooltip("Weight for separation from close boids")]
        public float separationFactor = 0.005f;
    // [Tooltip("Weight for avoiding obstacles. CURRENTLY UNUSED!")]
    //    public float obstacleFactor = 0.05f;
    [Tooltip("Weight for turning back into the boundary space")]
        public float turnFactor = 0.05f;
    [Tooltip("Weight for turning towards a target")]
        public float targetFactor = 0.75f;

    [Header("= TARGETS AND OBSTACLES =")]
    //[Tooltip("All gameobjects that are considered targets")]
        //public Transform[] targets = new Transform[0];
    [Tooltip("The total number of boids that a target can be associated with")]
        public int boidsPerTarget = 200;
    [Tooltip("All game objects (that have meshes attached to them) that are considered obstacles")]
        public Obstacle[] obstacles = new Obstacle[0];

    [Header("= DEBUG TOOLS =")]
    [Tooltip("Visualize the number of boids in each grid cell")]
        public bool visualizeBoidFrequency = false;
    [Tooltip("What color should grid cells be visualized with?")]
        public Color visualizeBoidFrequencyColor = Color.red;

    //private Target[] _targets;
    private Boid[] boids;
    private GridCell[] boidsGrid;
    
    private int boidUpdateNeighborCellsKernel;
    private int boidUpdateInfluenceVectorsKernel;
    private int boidUpdateStatusKernel;
    //private int updateKernel;

    private ComputeBuffer boidsBuffer;
    private ComputeBuffer boidsGridBuffer;
    //private ComputeBuffer targetsBuffer;

    void OnDrawGizmos() {
        Gizmos.color = Color.white;
        Vector3 minBoundVector = new Vector3(minBounds[0], minBounds[1], minBounds[2]);
        Vector3 d = new Vector3(
            maxBounds[0]-minBounds[0], 
            maxBounds[1]-minBounds[1], 
            maxBounds[2]-minBounds[2]
        );
        Gizmos.DrawWireCube(minBoundVector + d/2f, d);

        Gizmos.color = Color.red;
        Vector3 minInnerBoundVector = new Vector3(innerMinBounds[0], innerMinBounds[1], innerMinBounds[2]);
        Vector3 d2 = new Vector3(
            innerMaxBounds[0]-innerMinBounds[0], 
            innerMaxBounds[1]-innerMinBounds[1], 
            innerMaxBounds[2]-innerMinBounds[2]
        );
        Gizmos.DrawWireCube(minInnerBoundVector + d2/2f, d2);

        if (Application.isPlaying) {
            
            int bCount = boidsBuffer.count;
            Boid[] _boids = new Boid[bCount];
            boidsBuffer.GetData(_boids);

            GridCell[] visualGrid = new GridCell[boidsGridBuffer.count];
            boidsGridBuffer.GetData(visualGrid);
            
            for(int i = 0; i < bCount; i++) {
                Gizmos.color = Color.yellow;
                Gizmos.DrawSphere(_boids[i].position, 1f);
                Gizmos.DrawLine(_boids[i].position, _boids[i].position + _boids[i].velocity);
                Debug.Log($"Boid {i} has position {_boids[i].position.ToString()}");
                /*
                for (int x = _boids[i].minX; x <= _boids[i].maxX; x++) {
                    for (int y = _boids[i].minY; y <= _boids[i].maxY; y++) {
                        for (int z = _boids[i].minZ; z <= _boids[i].maxZ; z++) {
                            Gizmos.color = visualizeBoidFrequencyColor;
                            Gizmos.DrawLine(
                                visualGrid[GetProjectedHashIndex(_boids[i].hashPosition)].position, 
                                visualGrid[GetProjectedHashIndex(x,y,z)].position
                            );
                        }
                    }
                }
                */
            }
            if(visualizeBoidFrequency) {
                Color gridColor = visualizeBoidFrequencyColor;
                Vector3 boidGridSizeVector3 = new Vector3(boidGridCellSize[0], boidGridCellSize[1], boidGridCellSize[2]);
                for(int j = 0; j < visualGrid.Length; j++) {
                    GridCell cell = visualGrid[j];
                    if (cell.n <= 0) continue;
                    Debug.Log($"GIZMOS - Cell {j} has separation value of {cell.separation.ToString()}");
                    gridColor.a = Mathf.Clamp((float)cell.n / ((float)numBoids/1f), 0f, 1f);
                    Gizmos.color = gridColor;   
                    Gizmos.DrawCube(cell.position,boidGridSizeVector3);
                }
            }
        }
    }

    // Start is called before the first frame update
    private void Awake() {
        CalculateDimensions();
        //InitializeTargets();
        InitializeBoids();
        InitializeGrid();
        InitializeShaderVariables();
        InitializeBuffers();

        Debug.Log("=== UPDATINMG STARTOOOO ===");
    }

    private void CalculateDimensions() {
        dimensions = new float[3];
        dimensions[0] = Mathf.Abs(maxBounds[0] - minBounds[0]);
        dimensions[1] = Mathf.Abs(maxBounds[1] - minBounds[1]);
        dimensions[2] = Mathf.Abs(maxBounds[2] - minBounds[2]);
        boidGridCellSize = new float[3];
        boidGridCellSize[0] = dimensions[0] / (float)boidGridDiscretizations[0];
        boidGridCellSize[1] = dimensions[1] / (float)boidGridDiscretizations[1];
        boidGridCellSize[2] = dimensions[2] / (float)boidGridDiscretizations[2];
        obstacleGridSize = new float[3];
        obstacleGridSize[0] = dimensions[0] / (float)obstacleGridDiscretizations[0];
        obstacleGridSize[1] = dimensions[1] / (float)obstacleGridDiscretizations[1];
        obstacleGridSize[2] = dimensions[2] / (float)obstacleGridDiscretizations[2];
        Debug.Log($"Boid Grid Size: ({boidGridCellSize[0]}, {boidGridCellSize[1]}, {boidGridCellSize[2]})");
    }

    /*
    private void InitializeTargets() {
        if (targets.Length > 0) {
            _targets = new Target[targets.Length];
            for(int i = 0; i < targets.Length; i++) {
                Target _target = new Target();
                _target.position = targets[i].position;
                _target.numBoids = 0;
                _target.isActive = 1;
                _targets[i] = _target;
            }
        } else {
            _targets = new Target[1];
            Target _target = new Target();
            _target.position = Vector3.zero;
            _target.numBoids = 0;
            _target.isActive = 0;
            _targets[0] = _target;
        }
    }
    */

    private void InitializeBoids() {
        int totalNumBoidsAttachedToTarget = 0;
        boids = new Boid[numBoids];
        for(int i = 0; i < numBoids; i++) {
            Vector3 initPosition = new Vector3(
                Random.Range(minBounds[0],maxBounds[0]),
                Random.Range(minBounds[1],maxBounds[1]),
                Random.Range(minBounds[2],maxBounds[2])
            );
            Vector3 initVelocity = new Vector3(
                Random.Range(-1f,1f),
                Random.Range(-1f,1f),
                Random.Range(-1f,1f)
            );
            Vector3Int hashPosition = GetHashPosition(initPosition);

            Boid boid = new Boid();
            boid.px = initPosition.x;
            boid.py = initPosition.y;
            boid.pz = initPosition.z;
            boid.vx = initVelocity.x;
            boid.vy = initVelocity.y;
            boid.vz = initVelocity.z;
            boid.position = initPosition;
            boid.velocity = initVelocity;
            boid.hashPosition = hashPosition;
            boid.projectedHashPosition = GetProjectedHashIndex(boid.hashPosition);
            boid.separationInfluence = Vector3.zero;
            boid.alignmentInfluence = Vector3.zero;
            boid.cohesionInfluence = Vector3.zero;
            //boid.targetInfluence = Vector3.zero;
            boid.edgeInfluence = Vector3.zero;
            boid.targetIndex = -1;
            boid.numUpdates = 0;
            Debug.Log($"ORIGINAL HASH POSITION: {boid.hashPosition}");
            Debug.Log($"ORIGINAL PROJECTED HASH POSITION: {boid.projectedHashPosition}");
            /*
            if(_targets.Length > 0 && totalNumBoidsAttachedToTarget < _targets.Length * boidsPerTarget) {
                int targetIndex = Random.Range(0,_targets.Length);
                if (_targets[targetIndex].numBoids < boidsPerTarget) {
                    boid.targetIndex = targetIndex;
                    _targets[targetIndex].numBoids += 1;
                    totalNumBoidsAttachedToTarget += 1;
                }
            }
            */
            boids[i] = boid;
        }
    }

    private void InitializeGrid() {
        boidsGrid = new GridCell[boidGridDiscretizations[0] * boidGridDiscretizations[1] * boidGridDiscretizations[2]];
        for(int x = 0; x < boidGridDiscretizations[0]; x++) {
            for(int y = 0; y < boidGridDiscretizations[1]; y++) {
                for(int z = 0; z < boidGridDiscretizations[2]; z++) {
                    int gridIndex = GetProjectedHashIndex(x,y,z);
                    GridCell cell;
                    cell.separation = Vector3.zero;
                    cell.alignment = Vector3.zero;
                    cell.cohesion = Vector3.zero;
                    cell.position = new Vector3(
                        minBounds[0] + x*boidGridCellSize[0] + (boidGridCellSize[0]/2f),
                        minBounds[1] + y*boidGridCellSize[1] + (boidGridCellSize[1]/2f),
                        minBounds[2] + z*boidGridCellSize[2] + (boidGridCellSize[2]/2f)
                    );
                    cell.n = 0;
                    boidsGrid[gridIndex] = cell;
                }
            }
        }
        Debug.Log($"Number of grid cells: {boidsGrid.Length}");

        for(int i = 0; i < boids.Length; i++) {
            //Debug.Log(boids[i].position.ToString() + " | " +  boids[i].hashPosition.ToString() + " | " + projectedIndex.ToString());
            boidsGrid[boids[i].projectedHashPosition].separation += boids[i].position;
            boidsGrid[boids[i].projectedHashPosition].alignment += boids[i].velocity;
            boidsGrid[boids[i].projectedHashPosition].cohesion += boids[i].velocity;
            boidsGrid[boids[i].projectedHashPosition].n += 1;
        }

        for(int j = 0; j < boidsGrid.Length; j++) {
            if (boidsGrid[j].n > 0) Debug.Log($"Boid Grid Cell {j} n: {boidsGrid[j].n}");
        }
    }

    private void InitializeShaderVariables() {
        
        boidUpdateNeighborCellsKernel = boidShader.FindKernel("UpdateBoidNeighborCells");
        boidUpdateInfluenceVectorsKernel = boidShader.FindKernel("UpdateBoidInfluenceVectors");
        boidUpdateStatusKernel = boidShader.FindKernel("UpdateBoidStatus");
        //updateKernel = boidShader.FindKernel("UpdateBoids");

        boidShader.SetInt("numBoids", numBoids);
        boidShader.SetFloats("minBounds", minBounds);
        boidShader.SetFloats("maxBounds", maxBounds);
        boidShader.SetFloats("dimensions",dimensions);
        boidShader.SetFloats("boidGridCellSize", boidGridCellSize);
        boidShader.SetInts("boidGridDiscretizations",boidGridDiscretizations);

        UpdateShaderVariables();
    }

    private void UpdateShaderVariables(float deltaTime = -1f) {
        // Boid configurations
        boidShader.SetInt("visualRange", visualRange);
        boidShader.SetInt("closeRange", closeRange);
        boidShader.SetFloat("minSpeed", minSpeed);
        boidShader.SetFloat("maxSpeed", maxSpeed);

        // Weight configurations
        boidShader.SetFloat("cohesionFactor", cohesionFactor);
        boidShader.SetFloat("alignmentFactor", alignmentFactor);
        boidShader.SetFloat("separationFactor", separationFactor);
        // boidShader.SetFloat("obstacleFactor", obstacleFactor);
        boidShader.SetFloat("turnFactor", turnFactor);
        boidShader.SetFloat("targetFactor", targetFactor);

        // World Configurations
        boidShader.SetFloats("innerDim_min", innerMinBounds);
        boidShader.SetFloats("innerDim_max", innerMaxBounds);
        if (deltaTime != -1f) 
            boidShader.SetFloat("deltaTime", deltaTime);
    }

    private void InitializeBuffers() {
        boidsBuffer = new ComputeBuffer(
            numBoids, 
            sizeof(float)*24 + sizeof(int)*18
        );
        boidsBuffer.SetData(boids);

        boidsGridBuffer = new ComputeBuffer(
            boidGridDiscretizations[0] * boidGridDiscretizations[1] * boidGridDiscretizations[2], 
            sizeof(float)*12 + sizeof(int)
        );
        boidsGridBuffer.SetData(boidsGrid);
        
        /*
        targetsBuffer = new ComputeBuffer(
            _targets.Length, 
            sizeof(float)*3+sizeof(int)*2
        );
        targetsBuffer.SetData(_targets);
        */

        boidShader.SetBuffer(boidUpdateNeighborCellsKernel, "boids", boidsBuffer);

        boidShader.SetBuffer(boidUpdateInfluenceVectorsKernel, "boids", boidsBuffer);
        boidShader.SetBuffer(boidUpdateInfluenceVectorsKernel, "grid", boidsGridBuffer);
        //boidShader.SetBuffer(boidUpdateInfluenceVectorsKernel, "targets", targetsBuffer);

        boidShader.SetBuffer(boidUpdateStatusKernel, "boids", boidsBuffer);
        boidShader.SetBuffer(boidUpdateStatusKernel, "grid", boidsGridBuffer);
        //boidShader.SetBuffer(boidUpdateStatusKernel, "targets", targetsBuffer);
        
        /*
        boidShader.SetBuffer(updateKernel, "boids", boidsBuffer);
        boidShader.SetBuffer(updateKernel, "grid", boidsGridBuffer);
        boidShader.SetBuffer(updateKernel, "targets", targetsBuffer);
        */
    }

    /*
    private void UpdateTargets() {
        for(int i = 0; i < _targets.Length; i++) {
            if(_targets[i].isActive == 1)
                _targets[i].position = targets[i].position;
        }
    }
    private void UpdateBuffers() {
        targetsBuffer.SetData(_targets);
    }
    */

    // Update is called once per frame
    private void Update() {
            
            float deltaTime = Time.deltaTime;
            //UpdateTargets();
            UpdateShaderVariables(deltaTime);
            //UpdateBuffers();

            boidShader.Dispatch(boidUpdateNeighborCellsKernel, 200,1,1);
            boidShader.Dispatch(boidUpdateInfluenceVectorsKernel, 200,1,1);
            boidShader.Dispatch(boidUpdateStatusKernel, 200,1,1);
            //boidShader.Dispatch(updateKernel,200,1,1);

            /*
            GridCell[] visualGrid = new GridCell[boidsGridBuffer.count];
            boidsGridBuffer.GetData(visualGrid);    
            for(int i = 0; i < visualGrid.Length; i++) {
                GridCell cell = visualGrid[i];
                if (cell.n != 0) {
                    Debug.Log($"{i}: {cell.n}");
                }
            }
            */

            boidsBuffer.GetData(boids);
            for(int i = 0; i < boids.Length; i++) {
                //Debug.Log($"Boid {i} # Updates: {boids[i].numUpdates}");
                //Debug.Log($"Boid {i} Old hash position: {boids[i].oldHashPosition}");
                //Debug.Log($"Boid {i} New hash position: {boids[i].hashPosition}");
                //Debug.Log($"Boid {i} Old projected hash position: {boids[i].oldProjectedHashPosition}");
                //Debug.Log($"Boid {i} New projected hash position: {boids[i].projectedHashPosition}");

            }

    }

    private Vector3Int GetHashPosition(Vector3 position) {
        /*
        int hashX = Mathf.FloorToInt(((position.x - minBounds[0]) / dimensions[0]) * boidGridDiscretizations[0]);
        int hashY = Mathf.FloorToInt(((position.y - minBounds[1]) / dimensions[1]) * boidGridDiscretizations[1]);
        int hashZ = Mathf.FloorToInt(((position.z - minBounds[2]) / dimensions[2]) * boidGridDiscretizations[2]);
        */
        int hashX = Mathf.FloorToInt((position.x - minBounds[0])/boidGridCellSize[0]);
        int hashY = Mathf.FloorToInt((position.y - minBounds[1])/boidGridCellSize[1]);
        int hashZ = Mathf.FloorToInt((position.z - minBounds[2])/boidGridCellSize[2]);
        return new Vector3Int(hashX, hashY, hashZ);
    }

    private int GetProjectedHashIndex(int x, int y, int z) {
        return x + (y * boidGridDiscretizations[0]) + (z * boidGridDiscretizations[0] * boidGridDiscretizations[1]);
    }
    private int GetProjectedHashIndex(Vector3Int p) {
        return p.x + (p.y * boidGridDiscretizations[0]) + (p.z * boidGridDiscretizations[0] * boidGridDiscretizations[1]);
    }

    private void OnDestroy()
    {
        ReleaseBuffers();
    }

    private void ReleaseBuffers()
    {
        boidsBuffer.Dispose();
        boidsGridBuffer.Dispose();
        //targetsBuffer.Dispose();
    }
}
