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
        public Vector3 velocity;

        public int minX, maxX, minY, maxY, minZ, maxZ;
        public int numNeighbors, numCloseNeighbors;

        public int targetIndex;
    };

    struct GridCell {
        public Vector3 separation;
        public Vector3 alignment;
        public Vector3 cohesion;
        public int n;
    };

    // Private list of structs that represent targets. Sent into a custom targets buffer
    private struct Target {
        public Vector3 position;
        public int numBoids;
    }
    
    [Header("= CORE SETTINGS =")]
    [Tooltip("Total # of boids in the simulation")] 
        public int numBoids = 32;
    [Tooltip("The minimum limits of the world space boundaries of the simulation. Stored in XYZ format")] 
        public float[] minBounds = { -100f,-100f,-100f };
    [Tooltip("The maximum limits of the world space boundaries of the simulation. Stored in XYZ format")] 
        public float[] maxBounds = { 100f,100f,100f };
    // The world space size of the simulation
    private float[] dimensions;
    [Tooltip("The # of cells along the X, Y, and Z axes each")]
        public int[] numCells = { 10,10,10 };
    [Tooltip("The Compute Shader for Boid Behavior")]
        public ComputeShader boidShader;

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

    [Header("= WORLD CONFIGURATIONS =")]
    [Tooltip("The inner XYZ grid dimension where boids are restricted to stay in")]
        public float[] innerHashPos_min = { 2f, 2f, 2f };
    [Tooltip("The inner XYZ grid dimension where boids are restricted to stay in")]
        public float[] innerHashPos_max = { 8f, 8f, 8f };

    [Header("= TARGETS AND OBSTACLES =")]
    [Tooltip("All gameobjects that are considered targets")]
        public Transform[] targets = new Transform[0];
    [Tooltip("The total number of boids that a target can be associated with")]
        public int boidsPerTarget = 200;

    private Target[] _targets;
    private Boid[] boids;
    private GridCell[] boidsGrid;

    private int updateKernel;
    private ComputeBuffer boidsBuffer;
    private ComputeBuffer boidsGridBuffer;
    private ComputeBuffer targetsBuffer;

    [SerializeField] private float updateTimeDelay = 0.1f;

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
        Vector3 minInnerBoundVector = new Vector3(innerHashPos_min[0], innerHashPos_min[1], innerHashPos_min[2]);
        Vector3 d2 = new Vector3(
            innerHashPos_max[0]-innerHashPos_min[0], 
            innerHashPos_max[1]-innerHashPos_min[1], 
            innerHashPos_max[2]-innerHashPos_min[2]
        );
        Gizmos.DrawWireCube(minInnerBoundVector + d2/2f, d2);

        if (Application.isPlaying) {
            Gizmos.color = Color.yellow;
            int bCount = boidsBuffer.count;
            Boid[] _boids = new Boid[bCount];
            boidsBuffer.GetData(_boids);
            for(int i = 0; i < bCount; i++) {
                Gizmos.DrawSphere(_boids[i].position, 1f);
            }
        }
    }

    // Start is called before the first frame update
    private void Start() {
        CalculateDimensions();
        InitializeTargets();
        InitializeBoids();
        InitializeGrid();
        InitializeShaderVariables();
        InitializeBuffers();

        StartCoroutine(CustomUpdate());
    }

    private void CalculateDimensions() {
        dimensions = new float[3];
        dimensions[0] = Mathf.Abs(maxBounds[0] - minBounds[0]);
        dimensions[1] = Mathf.Abs(maxBounds[1] - minBounds[1]);
        dimensions[2] = Mathf.Abs(maxBounds[2] - minBounds[2]);
    }

    private void InitializeTargets() {
        _targets = new Target[targets.Length];
        for(int i = 0; i < targets.Length; i++) {
            Target _target = new Target();
            _target.position = targets[i].position;
            _target.numBoids = 0;
            _targets[i] = _target;
        }
    }

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
            boid.targetIndex = -1;
            if(_targets.Length > 0 && totalNumBoidsAttachedToTarget < _targets.Length * boidsPerTarget) {
                int targetIndex = Random.Range(0,_targets.Length);
                if (_targets[targetIndex].numBoids < boidsPerTarget) {
                    boid.targetIndex = targetIndex;
                    _targets[targetIndex].numBoids += 1;
                    totalNumBoidsAttachedToTarget += 1;
                }
            }
            boids[i] = boid;
        }
    }

    private void InitializeGrid() {
        boidsGrid = new GridCell[numCells[0] * numCells[1] * numCells[2]];
        for(int x = 0; x < numCells[0]; x++) {
            for(int y = 0; y < numCells[1]; y++) {
                for(int z = 0; z < numCells[2]; z++) {
                    int gridIndex = GetProjectedHashIndex(x,y,z);
                    GridCell cell;
                    cell.separation = Vector3.zero;
                    cell.alignment = Vector3.zero;
                    cell.cohesion = Vector3.zero;
                    cell.n = 0;
                    boidsGrid[gridIndex] = cell;
                }
            }
        }
        // Debug.Log("Number of grid cells: " + boidsGrid.Length.ToString());

        for(int i = 0; i < boids.Length; i++) {
            int projectedIndex = GetProjectedHashIndex(boids[i].hashPosition);
            //Debug.Log(boids[i].position.ToString() + " | " +  boids[i].hashPosition.ToString() + " | " + projectedIndex.ToString());
            boidsGrid[projectedIndex].separation += boids[i].position;
            boidsGrid[projectedIndex].alignment += boids[i].velocity;
            boidsGrid[projectedIndex].cohesion += boids[i].velocity;
            boidsGrid[projectedIndex].n += 1;
        }
    }

    private void InitializeShaderVariables() {

        updateKernel = boidShader.FindKernel("UpdateBoids");
        boidShader.SetInt("numBoids", numBoids);
        boidShader.SetFloats("minBounds", minBounds);
        boidShader.SetFloats("maxBounds", maxBounds);
        boidShader.SetFloats("dimensions",dimensions);
        boidShader.SetInts("numCells",numCells);

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
        boidShader.SetFloats("innerDim_min", innerHashPos_min);
        boidShader.SetFloats("innerDim_max", innerHashPos_max);
        if (deltaTime != -1f) 
            boidShader.SetFloat("deltaTime", deltaTime);
    }

    private void InitializeBuffers() {
        boidsBuffer = new ComputeBuffer(numBoids, sizeof(float)*12 + sizeof(int)*12);
        boidsBuffer.SetData(boids);

        boidsGridBuffer = new ComputeBuffer(numCells[0] * numCells[1] * numCells[2], sizeof(float)*9 + sizeof(int));
        boidsGridBuffer.SetData(boidsGrid);

        targetsBuffer = new ComputeBuffer(targets.Length, sizeof(float)*3+sizeof(int));
        targetsBuffer.SetData(_targets);

        boidShader.SetBuffer(updateKernel, "boids", boidsBuffer);
        boidShader.SetBuffer(updateKernel, "grid", boidsGridBuffer);
        boidShader.SetBuffer(updateKernel, "targets", targetsBuffer);
    }

    private void UpdateTargets() {
        for(int i = 0; i < _targets.Length; i++) {
            _targets[i].position = targets[i].position;
        }
    }
    private void UpdateBuffers() {
        targetsBuffer.SetData(_targets);
    }

    // Update is called once per frame
    private IEnumerator CustomUpdate() {
        while(true) {
            float deltaTime = Time.deltaTime;
            UpdateTargets();
            UpdateShaderVariables(deltaTime);
            UpdateBuffers();

            boidShader.Dispatch(updateKernel,100,1,1);

            /*
            boidsBuffer.GetData(boids);
            for(int i = 0; i < boids.Length; i++) {
                Boid boid = boids[i];
                Debug.Log(boid.numNeighbors.ToString() + "|" + boid.numCloseNeighbors.ToString());
            }
            */

            yield return new WaitForSeconds(updateTimeDelay);
        }
    }

    private Vector3Int GetHashPosition(Vector3 position) {
        int hashX = Mathf.FloorToInt(((position.x - minBounds[0]) / dimensions[0]) * numCells[0]);
        int hashY = Mathf.FloorToInt(((position.y - minBounds[1]) / dimensions[1]) * numCells[1]);
        int hashZ = Mathf.FloorToInt(((position.z - minBounds[2]) / dimensions[2]) * numCells[2]);
        return new Vector3Int(hashX, hashY, hashZ);
    }

    private int GetProjectedHashIndex(int x, int y, int z) {
        return x + (y * numCells[0]) + (z * numCells[0] * numCells[1]);
    }
    private int GetProjectedHashIndex(Vector3Int p) {
        return p.x + (p.y * numCells[0]) + (p.z * numCells[0] * numCells[1]);
    }

    private void OnDestroy()
    {
        ReleaseBuffers();
    }

    private void ReleaseBuffers()
    {
        boidsBuffer.Dispose();
        boidsGridBuffer.Dispose();
        targetsBuffer.Dispose();
    }
}
