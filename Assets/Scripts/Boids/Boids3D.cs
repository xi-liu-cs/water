using Unity.Mathematics;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Random = UnityEngine.Random;

public class Boids3D : Grid3D
{
    public static Boids3D current;

    struct Boid3D {
        public float3 position;
        public int3 gridIndices;
        public int projectedGridIndex;
    }

    [Header("== REFERENCES ==")]
    public ComputeShader boidsShader;
    public bool awaitInitialization = false;

    [Header("== BOID SETTINGS ==")]
    public int numBoids = 32;
    public float boidSize = 0.05f;
    public Color boidColor = Color.red;
    public float visualRange = 0.5f;
    public float innerRange = 0.15f;
    public float maxSpeed = 1.5f;
    public float minSpeed = 0.5f;
    public float cohesionFactor = 1;
    public float separationFactor = 30;
    public float alignmentFactor = 5;
    public float sphFactor = 1f;
    public float turnSpeed = 2f;
    public float dt = -1f;

    private Vector3[] velocities;
    private int[] offsets;
    private Boid3D[] boids;
    private int[] rearrangedBoids; 
    // To get the index for a grid cell, pass the grid cell's XYZ indices into `Boid3DHelpers.GetProjectedGridIndexFromXYZ()`
    private int[] grid;

    [Header("== GPU SETTINGS ==")]
    public bool useGPU = false;
    private int _CPU_LIMIT = 2048;
    private int _BLOCK_SIZE = 256;
    private int numBlocks;

    private int generateGridKernel;
    private int generateBoidsKernel;
    private int clearGridKernel;
    private int updateGridCellCountsKernel;
    private int prefixSumKernel;
    private int sumBlocksKernel;
    private int addSumsKernel;
    private int rearrangeBoidsKernel;
    private int updateBoidsKernel;

    public ComputeBuffer boidsBuffer;
    public ComputeBuffer boidVelocitiesBuffer;
    
    private ComputeBuffer boidOffsetsBuffer;
    private ComputeBuffer rearrangedBoidsBuffer;
    private ComputeBuffer gridBuffer;
    private ComputeBuffer gridOffsetBuffer;
    private ComputeBuffer gridSumsBuffer1, gridSumsBuffer2;

    [Header("== DEBUGGING CONFIGURATIONS ==")]
    public bool verbose = false;
    public bool visualizeBoidGridCells = false;
    public bool visualizeBoidVelocities = false;
    private Vector3 debugPos;
    
    //public MeshFilter[] environmentMeshes = new MeshFilter[0];

    void OnDrawGizmos() {
        if (showGridAxes) {
            Gizmos.color = gridAxesColor;
            Vector3 p1, p2, p3;
            for(int z = 0; z <= numGridCellsPerAxis.z; z++) {
                p1 = Grid3DHelpers.GetGridCellWorldPositionFromXYZIndices(numGridCellsPerAxis, origin, gridCellSize, 0, 0, z, false);
                p2 = Grid3DHelpers.GetGridCellWorldPositionFromXYZIndices(numGridCellsPerAxis, origin, gridCellSize, numGridCellsPerAxis.x, 0, z, false);
                p3 = Grid3DHelpers.GetGridCellWorldPositionFromXYZIndices(numGridCellsPerAxis, origin, gridCellSize, 0, numGridCellsPerAxis.y, z, false);
                Gizmos.DrawLine(p1,p2);
                Gizmos.DrawLine(p1,p3);
            }
            for(int x = 0; x <= numGridCellsPerAxis.x; x++) {
                p1 = Grid3DHelpers.GetGridCellWorldPositionFromXYZIndices(numGridCellsPerAxis, origin, gridCellSize, x, 0, 0, false);
                p2 = Grid3DHelpers.GetGridCellWorldPositionFromXYZIndices(numGridCellsPerAxis, origin, gridCellSize, x, 0, numGridCellsPerAxis.z, false);
                Gizmos.DrawLine(p1,p2);
            }
            for(int y = 0; y <= numGridCellsPerAxis.y; y++) {
                p1 = Grid3DHelpers.GetGridCellWorldPositionFromXYZIndices(numGridCellsPerAxis, origin, gridCellSize, 0, y, 0, false);
                p2 = Grid3DHelpers.GetGridCellWorldPositionFromXYZIndices(numGridCellsPerAxis, origin, gridCellSize, 0, y, numGridCellsPerAxis.z, false );
                Gizmos.DrawLine(p1,p2);
            }
        }
        if (showBounds) {
            Gizmos.color = gridBoundsColor;
            Gizmos.DrawWireCube(origin, bounds);
        }
        if (showInnerBounds) {
            Gizmos.color = gridInnerBoundsColor;
            Gizmos.DrawWireCube(origin, innerBounds);
        }

        if (Application.isPlaying) {
            base.DrawGridAxes();
            base.DrawGridInnerBounds();
            base.DrawGridBounds();
            if (useGPU) {
                Boid3D[] tempBoids = new Boid3D[numBoids];
                boidsBuffer.GetData(tempBoids);
                int[] tempGrid = new int[numGridCells];
                gridBuffer.GetData(tempGrid);
                float3[] tempVelocities = new float3[numBoids];
                boidVelocitiesBuffer.GetData(tempVelocities);
                for(int i = 0; i < numBoids; i++) {
                    Boid3D boid = tempBoids[i];
                    Gizmos.color = boidColor;
                    Vector3 pos = new Vector3(boid.position[0], boid.position[1], boid.position[2]);
                    Gizmos.DrawSphere(pos, boidSize);
                    if (visualizeBoidGridCells) {
                        Gizmos.color = new Vector4(0f,0f,1f,Mathf.Clamp((float)tempGrid[boid.projectedGridIndex]/2f,0f,1f));
                        Vector3 debugPosGridPos = Grid3DHelpers.GetGridCellWorldPositionFromGivenPosition(numGridCellsPerAxis, origin, gridCellSize, pos);
                        Gizmos.DrawCube(debugPosGridPos, new Vector3(gridCellSize, gridCellSize, gridCellSize));
                    }
                    if (visualizeBoidVelocities) {
                        Gizmos.color = Color.yellow;
                        Gizmos.DrawRay(pos, tempVelocities[i]);
                    }     
                }
            } else {
                for(int i = 0; i < numBoids; i++) {
                    Boid3D boid = boids[i];
                    Gizmos.color = boidColor;
                    Vector3 pos = new Vector3(boid.position[0], boid.position[1], boid.position[2]);
                    Gizmos.DrawSphere(pos, boidSize);
                    if (visualizeBoidGridCells) {
                        Gizmos.color = new Vector4(0f,0f,1f,0.5f);
                        Vector3 debugPosGridPos = Grid3DHelpers.GetGridCellWorldPositionFromGivenPosition(numGridCellsPerAxis, origin, gridCellSize, pos);
                        Gizmos.DrawCube(debugPosGridPos, new Vector3(gridCellSize, gridCellSize, gridCellSize));
                    }        
                }
            }
        }
    }

    private void Awake() {
        current = this;
        if (!awaitInitialization) Initialize();
    }

    public void Initialize() {
        visualRange = gridCellSize;
        
        // Initialize `numBlocks` for GPU's sake
        numBlocks = Mathf.CeilToInt((float)numGridCells / (float)_BLOCK_SIZE);
        if (verbose) Debug.Log($"Number of blocks: {numBlocks}");

        if (!useGPU && numBoids > _CPU_LIMIT) {
            Debug.Log($"WARNING - `numBoids` over the CPU limit. Will reduce down to clamp to {_CPU_LIMIT} boids");
            numBoids = _CPU_LIMIT;
        }

        //ProcessEnvironmentMeshes();

        if (useGPU) {
            RetrieveShaderKernels();
            InitializeShaderVariables();
            InitializeShaderBuffers();
            boidsShader.Dispatch(generateGridKernel, numBlocks,1,1);
            boidsShader.Dispatch(generateBoidsKernel, Mathf.CeilToInt((float)numBoids / (float)_BLOCK_SIZE), 1, 1);
        } else {
            GenerateGrid();
            GenerateBoids();
        }
    }

    /*
    private void ProcessEnvironmentMeshes() {
        if (environmentMeshes.Length == 0) return;
        foreach(MeshFilter meshFilter in environmentMeshes) {
            
            for(int i = 0; i < meshFilter.mesh.triangles.Length-2; i+=3) {
                Debug.Log($"{meshFilter.mesh.vertices[meshFilter.mesh.triangles[i]]}, {meshFilter.mesh.vertices[meshFilter.mesh.triangles[i+1]]}, {meshFilter.mesh.vertices[meshFilter.mesh.triangles[i+2]]}");
            }
        }
    }
    */

    private void GenerateGrid() {
        grid = new int[numGridCells];
        for(int x = 0; x < numGridCellsPerAxis.x; x++) {
            for(int y = 0; y < numGridCellsPerAxis.y; y++) {
                for(int z = 0; z < numGridCellsPerAxis.z; z++) {
                    int projectedIndex = Boid3DHelpers.GetProjectedGridIndexFromXYZ(numGridCellsPerAxis, x, y, z);
                    grid[projectedIndex] = 0;
                }
            }
        }
    }

    private void GenerateBoids() {
        boids = new Boid3D[numBoids];
        velocities = new Vector3[numBoids];
        rearrangedBoids = new int[numBoids];
        for(int i = 0; i < numBoids; i++) {
            Boid3D boid = new Boid3D();
            Vector3 pos =  new Vector3(
                Random.Range(-bounds.x/2f, bounds.x/2f),
                Random.Range(-bounds.y/2f, bounds.y/2f),
                Random.Range(-bounds.z/2f, bounds.z/2f)
            );
            Vector3 vel = new Vector3(
                Random.Range(-1f,1f),
                Random.Range(-1f,1f),
                Random.Range(-1f,1f)
            );
            Vector3 normalizedVel = vel /= vel.magnitude;
            boid.position = new(pos.x, pos.y, pos.z);
            velocities[i] = normalizedVel;
            Vector3Int xyzIndices = Grid3DHelpers.GetGridXYZIndices(numGridCellsPerAxis, origin, gridCellSize, pos);
            boid.gridIndices = new(xyzIndices.x, xyzIndices.y, xyzIndices.z);
            boid.projectedGridIndex = Grid3DHelpers.GetProjectedGridIndexFromXYZ(numGridCellsPerAxis, xyzIndices);
            boids[i] = boid;
        }
    }

    private void RetrieveShaderKernels() {
        generateGridKernel = boidsShader.FindKernel("GenerateGrid");
        generateBoidsKernel = boidsShader.FindKernel("GenerateBoids");
        clearGridKernel = boidsShader.FindKernel("ClearGrid");
        updateGridCellCountsKernel = boidsShader.FindKernel("UpdateGridCellCounts");
        prefixSumKernel = boidsShader.FindKernel("PrefixSum");
        sumBlocksKernel = boidsShader.FindKernel("SumBlocks");
        addSumsKernel = boidsShader.FindKernel("AddSums");
        rearrangeBoidsKernel = boidsShader.FindKernel("RearrangeBoids");
        updateBoidsKernel = boidsShader.FindKernel("UpdateBoids");
    }

    private void InitializeShaderVariables() {
        boidsShader.SetFloats("origin", originF);

        boidsShader.SetFloat("boundsX", innerBounds.x);
        boidsShader.SetFloat("boundsY", innerBounds.y);
        boidsShader.SetFloat("boundsZ", innerBounds.z);

        boidsShader.SetFloat("gridCellSize", gridCellSize);
        boidsShader.SetInt("numGridCells", numGridCells);
        boidsShader.SetFloat("dimensionsX", numGridCellsPerAxis.x);
        boidsShader.SetFloat("dimensionsY", numGridCellsPerAxis.y);
        boidsShader.SetFloat("dimensionsZ", numGridCellsPerAxis.z);

        boidsShader.SetInt("numBoids", numBoids);
        boidsShader.SetFloat("visualRange", gridCellSize);
        boidsShader.SetFloat("innerRange", innerRange);
        boidsShader.SetInt("randomSeed", Random.Range(0,int.MaxValue));

        boidsShader.SetInt("numBlocks", numBlocks);

        UpdateShaderVariables();
    }

    private void InitializeShaderBuffers() {
        boidsBuffer = new ComputeBuffer(numBoids, sizeof(float)*3 + sizeof(int)*4);
        boidVelocitiesBuffer = new ComputeBuffer(numBoids, sizeof(float)*3);
        boidOffsetsBuffer = new ComputeBuffer(numBoids, sizeof(int));
        rearrangedBoidsBuffer = new ComputeBuffer(numBoids, sizeof(int));

        gridBuffer = new ComputeBuffer(numGridCells, sizeof(int));
        //gridBuffer.SetData(grid);
        gridOffsetBuffer = new ComputeBuffer(numGridCells, sizeof(int));
        //int[] gridSums = new int[numBlocks];
        gridSumsBuffer1 = new ComputeBuffer(numBlocks, sizeof(int));
        //gridSumsBuffer1.SetData(gridSums);
        gridSumsBuffer2 = new ComputeBuffer(numBlocks, sizeof(int));
        //gridSumsBuffer2.SetData(gridSums);

        boidsShader.SetBuffer(generateGridKernel, "grid", gridBuffer);
        boidsShader.SetBuffer(generateBoidsKernel, "boids", boidsBuffer);
        boidsShader.SetBuffer(generateBoidsKernel, "velocities", boidVelocitiesBuffer);

        boidsShader.SetBuffer(clearGridKernel, "grid", gridBuffer);

        boidsShader.SetBuffer(updateGridCellCountsKernel, "boids", boidsBuffer);
        boidsShader.SetBuffer(updateGridCellCountsKernel, "grid", gridBuffer);
        boidsShader.SetBuffer(updateGridCellCountsKernel, "offsets", gridOffsetBuffer);

        boidsShader.SetBuffer(prefixSumKernel, "gridOffsetBufferIn", gridBuffer);
        boidsShader.SetBuffer(prefixSumKernel, "gridOffsetBuffer", gridOffsetBuffer);
        boidsShader.SetBuffer(prefixSumKernel, "gridSumsBuffer", gridSumsBuffer2);

        boidsShader.SetBuffer(addSumsKernel, "gridOffsetBuffer", gridOffsetBuffer);

        boidsShader.SetBuffer(rearrangeBoidsKernel, "boids", boidsBuffer);
        boidsShader.SetBuffer(rearrangeBoidsKernel, "rearrangedBoids", rearrangedBoidsBuffer);
        boidsShader.SetBuffer(rearrangeBoidsKernel, "gridOffsetBuffer", gridOffsetBuffer);
        boidsShader.SetBuffer(rearrangeBoidsKernel, "offsets", boidOffsetsBuffer);

        boidsShader.SetBuffer(updateBoidsKernel, "rearrangedBoids", rearrangedBoidsBuffer);
        boidsShader.SetBuffer(updateBoidsKernel, "boids", boidsBuffer);
        boidsShader.SetBuffer(updateBoidsKernel, "velocities", boidVelocitiesBuffer);
        boidsShader.SetBuffer(updateBoidsKernel, "gridOffsetBuffer", gridOffsetBuffer);
    }

    public void AddParticleBuffer() {
        boidsShader.SetBuffer(updateBoidsKernel, "particles", ParticleManager.current.particle_buffer);
        boidsShader.SetBuffer(updateBoidsKernel, "particleVelocities", ParticleManager.current.velocity_buffer);
    }

    // Update is called once per frame
    void Update() {
        if (useGPU) GPU_Process();
        else CPU_Process();
    }

    
    private void GPU_Process() {
        UpdateShaderVariables();
        // Clean up `grid` so that all values = 0
        boidsShader.Dispatch(clearGridKernel, numBlocks,1,1);
        // We should update each boid with their current indice and projected index, while also declaring offset and performing atomic addition
        boidsShader.Dispatch(updateGridCellCountsKernel, Mathf.CeilToInt((float)numBoids / (float)_BLOCK_SIZE),1,1);
        if (verbose) {
            DebugBufferBoid3D("Boid Initial Positions", Mathf.Min(100,numBoids), boidsBuffer);
            DebugBufferInt("Grid Values", Mathf.Min(500,numGridCells), gridBuffer);
        }

        boidsShader.Dispatch(prefixSumKernel, numBlocks, 1, 1);
        bool swap = false;
        for (int d = 1; d < numBlocks; d *= 2) {
            boidsShader.SetBuffer(sumBlocksKernel, "gridSumsBufferIn", swap ? gridSumsBuffer1 : gridSumsBuffer2);
            boidsShader.SetBuffer(sumBlocksKernel, "gridSumsBuffer", swap ? gridSumsBuffer2 : gridSumsBuffer1);
            boidsShader.SetInt("d", d);
            boidsShader.Dispatch(sumBlocksKernel, Mathf.CeilToInt((float)numBlocks / (float)_BLOCK_SIZE), 1, 1);
            swap = !swap;
        }
        boidsShader.SetBuffer(addSumsKernel, "gridSumsBufferIn", swap ? gridSumsBuffer1 : gridSumsBuffer2);
        boidsShader.Dispatch(addSumsKernel, numBlocks, 1, 1);
        if (verbose) {
            DebugBufferInt("Sums Buffer", Mathf.Min(500,numGridCells), gridOffsetBuffer);
        }

        boidsShader.Dispatch(rearrangeBoidsKernel, Mathf.CeilToInt((float)numBoids / (float)_BLOCK_SIZE), 1, 1);
        if (verbose) {
            DebugBufferInt("Rearranged Buffer", Mathf.Min(100,numBoids), rearrangedBoidsBuffer);
        }

        boidsShader.Dispatch(updateBoidsKernel, Mathf.CeilToInt((float)numBoids / (float)_BLOCK_SIZE), 1, 1);

        /*        
        if (verbose) {
            Boid3D[] tempBoids = new Boid3D[numBoids];
            boidsBuffer.GetData(tempBoids);
            int[] tempRearranged = new int[numBoids];
            rearrangedBoidsBuffer.GetData(tempRearranged);
            string top = "";
            string bottom = "";
            for(int j = 0; j < numBoids; j++) {
                top += $"I{j}\t|";
                bottom += $"B{tempBoids[tempRearranged[j]]}\t|";
            }
            Debug.Log("REARRANGED BOIDS:\n"+top+"\n"+bottom);
            
            int[] tempGrid = new int[numGridCells];
            gridOffsetBuffer.GetData(tempGrid);
            top = "";
            bottom="";
            for(int k = 0; k < numGridCells; k++) {
                top += $"{k}\t|";
                bottom += $"{tempGrid[k]}\t|";
            }
            Debug.Log("PREFIX SUMS:\n"+top+"\n"+bottom);
        }
        */
    }

    private void UpdateShaderVariables() {
        float deltaTime = (dt < 0) ? Time.deltaTime : dt;
        boidsShader.SetFloat("maxSpeed", maxSpeed);
        boidsShader.SetFloat("minSpeed", minSpeed);
        boidsShader.SetFloat("cohesionFactor", cohesionFactor);
        boidsShader.SetFloat("separationFactor", separationFactor);
        boidsShader.SetFloat("alignmentFactor", alignmentFactor);
        boidsShader.SetFloat("sphFactor", sphFactor);
        boidsShader.SetFloat("turnSpeed", turnSpeed);
        boidsShader.SetFloat("deltaTime", deltaTime);
    }

    private void CPU_Process() {
        ClearGrid();
        UpdateGridCellCounts();
        GenerateGridOffsets();
        RearrangeBoids();
        if (verbose) {
            for (int i = 0; i < numBoids; i++ ) {
                Debug.Log($"REARRANGED BOID {i}: Position={boids[rearrangedBoids[i]].position} | Grid XYZ={boids[rearrangedBoids[i]].gridIndices} | Grid Projected Index={boids[rearrangedBoids[i]].projectedGridIndex}");
            }
        }
        for(int i = 0; i < numBoids; i++) {
            int boidId = rearrangedBoids[i];
            Boid3D boid = boids[boidId];
            if (verbose) Debug.Log($"LOOKING AT REARRANGED BOID {i}");
            MergedBehaviours(ref boid, i);
            LimitSpeed(ref boid, i);
            KeepInBounds(ref boid, i);
            Vector3 currentPosition = new Vector3(boid.position[0], boid.position[1], boid.position[2]);
            currentPosition += velocities[i] * Time.deltaTime;
            boid.position = new(currentPosition.x, currentPosition.y, currentPosition.z);
            boids[i] = boid;
        }
    }

    private void ClearGrid() {
        for(int i = 0; i < grid.Length; i++) {
            grid[i] = 0;
        }
    }

    private void UpdateGridCellCounts() {
        for(int i = 0; i < numBoids; i++) {
            Vector3Int xyz = Grid3DHelpers.GetGridXYZIndices(numGridCellsPerAxis, origin, gridCellSize, boids[i].position);
            boids[i].gridIndices = new(xyz.x, xyz.y, xyz.z);
            boids[i].projectedGridIndex = Grid3DHelpers.GetProjectedGridIndexFromXYZ(numGridCellsPerAxis, xyz);
            offsets[i] = grid[boids[i].projectedGridIndex];
            grid[boids[i].projectedGridIndex] += 1;
        }
        if(verbose) {
            string top = "";
            string bottom = "";
            for(int j = 0; j < grid.Length; j++) {
                top += $"{j}\t|";
                bottom += $"{grid[j]}\t|";
            }
            Debug.Log("ORIGINAL OFFSETS:\n"+top+"\n"+bottom);
        }
    }

    private void GenerateGridOffsets() {
        for (int i = 1; i < grid.Length; i++) {
            grid[i] += grid[i - 1];
        }
        if (verbose) {
            string top = "";
            string bottom ="";
            for(int j = 0; j < grid.Length; j++) {
                top += $"{j}\t|";
                bottom += $"{grid[j]}\t|";
            }
            Debug.Log("PREFIX SUMS:\n"+top+"\n"+bottom);
        }
    }

    private void RearrangeBoids() {
        int[] temp = new int[numBoids];
        for (int i = 0; i < numBoids; i++) {
            int index = grid[boids[i].projectedGridIndex] - 1 - offsets[i];
            rearrangedBoids[index] = i;
            temp[index] = i;
        }
        if (verbose) {
            string top = "";
            string bottom = "";
            for(int j = 0; j < numBoids; j++) {
                top += $"I{j}\t|";
                bottom += $"B{temp[j]}\t|";
            }
            Debug.Log("REARRANGED BOIDS:\n"+top+"\n"+bottom);
        }
    }

    private void MergedBehaviours(ref Boid3D boid, int i) {
        float deltaTime = (dt < 0) ? Time.deltaTime : dt;
        Vector3 center = Vector3.zero;
        Vector3 close = Vector3.zero;
        Vector3 avgVel = Vector3.zero;
        Vector3 currentPosition = new Vector3(boid.position[0], boid.position[1], boid.position[2]);
        Vector3 currentVelocity = velocities[i];
        int neighbors = 0;

        // Remember: to iterate through the Z axis, we merely jump in intervals of X*Y.
        int zStep = numGridCellsPerAxis.x * numGridCellsPerAxis.y;
        // We first start at the current cell index
        // We iterate through the cells around the current cell by first iterating through the z-axis. min = current cell - X*Y. Max = cell + X*Y
        if (verbose) Debug.Log($"\tProjected Grid Index: {boid.projectedGridIndex}");
        for(int zIndex = boid.projectedGridIndex - zStep; zIndex <= boid.projectedGridIndex + zStep; zIndex += zStep) {
            if (verbose) Debug.Log($"\tZ index: {zIndex}");
            // Remember: to iterate through the Y axis, we simply add in steps of X. min = zIndex - X.
            for(int yIndex = zIndex - numGridCellsPerAxis.x; yIndex <= zIndex + numGridCellsPerAxis.x; yIndex += numGridCellsPerAxis.x) {
                // To iterate through the X, we can add and subtract by 1.
                int start = grid[yIndex-1];
                int end = grid[yIndex+1];
                if (verbose) {
                    Debug.Log($"\t\tY index: {yIndex}");
                    Debug.Log($"\t\t\tStart ({yIndex-1}): {grid[yIndex-1]} | End ({yIndex+1}): {grid[yIndex+1]}");
                }
                // If there aren't any boids in these cells, the start and end owuld be the same. This prevents the loop below from running
                // However, if there ARE boids, we can iterate through them. `i` in this case would iterate through `rearrangedBoids`
                for (int j = start; j < end; j++) {
                    Boid3D otherBoid = boids[rearrangedBoids[j]];
                    Vector3 otherPosition = new Vector3(otherBoid.position[0], otherBoid.position[1], otherBoid.position[2]);
                    Vector3 otherVelocity = velocities[rearrangedBoids[j]];
                    float distance = Vector3.Distance(currentPosition, otherPosition);
                    if (distance > 0 && distance < gridCellSize) {
                        if (distance < innerRange) {
                            close += currentPosition - otherPosition;
                        }
                        center += otherPosition;
                        avgVel += otherVelocity;
                        neighbors += 1;
                    }
                }
            }
        }
        // Assuming we actually encountered any neighbors, we have to perform the necessary calcualtions
        if (neighbors > 0) {
            center /= neighbors;
            avgVel /= neighbors;
            Vector3 centerEffect = (center - currentPosition) * (cohesionFactor * deltaTime);
            Vector3 alignmentEffect = (avgVel - currentVelocity) * (alignmentFactor * deltaTime);
            velocities[i] += centerEffect;
            velocities[i] += alignmentEffect;
        }
        Vector3 closeEffect = close * (separationFactor * deltaTime);
        velocities[i] += closeEffect;
    }

    private void LimitSpeed(ref Boid3D boid, int i) {
        Vector3 currentVelocity = velocities[i];
        float speed = currentVelocity.magnitude;
        if (speed > maxSpeed) currentVelocity = (currentVelocity / speed) * maxSpeed;
        else if (speed < minSpeed) currentVelocity = (currentVelocity / speed) * minSpeed;
        velocities[i] = currentVelocity;
    }

    private void KeepInBounds(ref Boid3D boid, int i) {
        float vX = velocities[i].x;
        float vY = velocities[i].y;
        float vZ = velocities[i].z;
        float deltaTime = (dt < 0) ? Time.deltaTime : dt;

        if (boid.position[0] < -bounds.x/2f) vX += deltaTime * turnSpeed;
        else if (boid.position[0] > bounds.x/2f) vX -= deltaTime * turnSpeed;

        if (boid.position[1] < -bounds.y/2f) vY += deltaTime * turnSpeed;
        else if (boid.position[1] > bounds.y/2f) vY -= deltaTime * turnSpeed;
    
        if (boid.position[2] < -bounds.z/2f) vZ += deltaTime * turnSpeed;
        else if (boid.position[2] > bounds.z/2f) vZ -= deltaTime * turnSpeed;

        velocities[i] = new Vector3(vX,vY,vZ);
    }

    private void OnDestroy() {
        ReleaseBuffers();
    }

    private void ReleaseBuffers() {
        gridBuffer.Release();
        boidsBuffer.Release();
        boidVelocitiesBuffer.Release();
        boidOffsetsBuffer.Release();
        rearrangedBoidsBuffer.Release();
        gridOffsetBuffer.Release();
        gridSumsBuffer1.Release();
        gridSumsBuffer2.Release();
    }

    private void DebugBufferInt(string debugText, int debugSize, ComputeBuffer b, bool getTotal = false) {
        int[] temp = new int[debugSize];
        b.GetData(temp);
        string top = "", bottom = "";
        int count = 0;
        for(int i = 0; i < debugSize; i++) {
            top += $"{i}\t|";
            bottom += $"{temp[i]}\t|";
            count += temp[i];
        }
        Debug.Log($"{debugText}:\n{top}\n{bottom}");
        if (getTotal) Debug.Log($"Total: {count}");
        temp = null;
    }
    private void DebugBufferFloat(string debugText, int debugSize, ComputeBuffer b) {
        float[] temp = new float[debugSize];
        b.GetData(temp);
        string top = "", bottom = "";
        for(int i = 0; i < debugSize; i++) {
            top += $"{i}\t|";
            bottom += $"{temp[i]}\t|";
        }
        Debug.Log($"{debugText}:\n{top}\n{bottom}");
        temp = null;
    }
    private void DebugBufferFloat3(string debugText, int debugSize, ComputeBuffer b) {
        float3[] temp = new float3[debugSize];
        b.GetData(temp);
        string top = "", bottom = "";
        for(int i = 0; i < debugSize; i++) {
            top += $"{i}\t|";
            bottom += $"{temp[i]}\t|";
        }
        Debug.Log($"{debugText}:\n{top}\n{bottom}");
        temp = null;
    }
    private void DebugBufferBoid3D(string debugText, int debugSize, ComputeBuffer b) {
        Boid3D[] temp = new Boid3D[debugSize];
        b.GetData(temp);
        string top = "";
        string bottom = "";
        for(int i = 0; i < debugSize; i++) {
            top += $"{i}\t|";
            bottom += $"{temp[i].position}\t|";
        }
        Debug.Log($"{debugText} positions:\n{top}\n{bottom}");
        temp = null;
    }
}
