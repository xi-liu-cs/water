using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(MeshRenderer))]
[RequireComponent(typeof(MeshFilter))]
public class MeshManager : MonoBehaviour
{

    struct Triangle {
        public Vector3 a;
        public Vector3 b;
        public Vector3 c;

        public Vector3 this [int i] {
            get {
                switch (i) {
                    case 0:
                        return a;
                    case 1:
                        return b;
                    default:
                        return c;
                }
            }
        }
    }

    [Header("== REFERENCES ==")]
    public ParticleManager particleManager = null;
    public Material fluidMaterial = null;
    private MeshRenderer fluidMeshRenderer = null;
    private Mesh fluidMesh = null;
    public ComputeShader densityShader = null;
    public ComputeShader triangleShader = null;

    [Header("== WORLD CONFIGURATIONS ==")]
    public float[] origin = new float[3] { 0f, 0f, 0f };
    public Vector3 originVector3 => new Vector3(origin[0], origin[1], origin[2]);
    private float[] bounds = new float[3] {100f,100f,100f}; // Set based on voxelSize and numVoxelsPerAxis
    private Vector3 boundsVector3 => new Vector3(bounds[0], bounds[1], bounds[2]);

    [Header("== Marching Cubes Configurations ==")]
    public float voxelSize = 0.1f;
    public int[] numPointsPerAxis = new int[3] { 50, 50, 50 };
    private float[] pointGridCellSizes;
    private int[] numVoxelsPerAxis;
    private int numVoxels;
    public float smoothingRadius;
    public float radius2 => Mathf.Pow(smoothingRadius,2);
    public float radius3 => Mathf.Pow(smoothingRadius,3);
    private int numTriangles;
    Triangle[] triangles;
    public int maxParticlesPerVoxel = 10;
    private int maxTriangleCount;
    private int numPoints;
    public float isoLevel = 8f;

    [Header("== GPU CONFIGURATIONS  ==")]
    const int threadGroupSize = 8;
    private int[] numThreadsPerAxis;
    int size_property = Shader.PropertyToID("size");
    int particle_buffer_property = Shader.PropertyToID("particle_buffer");

    [Header("== GPU KERNELS ==")]
    private int clear_cube_corner_neighbor_tracker_kernel;
    private int compute_neighbor_list_kernel;
    private int compute_density_kernel;
    private int compute_triangles_kernel;

    [Header("== GPU BUFFERS ==")]
    private ComputeBuffer pointsBuffer;
    private ComputeBuffer pointDensitiesBuffer;
    private ComputeBuffer cubeCornerNeighborListBuffer;
    private ComputeBuffer cubeCornerNeighborTrackerBuffer;
    private ComputeBuffer triangleBuffer;
    private ComputeBuffer triangleCountBuffer;

    [Header("== DEBUGGING CONFIGURATIONS ==")]
    public bool gizmos_densities = false;
    public bool gizmos_points = false;

    void OnDrawGizmos() {
        int[] voxelsPerAxis = new int[3] { numPointsPerAxis[0]-1, numPointsPerAxis[1]-1, numPointsPerAxis[2]-1 };
        Vector3 gizmos_bounds = new Vector3(
            voxelsPerAxis[0] * voxelSize,
            voxelsPerAxis[1] * voxelSize,
            voxelsPerAxis[2] * voxelSize
        );
        float[] gizmos_cellSizes = new float[3] {
            gizmos_bounds.x / voxelsPerAxis[0],
            gizmos_bounds.y / voxelsPerAxis[1],
            gizmos_bounds.z / voxelsPerAxis[2]
        };

        Gizmos.color = Color.red;
        Gizmos.DrawWireCube(originVector3, gizmos_bounds);

        Gizmos.color = Color.yellow;
        if (Application.isPlaying) {
            if (!gizmos_densities) return;
            int[] temp_numbers = new int[numPoints];
            cubeCornerNeighborTrackerBuffer.GetData(temp_numbers);
            for(int x = 0; x < numPointsPerAxis[0]; x++) {
                for(int y = 0; y < numPointsPerAxis[1]; y++) {
                    for(int z = 0; z < numPointsPerAxis[2]; z++) {
                        int proj_index = GetProjectedGridIndexFromXYZ(x,y,z);
                        float num = Mathf.Clamp(Mathf.FloorToInt(temp_numbers[proj_index]/1f),0f,1f);
                        Gizmos.color = new Color(1f,1f,0f,num);
                        Vector3 pos = new Vector3(
                            originVector3.x - (gizmos_bounds[0]/2f) + (x * gizmos_cellSizes[0]),
                            originVector3.y - (gizmos_bounds[1]/2f) + (y * gizmos_cellSizes[1]),
                            originVector3.z - (gizmos_bounds[2]/2f) + (z * gizmos_cellSizes[2])
                        );
                        Gizmos.DrawSphere(pos, 1f);
                    }
                }
            }
        } else {
            if (!gizmos_points) return;
            for(int x = 0; x < numPointsPerAxis[0]; x++) {
                for(int y = 0; y < numPointsPerAxis[1]; y++) {
                    for(int z = 0; z < numPointsPerAxis[2]; z++) {
                        Vector3 pos = new Vector3(
                            originVector3.x - (gizmos_bounds[0]/2f) + (x * gizmos_cellSizes[0]),
                            originVector3.y - (gizmos_bounds[1]/2f) + (y * gizmos_cellSizes[1]),
                            originVector3.z - (gizmos_bounds[2]/2f) + (z * gizmos_cellSizes[2])
                        );
                        Gizmos.DrawSphere(pos, 1f);
                    }
                }
            }
        }
    }

    int GetProjectedGridIndexFromXYZ(int x, int y, int z) {
        return x + (numPointsPerAxis[0] * y) + (numPointsPerAxis[0] * numPointsPerAxis[1] * z);
    }

    private void Awake() {
        particleManager.Initialize();
        fluidMeshRenderer = GetComponent<MeshRenderer>();
        fluidMesh = GetComponent<MeshFilter>().mesh;   

        InitializeVariables();
        InitializeShaderVariables();
        InitializeShaderKernels();
        InitializeShaderBuffers();        
    }

    private void InitializeVariables() {
        numPoints = numPointsPerAxis[0] * numPointsPerAxis[1] * numPointsPerAxis[2];
        numVoxelsPerAxis = new int[3] { numPointsPerAxis[0]-1, numPointsPerAxis[1]-1, numPointsPerAxis[2]-1 };
        numVoxels = numVoxelsPerAxis[0] * numVoxelsPerAxis[1] * numVoxelsPerAxis[2];
        maxTriangleCount = numVoxels * 5;

        bounds = new float[3] {
            numVoxelsPerAxis[0] * voxelSize,
            numVoxelsPerAxis[1] * voxelSize,
            numVoxelsPerAxis[2] * voxelSize
        };

        pointGridCellSizes = new float[3] {
            bounds[0] / numVoxelsPerAxis[0],
            bounds[1] / numVoxelsPerAxis[1],
            bounds[2] / numVoxelsPerAxis[2]
        };

        numThreadsPerAxis = new int[3] {
            Mathf.CeilToInt(numVoxelsPerAxis[0] / (float)threadGroupSize),
            Mathf.CeilToInt(numVoxelsPerAxis[1] / (float)threadGroupSize),
            Mathf.CeilToInt(numVoxelsPerAxis[2] / (float)threadGroupSize),
        };

        Debug.Log($"Bounds: ({bounds[0]},{bounds[1]},{bounds[2]})");
        Debug.Log($"Number of total points: {numPoints}");
        Debug.Log($"Number of total voxels: {numVoxels}");
        Debug.Log($"Calculated Grid Cell Sizes: ({pointGridCellSizes[0]},{pointGridCellSizes[1]},{pointGridCellSizes[2]})");
    }

    private void InitializeShaderVariables() {
        densityShader.SetFloats("origin", origin);
        densityShader.SetFloat("radius", smoothingRadius);
        densityShader.SetFloat("radius2", radius2);
        densityShader.SetFloat("radius3", radius3);
        densityShader.SetFloats("bounds", bounds);
        densityShader.SetInt("max_particles_per_cube", maxParticlesPerVoxel);
        densityShader.SetFloats("pointGridCellSizes", pointGridCellSizes);
        densityShader.SetInts("dimension", numPointsPerAxis);
        // Retrieved from ParticleManager
        densityShader.SetInt("numParticles", particleManager.numParticles);
        densityShader.SetFloat("particleRenderSize", particleManager.particleRenderSize);
    }

    private void InitializeShaderKernels() {
        clear_cube_corner_neighbor_tracker_kernel = densityShader.FindKernel("clear_cube_corner_neighbor_tracker");
        compute_neighbor_list_kernel = densityShader.FindKernel("compute_neighbor_list");
        compute_density_kernel = densityShader.FindKernel("compute_density");
    }
    
    private void InitializeShaderBuffers() {
        pointsBuffer = new ComputeBuffer(numPoints, sizeof(float) * 3);
        pointDensitiesBuffer = new ComputeBuffer(numPoints, sizeof(float));
        cubeCornerNeighborListBuffer = new ComputeBuffer(numPoints * maxParticlesPerVoxel, sizeof(int));
        cubeCornerNeighborTrackerBuffer = new ComputeBuffer(numPoints, sizeof(int));

        triangleBuffer = new ComputeBuffer(maxTriangleCount, sizeof(float)*9, ComputeBufferType.Append);
        triangleCountBuffer = new ComputeBuffer(1, sizeof(int), ComputeBufferType.Raw);

        densityShader.SetBuffer(clear_cube_corner_neighbor_tracker_kernel, "cube_corner_neighbor_tracker", cubeCornerNeighborTrackerBuffer);

        densityShader.SetBuffer(compute_neighbor_list_kernel, "particles", particleManager.particle_buffer);
        densityShader.SetBuffer(compute_neighbor_list_kernel, "cube_corner_neighbor_list", cubeCornerNeighborListBuffer);
        densityShader.SetBuffer(compute_neighbor_list_kernel, "cube_corner_neighbor_tracker", cubeCornerNeighborTrackerBuffer);

        densityShader.SetBuffer(compute_density_kernel, "particles", particleManager.particle_buffer);
        densityShader.SetBuffer(compute_density_kernel, "voxel_density", pointDensitiesBuffer);
        densityShader.SetBuffer(compute_density_kernel, "points", pointsBuffer);
        densityShader.SetBuffer(compute_density_kernel, "cube_corner_neighbor_list", cubeCornerNeighborListBuffer);
        densityShader.SetBuffer(compute_density_kernel, "cube_corner_neighbor_tracker", cubeCornerNeighborTrackerBuffer);

        triangleShader.SetBuffer(0, "triangles", triangleBuffer);
        triangleShader.SetInts("n_point_per_axis", numPointsPerAxis);
        triangleShader.SetFloat("isolevel", isoLevel);
        triangleShader.SetBuffer(0, "voxel_density", pointDensitiesBuffer);
        triangleShader.SetBuffer(0, "particles", particleManager.particle_buffer);
        triangleShader.SetBuffer(0, "points", pointsBuffer);


        fluidMaterial.SetFloat(size_property, particleManager.particleRenderSize);
        fluidMaterial.SetBuffer(particle_buffer_property, particleManager.particle_buffer);
        fluidMeshRenderer.material = fluidMaterial;
    }

    private void Update() {
        particleManager.UpdateParticles();
        UpdateDensities();
        UpdateTriangles();
        UpdateMesh();
    }

    private void UpdateDensities() {
        densityShader.Dispatch(
            clear_cube_corner_neighbor_tracker_kernel, 
            numThreadsPerAxis[0], 
            numThreadsPerAxis[1], 
            numThreadsPerAxis[2]
        );
        densityShader.Dispatch(
            compute_neighbor_list_kernel, 
            particleManager.NUM_THREADS_FOR_PARTICLES, 
            1,
            1
        );
        densityShader.Dispatch(
            compute_density_kernel, 
            numThreadsPerAxis[0], 
            numThreadsPerAxis[1], 
            numThreadsPerAxis[2]
        );
    }

    private void UpdateTriangles() {
        triangleBuffer.SetCounterValue(0);
        triangleShader.Dispatch(0, numThreadsPerAxis[0], numThreadsPerAxis[1], numThreadsPerAxis[2]);
        ComputeBuffer.CopyCount(triangleBuffer, triangleCountBuffer, 0);
        int[] triCountArray = { 0 };
        triangleCountBuffer.GetData(triCountArray);
        numTriangles = triCountArray[0];
        triangles = new Triangle[numTriangles];
        triangleBuffer.GetData(triangles, 0, 0, numTriangles);
        Debug.Log($"Number of triangles / total max count: {(float)numTriangles / (float)maxTriangleCount}");
    }

    private void UpdateMesh() {
        fluidMesh.Clear();
        Vector3[] vertices = new Vector3[numTriangles * 3];
        int[] meshTriangles = new int[numTriangles * 3];

        for(int i = 0; i < numTriangles; i++) {
            for(int j = 0; j < 3; j++) {
                meshTriangles[i * 3 + j] = i * 3 + j; // assign index
                vertices[i * 3 + j] = triangles[i][j];
            }
        }
        fluidMesh.vertices = vertices;
        fluidMesh.triangles = meshTriangles;
        fluidMesh.RecalculateNormals();
    }

    void OnDestroy() {
        pointsBuffer.Dispose();
        pointDensitiesBuffer.Dispose();
        cubeCornerNeighborListBuffer.Dispose();
        cubeCornerNeighborTrackerBuffer.Dispose();
        triangleBuffer.Dispose();
        triangleCountBuffer.Dispose();
    }
}
