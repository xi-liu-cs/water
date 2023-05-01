using UnityEngine;

namespace MarchingCubes {

sealed class VolumeDataVisualizer : MonoBehaviour
{
    #region Editable attributes

    [SerializeField] TextAsset _volumeData = null;
    [SerializeField] Vector3Int _dimension = new Vector3Int(256, 256, 113);
    [SerializeField] float _gridScale = 4.0f / 256;
    [SerializeField] int _triangleBudget = 65536 * 16;

    #endregion

    #region Project asset references

    [SerializeField] ComputeShader _converterCompute = null;
    [SerializeField] ComputeShader _builderCompute = null;
    [SerializeField] fluid_gpu fluid_cs;
    [SerializeField] Material material;
    #endregion

    #region Target isovalue

    public float TargetValue { get; set; } = 0.4f;
    float _builtTargetValue;

    #endregion

    #region Private members

    int VoxelCount => _dimension.x * _dimension.y * _dimension.z;

    ComputeBuffer _voxelBuffer;
    MeshBuilder _builder;
    int volume_data_visualizer_kernel;

    #endregion

    #region MonoBehaviour implementation

    void Start()
    {
        fluid_cs.Initialize();
        gameObject.AddComponent<MeshFilter>();
        gameObject.AddComponent<MeshRenderer>();
        /* _voxelBuffer = new ComputeBuffer(VoxelCount, sizeof(float)); */
        volume_data_visualizer_kernel = _converterCompute.FindKernel("NoiseFieldGenerator");
        _voxelBuffer = fluid_cs.density_buffer;
        _builder = new MeshBuilder(_dimension, _triangleBudget, _builderCompute, material);

        // Voxel data conversion (ushort -> float)
        /* using var readBuffer = new ComputeBuffer(VoxelCount / 2, sizeof(uint));
        readBuffer.SetData(_volumeData.bytes); */

        _converterCompute.SetInts("Dims", _dimension);
        /* _converterCompute.SetBuffer(0, "Source", readBuffer); */
        _converterCompute.SetBuffer(volume_data_visualizer_kernel, "Voxels", _voxelBuffer);
        _converterCompute.DispatchThreads(volume_data_visualizer_kernel, _dimension);
    }

    void OnDestroy()
    {
        _voxelBuffer.Dispose();
        _builder.Dispose();
    }

    void Update()
    {
        fluid_cs.UpdateParticles();
        // Rebuild the isosurface only when the target value has been changed.
        if (TargetValue == _builtTargetValue) return;

        _builder.BuildIsosurface(_voxelBuffer, TargetValue, _gridScale);
        GetComponent<MeshFilter>().sharedMesh = _builder.Mesh;

        _builtTargetValue = TargetValue;
    }

    #endregion
}

} // namespace MarchingCubes
