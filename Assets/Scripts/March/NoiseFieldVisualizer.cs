using System;
using UnityEngine;

namespace MarchingCubes {

sealed class NoiseFieldVisualizer : MonoBehaviour
{
    #region Editable attributes

    [SerializeField] Vector3Int _dimension = new Vector3Int(64, 32, 64);
    [SerializeField] float _gridScale = 4.0f / 64;
    [SerializeField] int _triangleBudget = 65536;
    [SerializeField] float _targetValue = 0;

    #endregion

    #region Project asset references

    [SerializeField] ComputeShader _volumeCompute = null;
    [SerializeField] ComputeShader _builderCompute = null;
    [SerializeField] fluid_gpu fluid_cs;
    [SerializeField] Material material;

    #endregion

    #region Private members

    int VoxelCount => _dimension.x * _dimension.y * _dimension.z;

    ComputeBuffer _voxelBuffer;
    MeshBuilder _builder;
    int noise_field_visualizer_kernel;

    #endregion

    #region MonoBehaviour implementation

    void Start()
    {
        fluid_cs.Awake();
        gameObject.AddComponent<MeshFilter>();
        gameObject.AddComponent<MeshRenderer>();
        /* _voxelBuffer = new ComputeBuffer(VoxelCount, sizeof(float)); */
        noise_field_visualizer_kernel = _volumeCompute.FindKernel("NoiseFieldGenerator");
        /* _voxelBuffer = fluid_cs.noise_density_buffer; */
        _dimension = fluid_cs.dimension;
        _gridScale = 1f;
        _voxelBuffer = fluid_cs.density_buffer;
        _builder = new MeshBuilder(_dimension, _triangleBudget, _builderCompute, material);
    }

    void OnDestroy()
    {
        _voxelBuffer.Dispose();
        _builder.Dispose();
    }

    void Update()
    {
        fluid_cs.Update();
        /* float[] b = new float[1000];
        _voxelBuffer.GetData(b);
        for(int i = 0; i < 1000; ++i) Debug.Log(String.Format("b[{0}] = {1}", i, b[i])); */
        // Noise field update
        _volumeCompute.SetInts("dimension", _dimension);
        _volumeCompute.SetFloat("scale", _gridScale);
        _volumeCompute.SetFloat("time", Time.time);
        _volumeCompute.SetFloat("max_density", fluid_cs.max_density);
        _volumeCompute.SetBuffer(noise_field_visualizer_kernel, "particles", fluid_cs.particle_buffer);
        _volumeCompute.SetBuffer(noise_field_visualizer_kernel, "voxels", _voxelBuffer);
        _volumeCompute.DispatchThreads(noise_field_visualizer_kernel, _dimension);

        // Isosurface reconstruction
        _builder.BuildIsosurface(_voxelBuffer, _targetValue, _gridScale);
        GetComponent<MeshFilter>().sharedMesh = _builder.Mesh;
    }

    #endregion
}

} // namespace MarchingCubes
