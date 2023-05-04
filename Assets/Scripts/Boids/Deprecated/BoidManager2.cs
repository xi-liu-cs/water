using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using Random = Unity.Mathematics.Random;
using System.Runtime.InteropServices;

public class BoidManager2 : MonoBehaviour
{

    public struct BoidS {
        public Vector3 position;
        public Vector3 forward;
    } 

    public int numBoids = 32;
    private int boidCountPoT;
    public float3 dimensions = new (50.0f,50.0f,50.0f);
    private int prefixSumBlockSize = 32;

    public GraphicsBuffer boidsBuffer;
    public GraphicsBuffer boidsPrefixSumBuffer;

    public ComputeShader parallelShader;
    public ComputeShader steerShader;

    private int parallelMainKernel;
    private int steerMainKernel;

    public Transform boidTarget = null;

    //[Range(0f,0.5f)]
    //public float turnFactor = 0.2f;
    [Range(0f,0.01f)]
    public float alignmentFactor = 0.0005f;
    [Range(0f,0.1f)]
    public float separationFactor = 0.05f;
    public float targetFactor = 0.5f;
    //[Range(0f,0.1f)]
    //public float matchingFactor = 0.05f;
    //public Vector2 speedLimits = new Vector2(3f,6f);

    public float moveSpeed = 10f;

    void OnDrawGizmos() {
        if (Application.isPlaying) {
            Gizmos.color = Color.yellow;
            int bCount = boidsBuffer.count;
            BoidS[] boids = new BoidS[bCount];
            boidsBuffer.GetData(boids);
            for(int i = 0; i < bCount; i++) {
                Gizmos.DrawSphere(boids[i].position, 1f);
            }
        }
    }

    // Start is called before the first frame update
    void Start() {
        // Get closest minimum power of 2 greater than or equal to `numBoids`
        boidCountPoT = math.ceilpow2(numBoids);
        Debug.Log($"Boid Count: {numBoids} - Boud Count POT: {boidCountPoT}");

        // Prepping the boid buffer. Note that we use `boidCountPoT`
        InitializeBuffers();

        // Initialize the Shaders
        InitializeShaders();
    }

    private void InitializeBuffers() {
        var random = new Random(256);
        BoidS[] boidArray = new BoidS[numBoids];

        for(int i = 0; i < boidArray.Length; i++) {
            boidArray[i] = new BoidS {
                position = random.NextFloat3(-dimensions, dimensions),
                forward = math.rotate(random.NextQuaternionRotation(), Vector3.forward),
            };
        }
       
        boidsBuffer = new GraphicsBuffer(
            GraphicsBuffer.Target.Structured, 
            boidArray.Length, 
            Marshal.SizeOf<BoidS>()
        );
        boidsBuffer.SetData(boidArray);

        boidsPrefixSumBuffer = new GraphicsBuffer(
            GraphicsBuffer.Target.Structured, 
            boidCountPoT,
            Marshal.SizeOf<BoidS>()
        );
    }

    private void InitializeShaders() {
        // Initialize the parallel shader
        parallelShader.SetInt("numBoids",boidCountPoT);
        parallelMainKernel = parallelShader.FindKernel("main");

        // Initialize the steer shader
        steerMainKernel = steerShader.FindKernel("main");
        steerShader.SetBuffer(steerMainKernel, "boidsBuffer", boidsBuffer);
        steerShader.SetBuffer(steerMainKernel, "boidsPrefixSumBuffer", boidsPrefixSumBuffer);
        steerShader.SetInt("numBoids", boidCountPoT);
    }

    // Update is called once per frame
    void Update()
    {
        UpdateAggregation();
        UpdateSteering();
    }

    private void UpdateAggregation() {
        parallelShader.GetKernelThreadGroupSizes(parallelMainKernel, out var x, out var y, out var z);
        Debug.Log(x);
        parallelShader.SetBuffer(parallelMainKernel, "boidsBuffer", boidsBuffer);
        parallelShader.SetBuffer(parallelMainKernel, "boidsPrefixSumBuffer", boidsPrefixSumBuffer);
        for (var n = boidCountPoT; n >= prefixSumBlockSize; n /= prefixSumBlockSize)
        {
            parallelShader.Dispatch(parallelMainKernel, (int) (n / x), 1, 1);
            parallelShader.SetBuffer(parallelMainKernel, "boidsBuffer", boidsPrefixSumBuffer);
        }
    }

    private void UpdateSteering() {
        Vector3 boidTargetPos = boidTarget != null
            ? boidTarget.position
            : transform.position;
        steerShader.SetFloat("deltaTime", Time.deltaTime);
        steerShader.SetInt("numBoids", boidCountPoT);

        steerShader.SetFloat("separationWeight", separationFactor);
        steerShader.SetFloat("alignmentWeight", alignmentFactor);
        steerShader.SetFloat("targetWeight", targetFactor);
        steerShader.SetFloat("moveSpeed", moveSpeed);
        steerShader.SetVector("targetPosition", boidTargetPos);

        steerShader.GetKernelThreadGroupSizes(steerMainKernel, out var x, out var y, out var z);
        steerShader.Dispatch(steerMainKernel, (int) (boidCountPoT / x), 1, 1);
    }

    private void OnDestroy()
    {
        ReleaseBuffers();
    }

    private void ReleaseBuffers()
    {
        boidsBuffer.Dispose();
        boidsPrefixSumBuffer.Dispose();
    }
}
