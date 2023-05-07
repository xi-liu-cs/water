using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using System.Linq;

public class PointCloudObstacleManager : MonoBehaviour
{

    public List<PointCloudObstacle> obstacles;

    [SerializeField]
    private float3[] _points = new float3[0];
    public float3[] points { get => _points; set {} }
    public int numPoints = 0;
    public ComputeBuffer pointsBuffer;

    public bool waitForInitialization = false;

    // Start is called before the first frame update
    void Awake() {
        if (!waitForInitialization) Initialize();
    }

    public void Initialize() {
    }

    // Update is called once per frame
    void Update() {
        GetPoints();
        pointsBuffer = new ComputeBuffer(_points.Length, sizeof(float)*3);
    }

    public void GetPoints() {
        numPoints = 0;
        foreach(PointCloudObstacle obstacle in obstacles) {
            numPoints += obstacle.points.Count;
        }

        if (numPoints == 0) {
            _points = new float3[1] { new(0f,0f,0f) };
        }
        else {
            _points = new float3[numPoints];
            foreach(PointCloudObstacle obstacle in obstacles) {
                _points.Concat(obstacle.pointsF3);
            }
        }
    }
}
