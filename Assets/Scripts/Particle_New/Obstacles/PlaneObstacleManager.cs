using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PlaneObstacleManager : MonoBehaviour
{   
    // We store a global reference to this specific script via the `current` method.
    public static PlaneObstacleManager current;

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
}
