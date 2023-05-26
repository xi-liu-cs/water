using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

[RequireComponent(typeof(MeshFilter))]
[RequireComponent(typeof(MeshRenderer))]
[ExecuteInEditMode]
public class PlaneObstacle : MonoBehaviour
{

    [System.Serializable]
    public struct Obs {
        public int index;
        // These are all in world scale!
        public float3 position;
        public float4 rotation;
        public float3 scale;
        public int verticesStart, verticesCount;
        public int hasChanged;
    }
    [System.Serializable]
    public struct ObsPlane {
        public int obstacleIndex;
        // these are all in local scale!
        public float3 vertex1, vertex2, vertex3, centroid;
        public float3 normalVector;
    }
    
    private MeshFilter meshFilter;
    private  MeshRenderer meshRenderer;

    public Obs obs;
    public List<float3> vertices = new List<float3>();   // local scale
    public List<float3> normals = new List<float3>();   // local scale

    public int[] oldToNewVertMap;
    public List<ObsPlane> obstaclePlanes = new List<ObsPlane>();
    public List<int> vertexTriangleMap;
    public int3[] vertexTriangleCount;

    public bool awaitInitialization = false;

    private void Awake() {
        meshFilter = GetComponent<MeshFilter>();
        meshRenderer = GetComponent<MeshRenderer>();

        obs = new Obs();
        
        if (!awaitInitialization) Initialize();
    }

    public void Initialize(int index = 0, int verticesStart = 0) {
        // This step is to initialize details about htis obstacle that are pertinent to doing calculations such as:
        //  1) converting local vertices to world scale
        //  2) Identifying which vertices are associated with this particle inside of a global `vertices` array
        obs.index = index;
        UpdateObs();
        // Handle mesh vertices stuff
        obs.verticesStart = verticesStart;
        GenerateVertices();
        obs.verticesCount = vertices.Count;        
        // Now we need to generate all the planes, one for each mesh triangle.
        GeneratePlanes();
    }

    public void GenerateVertices() {
        // Because some meshes are inherently stupidly written, we need to create a new set of vertices in local space.
        // For example, the default Cube mesh in Unity has 24 vertices despite having 12 triangles. Why? I dunno. 
        // But what's important is that we take into account for these scenarios.

        // Initialize some local-global variables
        Mesh mesh = meshFilter.sharedMesh;
        int vertIndex = 0;

        // First, we initialize a map so that if we call any verts from the mesh, we'll have a map to the new vertices map        
        oldToNewVertMap = new int[mesh.vertices.Length];
        // Second, we initialize our new vertex list and our normals list
        vertices = new List<float3>();
        normals = new List<float3>();
        List<int> normalCounts = new List<int>();
        
        // We iterate through the existing vertices array
        for(int i = 0; i < mesh.vertices.Length; i++) {
            // Get the float3
            float3 newVert = new(
                mesh.vertices[i].x,
                mesh.vertices[i].y,
                mesh.vertices[i].z
            );
            float3 norm = new(
                mesh.normals[i].x,
                mesh.normals[i].y,
                mesh.normals[i].z
            );
            // If it doesn't contain our vert, we add it
            if (!vertices.Contains(newVert)) {
                vertices.Add(newVert);
                normals.Add(norm);
                vertIndex = vertices.IndexOf(newVert);
                normalCounts.Add(1);
            } else {
                // Since this one does exist, we merely grab the vert index and update the normal associated with it
                vertIndex = vertices.IndexOf(newVert);
                normals[vertIndex] = normals[vertIndex] + norm;
                normalCounts[vertIndex] = normalCounts[vertIndex] + 1;
            }
            // Get the index, then set the mapping
            oldToNewVertMap[i] = vertIndex;
        }

        // As a last step, we normalize the normals
        for(int i = 0; i < normals.Count; i++) {
            normals[i] = normals[i] / (float)normalCounts[i];
        }
    }

    public void GeneratePlanes() {
        // Initialize some local-global variables
        Mesh mesh = meshFilter.sharedMesh;
        float3 normalVector;
        ObsPlane newPlane;
        
        // Initialize the ObsPlane list
        obstaclePlanes = new List<ObsPlane>();
        int triangleIndex;

        // Initialize a temporary dictionary where keys = vertices, and values = lists of ints
        Dictionary<int, List<int>> vertexToTriangleList = new Dictionary<int, List<int>>();

        // Iterate through all triangles
        Debug.Log($"Should have {mesh.triangles.Length / 3} planes...");
        for(int t = 0; t < mesh.triangles.Length; t+=3) {
            // Initialize new ObsPlane struct
            newPlane = new ObsPlane();
            newPlane.obstacleIndex = obs.index;

            // Set the vertices
            newPlane.vertex1 = vertices[oldToNewVertMap[mesh.triangles[t]]];
            newPlane.vertex2 = vertices[oldToNewVertMap[mesh.triangles[t+1]]];
            newPlane.vertex3 = vertices[oldToNewVertMap[mesh.triangles[t+2]]];

            // Calculate centroid based on vertices
            newPlane.centroid = (newPlane.vertex1 + newPlane.vertex2 + newPlane.vertex3) / 3f;
            //centroid = (mesh.vertices[mesh.triangles[t]] + mesh.vertices[mesh.triangles[t+1]] + mesh.vertices[mesh.triangles[t+2]])/3f;
            
            // Calculate the normal vector
            newPlane.normalVector = (normals[oldToNewVertMap[mesh.triangles[t]]] + normals[oldToNewVertMap[mesh.triangles[t+1]]] + normals[oldToNewVertMap[mesh.triangles[t+2]]]) / 3f;
            // normalVector = (mesh.normals[mesh.triangles[t]] + mesh.normals[mesh.triangles[t+1]] + mesh.normals[mesh.triangles[t+2]])/3f;
            
            // Add the plane to the list of planes
            obstaclePlanes.Add(newPlane);
            triangleIndex = obstaclePlanes.IndexOf(newPlane);

            // For each vertex, we add the triangle plane's index
            if(!vertexToTriangleList.ContainsKey(oldToNewVertMap[mesh.triangles[t]])) {
                vertexToTriangleList.Add(oldToNewVertMap[mesh.triangles[t]], new List<int>());
            }
            vertexToTriangleList[oldToNewVertMap[mesh.triangles[t]]].Add(triangleIndex);
            if(!vertexToTriangleList.ContainsKey(oldToNewVertMap[mesh.triangles[t+1]])) {
                vertexToTriangleList.Add(oldToNewVertMap[mesh.triangles[t+1]], new List<int>());
            }
            vertexToTriangleList[oldToNewVertMap[mesh.triangles[t+1]]].Add(triangleIndex);
            if(!vertexToTriangleList.ContainsKey(oldToNewVertMap[mesh.triangles[t+2]])) {
                vertexToTriangleList.Add(oldToNewVertMap[mesh.triangles[t+2]], new List<int>());
            }
            vertexToTriangleList[oldToNewVertMap[mesh.triangles[t+2]]].Add(triangleIndex);
        }
        /*
        Debug.Log($"Ended up with {obstaclePlanes.Count} planes...");
        for(int i = 0; i < vertexToTriangleList.Length; i++) {
            string s = "";
            foreach(int ind in vertexToTriangleList[i]) {
                s += ind.ToString() + ",";
            }
            Debug.Log($"{i}: {s}");
        }
        */

        // We now need to generate vertexTriangleMap and vertexTriangleCount
        vertexTriangleMap = new List<int>();
        vertexTriangleCount = new int3[vertices.Count];
        for(int i = 0; i < vertices.Count; i++) {
            vertexTriangleCount[i] = new(obs.index, vertexTriangleMap.Count, vertexToTriangleList[i].Count);
            vertexTriangleMap.AddRange(vertexToTriangleList[i]);
        }
    }
    
    private void UpdateObs() {
        obs.position = new(transform.position.x, transform.position.y, transform.position.z);
        obs.rotation = new(transform.rotation.x, transform.rotation.y, transform.rotation.z, transform.rotation.w);
        obs.scale = new(transform.lossyScale.x, transform.lossyScale.y, transform.lossyScale.z);
        if (transform.hasChanged) {
            Debug.Log($"{gameObject.name} has changed!");
            obs.hasChanged = 1;
            transform.hasChanged = false;
        } else {
            obs.hasChanged = 0;
        }
    }

    private void Update() {
        UpdateObs();
    }

    /*
    
    
    public enum ObstacleType {
        Static,
        OnTransformChange
    }

    

    public Obs obs;
    public ObstacleType obstacleType = ObstacleType.Static;
    public List<ObsPlane> obstaclePlanes = new List<ObsPlane>();
    public float dotBetweenParticleAndNormal = 0f;

    public Transform particle = null;
    public bool isIntersecting = false;
    public Vector3 projectionPoint = Vector3.zero;
    public Vector3 worldCentroid = Vector3.zero;

    [SerializeField] private bool printPlaneGizmos = false;

    void OnDrawGizmos() {
        if (obstaclePlanes.Count == 0) return;
        Gizmos.color = Color.black;
        Gizmos.DrawSphere(worldCentroid, 0.025f);
        Gizmos.color = Color.yellow;
        Gizmos.DrawSphere(projectionPoint, 0.01f);
        if (!printPlaneGizmos) return;
        int max = Mathf.Min(obstaclePlanes.Count,1000);
        for(int i = 0; i < max; i++) {
            ObsPlane p = obstaclePlanes[i];
            Gizmos.color = (isIntersecting) ? Color.yellow : Color.red;
            Gizmos.DrawSphere(LocalPointToWorldPoint(obs, p.centroid),0.01f);
            Gizmos.color = Color.blue;
            Gizmos.DrawLine(LocalPointToWorldPoint(obs, p.centroid), LocalPointToWorldPoint(obs, p.centroid + p.normalVector));
        }
    }

    public void Initialize() {
        
    }

    
    private void Awake() {
        meshFilter = GetComponent<MeshFilter>();
        meshRenderer = GetComponent<MeshRenderer>();
        obs = new Obs();
        UpdateObs();
        GeneratePlanes();
    }

    private void Update() {
        // Update `obs` if our transform has changed in some way
        if (transform.hasChanged) {
            UpdateObs();
            if (obstacleType != ObstacleType.Static) GeneratePlanes();
            transform.hasChanged = false;
        }

        // Exit early if we don't have a particle to begin with, or if we don't have any obstacle planes to consider
        if (particle == null || obstaclePlanes.Count == 0) return;

        // First, find closest centroid to current target
        ObsPlane closest = obstaclePlanes[FindClosestPlane(particle.position)];

        // Secondly, convert centroid to world space
        worldCentroid = LocalPointToWorldPoint(obs, closest.centroid);

        // Secondly, Calculate dot product between normal vector and vector b/w centroid and target
        Vector3 targetVector = (particle.position - worldCentroid).normalized;
        dotBetweenParticleAndNormal = Vector3.Dot(targetVector, LocalVectorToWorldVector(obs, closest.normalVector));
        
        // Thirdly, calculate projection onto plane
        projectionPoint = ClosestPointOnPlane(
            worldCentroid, 
            LocalVectorToWorldVector(obs, closest.normalVector), 
            particle.position
        );
        
        float d = DistanceFromPlane(worldCentroid, LocalVectorToWorldVector(obs, closest.normalVector), particle.position);

        // Finally, check if intersecting or not
        isIntersecting = dotBetweenParticleAndNormal <= 0f 
            && ObstacleHelper.PointInTriangle(
                projectionPoint, 
                LocalPointToWorldPoint(obs, closest.vertex1), 
                LocalPointToWorldPoint(obs, closest.vertex2), 
                LocalPointToWorldPoint(obs, closest.vertex3)
        );
    }

    

    

    private int FindClosestPlane(Vector3 t) {
        int closest = 0;
        float distance, closestDistance = Mathf.Infinity;
        for(int i = 0; i < obstaclePlanes.Count; i++) {
            distance = Vector3.Distance(t, transform.TransformPoint(obstaclePlanes[i].centroid));
            if (distance < closestDistance) {
                closestDistance = distance;
                closest = i;
            }
        }
        return closest;
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
    public static Vector3 LocalPointToWorldPoint(Obs obstacle, Vector3 localPoint) {
        Vector3 s = new Vector3(localPoint.x * obstacle.scale.x, localPoint.y * obstacle.scale.y, localPoint.z * obstacle.scale.z);
        return RotMultVec3(obstacle.rotation, s) + obstacle.position;
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
    public static Vector3 LocalVectorToWorldVector(Vector3 pos, Vector4 rot, Vector3 scale, Vector3 v) {
        Vector3 start = LocalPointToWorldPoint(pos, rot, scale, Vector3.zero);
        Vector3 end = LocalPointToWorldPoint(pos, rot, scale, v);
        return end - start;
    }
    public static Vector3 LocalVectorToWorldVector(Obs obstacle, Vector3 v) {
        Vector3 start = LocalPointToWorldPoint(obstacle, Vector3.zero);
        Vector3 end = LocalPointToWorldPoint(obstacle, v);
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
    */
}
