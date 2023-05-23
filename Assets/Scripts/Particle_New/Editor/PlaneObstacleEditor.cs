using UnityEngine;
#if UNITY_EDITOR
    using UnityEditor;
#endif

[CustomEditor(typeof(PlaneObstacle))]
public class PlaneObstacleEditor : Editor
{
    public override void OnInspectorGUI() {
        PlaneObstacle obstacle = (PlaneObstacle)target;
        DrawDefaultInspector();
        /*
        if(GUILayout.Button("Generate Planes from Mesh Triangles")) {
            if (obstacle.GeneratePlanes()) Debug.Log("Planes generated");
        }
        */
    }
}
