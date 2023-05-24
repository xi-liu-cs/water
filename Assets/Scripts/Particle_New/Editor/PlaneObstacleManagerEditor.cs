using UnityEngine;
#if UNITY_EDITOR
    using UnityEditor;
#endif

[CustomEditor(typeof(PlaneObstacleManager))]
public class PlaneObstacleManagerEditor : Editor
{
    public override void OnInspectorGUI() {
        PlaneObstacleManager manager = (PlaneObstacleManager)target;
        DrawDefaultInspector();
        
        if (GUILayout.Button("Preprocess Planes")) {
            manager.PreprocessPlanes();
        }
        if (GUILayout.Button("Reorder Vertices")) {
            manager.ReorderVertices();
        }
        if (GUILayout.Button("Check If Intersecting")) {
            manager.DebugCheckIfIntersecting();
        }
        /*
        if(GUILayout.Button("Update Buffers")) {
            manager.UpdateBuffers();
        }
        if(GUILayout.Button("Dispose Buffers")) {
            manager.DisposeBuffers();
        }
        */
    }
}
