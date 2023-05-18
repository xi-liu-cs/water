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
        if(GUILayout.Button("Update Buffers")) {
            manager.UpdateBuffers();
        }
        if(GUILayout.Button("Dispose Buffers")) {
            manager.DisposeBuffers();
        }
    }
}
