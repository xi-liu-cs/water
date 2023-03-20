using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class fog_start : MonoBehaviour
{
    public bool AllowFog;
    private bool FogOn;
    public GameObject waterLevel;
    
    void Start()
    {
        RenderSettings.fogColor = new Vector4(0.15f, 0.25f, 0.45f, 1f);
        RenderSettings.fogDensity = 0.008f;
    }

    void Update()
    {/* SceneView.currentDrawingSceneView.camera.transform.position.y */ 
        if(transform.position.y > waterLevel.transform.position.y) /* when vehicle is above water */
            AllowFog = false;
        if(transform.position.y < waterLevel.transform.position.y) /* when vehicle is under water */
            AllowFog = true;
    }

    private void OnPreRender()
    {
        FogOn = RenderSettings.fog;
        RenderSettings.fog = AllowFog; /* if the variable Allowfog is true, it switches on the fog and similarly switches it off when the variable is false */
    }

    private void OnPostRender()
    {
        RenderSettings.fog = FogOn;
    }
}
