using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class fog_start : MonoBehaviour
{
    public bool AllowFog;

    private bool FogOn;

    public GameObject waterLevel;
    // Start is called before the first frame update
    void Start()
    {
        RenderSettings.fogColor = new Vector4(0.15f, 0.25f, 0.45f, 1f);
    }

    // Update is called once per frame
    void Update()
    {/* SceneView.currentDrawingSceneView.camera.transform.position.y */ 
        if(transform.position.y > waterLevel.transform.position.y) // When vehicle is above water, fog is turned off.
            AllowFog = false;
        if(transform.position.y < waterLevel.transform.position.y) // And its turned on when the vehicle goes under water.. 
            AllowFog = true;
    }

    private void OnPreRender()
    {
        FogOn = RenderSettings.fog;
        RenderSettings.fog = AllowFog; // If the variable Allowfog is true, it switches on the fog and similarly switches it off when the variable is false
    }

    private void OnPostRender()
    {
        RenderSettings.fog = FogOn;
    }
}
