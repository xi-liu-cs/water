using System.Collections;
using System.Collections.Generic;
using UnityEngine;

class sea_world_color : MonoBehaviour
{
    public Material mat;
    [Range(0, 1)]
    public float fog_distance_multiplier = 1;
    public Vector4 shader_param;
    mesh_generator mesh_gen;
    Camera cam;
    public Gradient gradient;
    public float normal_offset_weight;
    Texture2D texture;
    int texture_resolution = 50;

    void init()
    {
        if(texture == null || texture.width != texture_resolution)
            texture = new Texture2D(texture_resolution, 1, TextureFormat.RGBA32, false);
    }

    void Update()
    {
        init();
        update_texture();
        if(mesh_gen == null) mesh_gen = FindObjectOfType<mesh_generator>();
        if(cam == null) cam = FindObjectOfType<Camera>();
        mat.SetTexture("ramp", texture);
        mat.SetVector("param", shader_param);
        RenderSettings.fogColor = cam.backgroundColor;
        RenderSettings.fogEndDistance = mesh_gen.view_distance * fog_distance_multiplier;
    }

    void update_texture()
    {
        if(gradient != null)
        {
            Color[] colors = new Color[texture.width];
            for(int i = 0; i < texture_resolution; ++i)
            {
                Color gradient_col = gradient.Evaluate(i / (texture_resolution - 1));
                colors[i] = gradient_col;
            }
            texture.SetPixels(colors);
            texture.Apply();
        }
    }
}