using System.Collections;
using System.Collections.Generic;
using UnityEngine;

class noise_density : density_generator
{
    [Header("noise")]
    public int seed,
    n_octave = 4;
    public float lacunarity = 2,
    persistence = .5f,
    noise_scale = 1,
    noise_weight = 1,
    floor_offset = 1,
    weight_multiplier = 1,
    hard_floor_height,
    hard_floor_weight;
    public bool close_edge;
    public Vector4 shader_param;

    public override ComputeBuffer generate(ComputeBuffer point_buffer, int n_point_per_axis, float bound_size, Vector3 world_bound, Vector3 center, Vector3 offset, float spacing)
    {
        buffer_release = new List<ComputeBuffer>();
        System.Random rand = new System.Random(seed);
        Vector3[] offsets = new Vector3[n_octave];
        float offset_range = 1000;
        for(int i = 0; i < n_octave; ++i)
        {
            float a = 2 * (float)rand.NextDouble() - 1;
            offsets[i] = new Vector3(a, a, a) * offset_range;
        }
        ComputeBuffer offset_buffer = new ComputeBuffer(offsets.Length, 3 * sizeof(float));
        offset_buffer.SetData(offsets);
        buffer_release.Add(offset_buffer);
        density_shader.SetVector("center", new Vector4(center.x, center.y, center.z));
        density_shader.SetInt("octave", Mathf.Max(1, n_octave));
        density_shader.SetFloat("lacunarity", lacunarity);
        density_shader.SetFloat("persistence", persistence);
        density_shader.SetFloat("noise_scale", noise_scale);
        density_shader.SetFloat("noise_weight", noise_weight);
        density_shader.SetBool("close_edge", close_edge);
        density_shader.SetBuffer(density_kernel, "offsets", offset_buffer);
        density_shader.SetFloat("floor_offset", floor_offset);
        density_shader.SetFloat("weight_multiplier", weight_multiplier);
        density_shader.SetFloat("hard_floor", hard_floor_height);
        density_shader.SetFloat("hard_floor_weight", hard_floor_weight);
        density_shader.SetVector("param", shader_param);
        return base.generate(point_buffer, n_point_per_axis, bound_size, world_bound, center, offset, spacing);
    }
}