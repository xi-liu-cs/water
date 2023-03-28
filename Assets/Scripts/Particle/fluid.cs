using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class fluid : MonoBehaviour
{
    [Header("particle")]
    public Mesh particle_mesh;
    public Material material;
    public float mass = 4f,
    viscosity = 2.5f,
    particle_size = 4f,
    radius = 1f,
    gas_constant = 2000f,
    dt = 0.0008f,
    rest_density = 1f,
    damping = -0.5f;
    Vector3 g = new Vector3(0f, -9.81f, 0f);

    [Header("fluid")]
    public int n_particle = 1000,
    dimension = 10,
    maximum_particle_per_cell = 500;

    struct particle
    {
        Vector3 position;
        Vector4 color;
    }
}
