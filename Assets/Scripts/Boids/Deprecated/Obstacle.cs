using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[System.Serializable]
public struct ObstacleContainer
{
    public Vector3 position;
    public Vector3 normal;
    public Vector3 size;
    public Vector3 x;
    public Vector3 y;

    public static int GetSize => sizeof(float) * 15;
}

[ExecuteInEditMode]
public class Obstacle : MonoBehaviour
{
    public ObstacleContainer obstacle;

	private void Update() {
		if (transform.hasChanged) {
            obstacle.position = transform.position;
            obstacle.normal = transform.up;
            obstacle.size = 10f * new Vector3(
                transform.lossyScale.z, 
                transform.lossyScale.x, 
                transform.lossyScale.y / 5
            );

            obstacle.x = transform.forward;
            obstacle.y = transform.right;
		}
	}
}
