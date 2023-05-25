using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.MLAgents;
using Unity.MLAgents.Sensors;
using Unity.MLAgents.Actuators;

public class rov : Agent
{
    Rigidbody rigid_body;
    public Transform target;
    public Vector3 initial_position;
    public float speed = 1f,
    sensitivity = 1f,
    force_multiplier = 10f;
    public float[] bound = {-50, 50, -50, 50, -50, 50};
    
    void Start()
    {
        rigid_body = GetComponent<Rigidbody>();
        initial_position = gameObject.transform.position;
    }

    bool is_outside(Vector3 position)
    {
        return position.x < bound[0] || position.x > bound[1]
        || position.y < bound[2] || position.y > bound[3]
        || position.z < bound[4] || position.z > bound[5];
    }
    
    public override void OnEpisodeBegin()
    {
        transform.position = initial_position;
    }

    public override void CollectObservations(VectorSensor sensor)
    {
        sensor.AddObservation(target.localPosition);
        sensor.AddObservation(this.transform.localPosition);
        sensor.AddObservation(rigid_body.velocity);
    }

    public override void OnActionReceived(ActionBuffers actions)
    {
        float dx = actions.ContinuousActions[0],
        dy = actions.ContinuousActions[1],
        dz = actions.ContinuousActions[2];
        transform.position += new Vector3(dx, dy, dz) * Time.deltaTime * speed;

        float distance_to_target = Vector3.Distance(this.transform.localPosition, target.localPosition);
        if(distance_to_target < 1f)
        {
            SetReward(1f);
            EndEpisode();
        }
        if(is_outside(this.transform.localPosition))
        {
            SetReward(-1f);
            EndEpisode();
        }
    }

    public override void Heuristic(in ActionBuffers actions_out)
    {/* in inspector behavior parameters script, change to heuristic only */
        ActionSegment<float> continuous_actions_out = actions_out.ContinuousActions;
        continuous_actions_out[0] = -Input.GetAxisRaw("Horizontal"); /* negative since the scene is reflected */
        continuous_actions_out[1] = Input.GetAxisRaw("Vertical");
        continuous_actions_out[2] = Input.GetAxisRaw("FrontBack");
    }
}