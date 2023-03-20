using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.MLAgents;
using Unity.MLAgents.Sensors;
using Unity.MLAgents.Actuators;

public class rov : Agent
{
    Rigidbody rigid_body;
    void Start()
    {
        rigid_body = GetComponent<Rigidbody>();
    }

    int overlap = 0;
    public bool is_overlap(){ return overlap > 0; }
    void OnTriggerEnter(Collider other){ ++overlap; }
    void OnTriggerExit(Collider other){ --overlap; }
    
    public Transform target;
    public Vector3 initial_position;
    public override void OnEpisodeBegin()
    {

    }

    public override void CollectObservations(VectorSensor sensor)
    {
        sensor.AddObservation(target.localPosition);
        sensor.AddObservation(this.transform.localPosition);
        sensor.AddObservation(rigid_body.velocity.x);
        sensor.AddObservation(rigid_body.velocity.y);
        sensor.AddObservation(rigid_body.velocity.z);
    }

    public float force_multiplier = 10;
    public override void OnActionReceived(ActionBuffers action_buffers)
    {
        Vector3 control_signal = Vector3.zero;
        control_signal.x = action_buffers.ContinuousActions[0];
        control_signal.y = action_buffers.ContinuousActions[1];
        control_signal.z = action_buffers.ContinuousActions[2];
        rigid_body.AddForce(force_multiplier * control_signal);

        float distance_to_target = Vector3.Distance(this.transform.localPosition, target.localPosition);
        if(distance_to_target < 1.42f)
        {
            SetReward(1.0f);
            EndEpisode();
        }
        if(is_overlap()) /* if the rov collides with an object */
        {
            SetReward(-1.0f);
        }
    }

    public float speed = 1f,
    sensitivity = 10f;
    public override void Heuristic(in ActionBuffers actions_out)
    {
        ActionSegment<float> continuous_actions_out = actions_out.ContinuousActions;
        if(Input.GetKey(KeyCode.LeftArrow) || Input.GetKey(KeyCode.A))
        {/* transform.Translate(new Vector3(-speed * Time.deltaTime, 0, 0)); */
            continuous_actions_out[0] -= speed * Time.deltaTime;
        }
        if(Input.GetKey(KeyCode.RightArrow) || Input.GetKey(KeyCode.D))
        {/* transform.Translate(new Vector3(speed * Time.deltaTime, 0, 0)); */
            continuous_actions_out[0] += speed * Time.deltaTime;
        }
        if(Input.GetKey(KeyCode.UpArrow) || Input.GetKey(KeyCode.W))
        {/* transform.Translate(new Vector3(0, speed * Time.deltaTime, 0)); */
            continuous_actions_out[1] += speed * Time.deltaTime;
        }
        if(Input.GetKey(KeyCode.DownArrow) || Input.GetKey(KeyCode.S))
        {/* transform.Translate(new Vector3(0, -speed * Time.deltaTime, 0)); */
            continuous_actions_out[1] -= speed * Time.deltaTime;
        }
        if(Input.GetKey(KeyCode.Alpha1) || Input.GetKey(KeyCode.Keypad1))
        {/* transform.Translate(new Vector3(0, 0, speed * Time.deltaTime)); */
            continuous_actions_out[2] += speed * Time.deltaTime;
        }
        if(Input.GetKey(KeyCode.Alpha2) || Input.GetKey(KeyCode.Keypad2))
        {/* transform.Translate(new Vector3(0, 0, -speed * Time.deltaTime)); */
            continuous_actions_out[2] -= speed * Time.deltaTime;
        }
    }
}
