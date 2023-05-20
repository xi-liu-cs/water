Shader "Instanced/particle_opacity"
{
	Properties
	{
		_MainTex("Albedo (RGB)", 2D) = "white" {}
		_Glossiness("Smoothness", Range(0,1)) = 0.5
		_Metallic("Metallic", Range(0,1)) = 0.0
	}
	SubShader
	{
		Tags {"Queue" = "Transparent" "RenderType" = "Transparent" }
		LOD 200

		CGPROGRAM
		#include "../../Noise/noise.cginc"
		#pragma surface surf Standard addshadow fullforwardshadows alpha:fade
		#pragma multi_compile_instancing
		#pragma instancing_options procedural:setup

		sampler2D _MainTex;
		float size;

		struct Input
		{
			float2 uv_MainTex;
		};

		struct particle {
			float3 position;
			float speed;
			int isBoid;
			int touchedByBoid;
			float3 boidInfluence;
		};

		#ifdef UNITY_PROCEDURAL_INSTANCING_ENABLED
			StructuredBuffer<particle> particle_buffer;
		#endif

	void setup()
	{
	#ifdef UNITY_PROCEDURAL_INSTANCING_ENABLED
		float3 pos = particle_buffer[unity_InstanceID].position;
		unity_ObjectToWorld._11_21_31_41 = float4(size, 0, 0, 0);
		unity_ObjectToWorld._12_22_32_42 = float4(0, size, 0, 0);
		unity_ObjectToWorld._13_23_33_43 = float4(0, 0, size, 0);
		unity_ObjectToWorld._14_24_34_44 = float4(pos.xyz, 1);
		unity_WorldToObject = unity_ObjectToWorld;
		unity_WorldToObject._14_24_34 *= -1;
		unity_WorldToObject._11_22_33 = 1.0f / unity_WorldToObject._11_22_33;
	#endif
	}

	half _Glossiness;
	half _Metallic;

	void surf(Input IN, inout SurfaceOutputStandard o)
	{
		float r = 0.25;
		float g = 0.5;
		float b = 1.0;
		float a = 0.01;
		/*
		#ifdef UNITY_PROCEDURAL_INSTANCING_ENABLED
			//a = clamp(particle_buffer[unity_InstanceID].speed / 50.0, 0.02, 0.5);
			a = 0.01;
			if (particle_buffer[unity_InstanceID].touchedByBoid > 0) {
				r = 1;
				g = 0;
				b = 0;
				a = 1;
			}
		#endif
		*/
		float4 c = float4(r,g,b,a);
		/* #ifdef UNITY_PROCEDURAL_INSTANCING_ENABLED
		c = particle_buffer[unity_InstanceID].color;
		#endif */
		c *= tex2D(_MainTex, IN.uv_MainTex);
		o.Albedo = c.rgb + noise(c.rgb);
		o.Metallic = _Metallic;
		o.Smoothness = _Glossiness;
		o.Alpha = c.a;
	}
	ENDCG
	}
	FallBack "Diffuse"
}