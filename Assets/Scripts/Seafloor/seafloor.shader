Shader "Custom/seafloor"
{
    Properties
    {
        _Color ("Color", Color) = (1,1,1,1)
        _MainTex ("Albedo (RGB)", 2D) = "white" {}
        _Detail ("Detail", 2D) = "gray" {}
        _Glossiness ("Smoothness", Range(0,1)) = 0.5
        _Metallic ("Metallic", Range(0,1)) = 0.0
        _Tess("Tessellation", Range(1,16)) = 4
        _NoiseScale("Noise Scale", float) = 1
        _NoiseFrequency("Noise Frequency", float) = 1
        _NoiseOffset("Noise Offset", Vector) = (0, 0, 0, 0)
    }
    SubShader
    {
        Tags { "RenderType" = "Opaque" }
        LOD 200

        CGPROGRAM
        #include "../Noise/noise.cginc"
        // Physically based Standard lighting model, and enable shadows on all light types
        #pragma surface surf Standard fullforwardshadows tessellate:tess vertex:vert

        // Use shader model 3.0 target, to get nicer looking lighting
        #pragma target 4.6

        struct appdata
        {
            float4 vertex : POSITION;
            float3 normal : NORMAL;
            float4 tangent : TANGENT;
            float2 texcoord : TEXCOORD0;
        };

        sampler2D _MainTex,
        _Detail;

        struct Input
        {
            float2 uv_MainTex,
            uv_Detail;
        };

        half _Glossiness;
        half _Metallic;
        fixed4 _Color;
        float _Tess,
        _NoiseScale,
        _NoiseFrequency;
        float4 _NoiseOffset;
        float4 tess(){ return _Tess; }

        // Add instancing support for this shader. You need to check 'Enable Instancing' on materials that use the shader.
        // See https://docs.unity3d.com/Manual/GPUInstancing.html for more information about instancing.
        // #pragma instancing_options assumeuniformscaling
        UNITY_INSTANCING_BUFFER_START(Props)
            // put more per-instance properties here
        UNITY_INSTANCING_BUFFER_END(Props)

        void vert(inout appdata v)
        {
            float3 v1 = v.vertex.xyz,
            bitangent = cross(v.normal, v.tangent.xyz),
            v2 = v1 + 0.01 * v.tangent.xyz,
            v3 = v1 + 0.01 * bitangent;
            float ns1 = _NoiseScale * noise(_NoiseFrequency * float3(v1.x + _NoiseOffset.x, v1.y + _NoiseOffset.y, v1.z + _NoiseOffset.z));
            v1.xyz += ((ns1 + 1) / 2) * v.normal;
            float ns2 = _NoiseScale * noise(_NoiseFrequency * float3(v2.x + _NoiseOffset.x, v2.y + _NoiseOffset.y, v2.z + _NoiseOffset.z));
            v2.xyz += ((ns2 + 1) / 2) * v.normal;
            float ns3 = _NoiseScale * noise(_NoiseFrequency * float3(v3.x + _NoiseOffset.x, v3.y + _NoiseOffset.y, v3.z + _NoiseOffset.z));
            v3.xyz += ((ns2 + 1) / 2) * v.normal;
            v.normal = normalize(cross(v2 - v1, v3 - v1));
            v.vertex.xyz = v1;
        }

        void surf(Input IN, inout SurfaceOutputStandard o)
        {
            // Albedo comes from a texture tinted by color
            fixed4 c = tex2D(_MainTex, IN.uv_MainTex) * _Color;
            o.Albedo = c.rgb;
            o.Albedo *= tex2D(_Detail, IN.uv_Detail).rgb * 2;
            // Metallic and smoothness come from slider variables
            o.Metallic = _Metallic;
            o.Smoothness = _Glossiness;
            o.Alpha = c.a;
        }
        ENDCG
    }
    FallBack "Diffuse"
}
