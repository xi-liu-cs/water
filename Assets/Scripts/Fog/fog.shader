Shader "custom/fog"
{
	Properties
	{
		[HDR]_Albedo("Albedo", Color) = (1,1,1,1)
		_SimpleNoiseScale("Simple Noise Scale", Float) = 20
		_SimplexNoiseScale("Simplex Noise Scale", Float) = 4
		_VoronoiScale("Voronoi Scale", Float) = 5
		_SimpleNoiseAnimation("Simple Noise Animation", Vector) = (0,0,0,0)
		_SimplexNoiseAnimation("Simplex Noise Animation", Vector) = (0,0,0.02,0)
		_VoronoiNoiseAnimation("Voronoi Noise Animation", Vector) = (0,0,0,0)
		_SimpleNoiseAmount("Simple Noise Amount", Range( 0 , 1)) = 0.25
		_SimplexNoiseAmount("Simplex Noise Amount", Range( 0 , 1)) = 0.25
		_VoronoiNoiseAmount("Voronoi Noise Amount", Range( 0 , 1)) = 0.5
		_SimpleNoiseRemap("Simple Noise Remap", Range( 0 , 1)) = 0
		_SimplexNoiseRemap("Simplex Noise Remap", Range( 0 , 1)) = 0
		_VoronoiNoiseRemap("Voronoi Noise Remap", Range( 0 , 1)) = 0
		_CombinedNoiseRemap("Combined Noise Remap", Range( 0 , 1)) = 0
		_SurfaceDepthFade("Surface Depth Fade", Float) = 0
		_CameraDepthFadeRange("Camera Depth Fade Range", Float) = 0
		_CameraDepthFadeOffset("Camera Depth Fade Offset", Float) = 0
		[HideInInspector] _texcoord( "", 2D ) = "white" {}
		[HideInInspector] __dirty( "", Int ) = 1
	}

	SubShader
	{
		Tags{ "RenderType" = "Transparent"  "Queue" = "Transparent+0" "IgnoreProjector" = "True" "IsEmissive" = "true"  }
		Cull Back
		CGINCLUDE
		#include "UnityShaderVariables.cginc"
		#include "UnityCG.cginc"
		#include "UnityPBSLighting.cginc"
		#include "Lighting.cginc"
		#pragma target 3.5
		#undef TRANSFORM_TEX
		#define TRANSFORM_TEX(tex,name) float4(tex.xy * name##_ST.xy + name##_ST.zw, tex.z, tex.w)
		struct Input
		{
			float4 vertexColor : COLOR;
			float4 uv_texcoord;
			float4 screenPos;
			float eyeDepth;
		};

		uniform float4 _Albedo;
		uniform float2 _SimpleNoiseAnimation;
		uniform float _SimpleNoiseScale;
		uniform float _SimpleNoiseRemap;
		uniform float _SimpleNoiseAmount;
		uniform float3 _SimplexNoiseAnimation;
		uniform float _SimplexNoiseScale;
		uniform float _SimplexNoiseRemap;
		uniform float _SimplexNoiseAmount;
		uniform float _VoronoiScale;
		uniform float3 _VoronoiNoiseAnimation;
		uniform float _VoronoiNoiseRemap;
		uniform float _VoronoiNoiseAmount;
		uniform float _CombinedNoiseRemap;
		UNITY_DECLARE_DEPTH_TEXTURE( _CameraDepthTexture );
		uniform float4 _CameraDepthTexture_TexelSize;
		uniform float _SurfaceDepthFade;
		uniform float _CameraDepthFadeRange;
		uniform float _CameraDepthFadeOffset;


		inline float noise_randomValue (float2 uv) { return frac(sin(dot(uv, float2(12.9898, 78.233)))*43758.5453); }

		inline float noise_interpolate (float a, float b, float t) { return (1.0-t)*a + (t*b); }

		inline float valueNoise (float2 uv)
		{
			float2 i = floor(uv);
			float2 f = frac( uv );
			f = f* f * (3.0 - 2.0 * f);
			uv = abs( frac(uv) - 0.5);
			float2 c0 = i + float2( 0.0, 0.0 );
			float2 c1 = i + float2( 1.0, 0.0 );
			float2 c2 = i + float2( 0.0, 1.0 );
			float2 c3 = i + float2( 1.0, 1.0 );
			float r0 = noise_randomValue( c0 );
			float r1 = noise_randomValue( c1 );
			float r2 = noise_randomValue( c2 );
			float r3 = noise_randomValue( c3 );
			float bottomOfGrid = noise_interpolate( r0, r1, f.x );
			float topOfGrid = noise_interpolate( r2, r3, f.x );
			float t = noise_interpolate( bottomOfGrid, topOfGrid, f.y );
			return t;
		}


		float SimpleNoise(float2 UV)
		{
			float t = 0.0;
			float freq = pow( 2.0, float( 0 ) );
			float amp = pow( 0.5, float( 3 - 0 ) );
			t += valueNoise( UV/freq )*amp;
			freq = pow(2.0, float(1));
			amp = pow(0.5, float(3-1));
			t += valueNoise( UV/freq )*amp;
			freq = pow(2.0, float(2));
			amp = pow(0.5, float(3-2));
			t += valueNoise( UV/freq )*amp;
			return t;
		}


		float3 mod3D289( float3 x ) { return x - floor( x / 289.0 ) * 289.0; }

		float4 mod3D289( float4 x ) { return x - floor( x / 289.0 ) * 289.0; }

		float4 permute( float4 x ) { return mod3D289( ( x * 34.0 + 1.0 ) * x ); }

		float4 taylorInvSqrt( float4 r ) { return 1.79284291400159 - r * 0.85373472095314; }

		float snoise( float3 v )
		{
			const float2 C = float2( 1.0 / 6.0, 1.0 / 3.0 );
			float3 i = floor( v + dot( v, C.yyy ) );
			float3 x0 = v - i + dot( i, C.xxx );
			float3 g = step( x0.yzx, x0.xyz );
			float3 l = 1.0 - g;
			float3 i1 = min( g.xyz, l.zxy );
			float3 i2 = max( g.xyz, l.zxy );
			float3 x1 = x0 - i1 + C.xxx;
			float3 x2 = x0 - i2 + C.yyy;
			float3 x3 = x0 - 0.5;
			i = mod3D289( i);
			float4 p = permute( permute( permute( i.z + float4( 0.0, i1.z, i2.z, 1.0 ) ) + i.y + float4( 0.0, i1.y, i2.y, 1.0 ) ) + i.x + float4( 0.0, i1.x, i2.x, 1.0 ) );
			float4 j = p - 49.0 * floor( p / 49.0 );  // mod(p,7*7)
			float4 x_ = floor( j / 7.0 );
			float4 y_ = floor( j - 7.0 * x_ );  // mod(j,N)
			float4 x = ( x_ * 2.0 + 0.5 ) / 7.0 - 1.0;
			float4 y = ( y_ * 2.0 + 0.5 ) / 7.0 - 1.0;
			float4 h = 1.0 - abs( x ) - abs( y );
			float4 b0 = float4( x.xy, y.xy );
			float4 b1 = float4( x.zw, y.zw );
			float4 s0 = floor( b0 ) * 2.0 + 1.0;
			float4 s1 = floor( b1 ) * 2.0 + 1.0;
			float4 sh = -step( h, 0.0 );
			float4 a0 = b0.xzyw + s0.xzyw * sh.xxyy;
			float4 a1 = b1.xzyw + s1.xzyw * sh.zzww;
			float3 g0 = float3( a0.xy, h.x );
			float3 g1 = float3( a0.zw, h.y );
			float3 g2 = float3( a1.xy, h.z );
			float3 g3 = float3( a1.zw, h.w );
			float4 norm = taylorInvSqrt( float4( dot( g0, g0 ), dot( g1, g1 ), dot( g2, g2 ), dot( g3, g3 ) ) );
			g0 *= norm.x;
			g1 *= norm.y;
			g2 *= norm.z;
			g3 *= norm.w;
			float4 m = max( 0.6 - float4( dot( x0, x0 ), dot( x1, x1 ), dot( x2, x2 ), dot( x3, x3 ) ), 0.0 );
			m = m* m;
			m = m* m;
			float4 px = float4( dot( x0, g0 ), dot( x1, g1 ), dot( x2, g2 ), dot( x3, g3 ) );
			return 42.0 * dot( m, px);
		}


		float2 voronoihash2( float2 p )
		{
			
			p = float2( dot( p, float2( 127.1, 311.7 ) ), dot( p, float2( 269.5, 183.3 ) ) );
			return frac( sin( p ) *43758.5453);
		}


		float voronoi2( float2 v, float time, inout float2 id, inout float2 mr, float smoothness, inout float2 smoothId )
		{
			float2 n = floor( v );
			float2 f = frac( v );
			float F1 = 8.0;
			float F2 = 8.0; float2 mg = 0;
			for ( int j = -1; j <= 1; j++ )
			{
				for ( int i = -1; i <= 1; i++ )
			 	{
			 		float2 g = float2( i, j );
			 		float2 o = voronoihash2( n + g );
					o = ( sin( time + o * 6.2831 ) * 0.5 + 0.5 ); float2 r = f - g - o;
					float d = 0.5 * dot( r, r );
			 		if( d<F1 ) {
			 			F2 = F1;
			 			F1 = d; mg = g; mr = r; id = o;
			 		} else if( d<F2 ) {
			 			F2 = d;
			
			 		}
			 	}
			}
			return (F2 + F1) * 0.5;
		}


		void vertexDataFunc( inout appdata_full v, out Input o )
		{
			UNITY_INITIALIZE_OUTPUT( Input, o );
			o.eyeDepth = -UnityObjectToViewPos( v.vertex.xyz ).z;
		}

		inline half4 LightingUnlit( SurfaceOutput s, half3 lightDir, half atten )
		{
			return half4 ( 0, 0, 0, s.Alpha );
		}

		void surf( Input i , inout SurfaceOutput o )
		{
			float4 Albedo81 = ( _Albedo * i.vertexColor );
			o.Emission = Albedo81.rgb;
			float particleStableRandom43 = i.uv_texcoord.z;
			float2 uvs_TexCoord3 = i.uv_texcoord;
			uvs_TexCoord3.xy = i.uv_texcoord.xy + ( ( _SimpleNoiseAnimation * _Time.y ) + ( particleStableRandom43 * 10.0 ) );
			float simpleNoise1 = SimpleNoise( uvs_TexCoord3.xy*_SimpleNoiseScale );
			float SimpleNoise18 = saturate( (0.0 + (simpleNoise1 - _SimpleNoiseRemap) * (1.0 - 0.0) / (1.0 - _SimpleNoiseRemap)) );
			float lerpResult34 = lerp( 1.0 , SimpleNoise18 , _SimpleNoiseAmount);
			float simplePerlin3D19 = snoise( ( float3( i.uv_texcoord.xy ,  0.0 ) + ( _SimplexNoiseAnimation * _Time.y ) + ( particleStableRandom43 * 20.0 ) )*_SimplexNoiseScale );
			simplePerlin3D19 = simplePerlin3D19*0.5 + 0.5;
			float SimplexNoise25 = saturate( (0.0 + (simplePerlin3D19 - _SimplexNoiseRemap) * (1.0 - 0.0) / (1.0 - _SimplexNoiseRemap)) );
			float lerpResult33 = lerp( 1.0 , SimplexNoise25 , _SimplexNoiseAmount);
			float mulTime6 = _Time.y * _VoronoiNoiseAnimation.z;
			float time2 = mulTime6;
			float2 voronoiSmoothId2 = 0;
			float2 uvs_TexCoord7 = i.uv_texcoord;
			uvs_TexCoord7.xy = i.uv_texcoord.xy + ( (_VoronoiNoiseAnimation).xy * _Time.y );
			float2 coords2 = uvs_TexCoord7.xy * _VoronoiScale;
			float2 id2 = 0;
			float2 uv2 = 0;
			float voroi2 = voronoi2( coords2, time2, id2, uv2, 0, voronoiSmoothId2 );
			float VoronoiNoise12 = saturate( (0.0 + (voroi2 - _VoronoiNoiseRemap) * (1.0 - 0.0) / (1.0 - _VoronoiNoiseRemap)) );
			float lerpResult5 = lerp( 1.0 , VoronoiNoise12 , _VoronoiNoiseAmount);
			float Noise36 = ( lerpResult34 * lerpResult33 * lerpResult5 );
			float RemappedNoise60 = saturate( (0.0 + (Noise36 - _CombinedNoiseRemap) * (1.0 - 0.0) / (1.0 - _CombinedNoiseRemap)) );
			float4 ase_screenPos = float4( i.screenPos.xyz , i.screenPos.w + 0.00000000001 );
			float4 ase_screenPosNorm = ase_screenPos / ase_screenPos.w;
			ase_screenPosNorm.z = ( UNITY_NEAR_CLIP_VALUE >= 0 ) ? ase_screenPosNorm.z : ase_screenPosNorm.z * 0.5 + 0.5;
			float screenDepth50 = LinearEyeDepth(SAMPLE_DEPTH_TEXTURE( _CameraDepthTexture, ase_screenPosNorm.xy ));
			float distanceDepth50 = abs( ( screenDepth50 - LinearEyeDepth( ase_screenPosNorm.z ) ) / ( _SurfaceDepthFade ) );
			float SurfaceDepthFade55 = saturate( distanceDepth50 );
			float cameraDepthFade87 = (( i.eyeDepth -_ProjectionParams.y - _CameraDepthFadeOffset ) / _CameraDepthFadeRange);
			float CameraDepthFade88 = saturate( cameraDepthFade87 );
			float2 uvs_TexCoord76 = i.uv_texcoord;
			uvs_TexCoord76.xy = i.uv_texcoord.xy * float2( 2,2 ) + float2( -1,-1 );
			float RadialMask79 = saturate( ( 1.0 - length( uvs_TexCoord76.xy ) ) );
			o.Alpha = ( RemappedNoise60 * SurfaceDepthFade55 * CameraDepthFade88 * RadialMask79 * Albedo81.a );
		}

		ENDCG
		CGPROGRAM
		#pragma surface surf Unlit alpha:fade keepalpha fullforwardshadows vertex:vertexDataFunc 

		ENDCG
		Pass
		{
			Name "ShadowCaster"
			Tags{ "LightMode" = "ShadowCaster" }
			ZWrite On
			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			#pragma target 3.5
			#pragma multi_compile_shadowcaster
			#pragma multi_compile UNITY_PASS_SHADOWCASTER
			#pragma skip_variants FOG_LINEAR FOG_EXP FOG_EXP2
			#include "HLSLSupport.cginc"
			#if ( SHADER_API_D3D11 || SHADER_API_GLCORE || SHADER_API_GLES || SHADER_API_GLES3 || SHADER_API_METAL || SHADER_API_VULKAN )
				#define CAN_SKIP_VPOS
			#endif
			#include "UnityCG.cginc"
			#include "Lighting.cginc"
			#include "UnityPBSLighting.cginc"
			sampler3D _DitherMaskLOD;
			struct v2f
			{
				V2F_SHADOW_CASTER;
				float4 customPack1 : TEXCOORD1;
				float1 customPack2 : TEXCOORD2;
				float3 worldPos : TEXCOORD3;
				float4 screenPos : TEXCOORD4;
				half4 color : COLOR0;
				UNITY_VERTEX_INPUT_INSTANCE_ID
				UNITY_VERTEX_OUTPUT_STEREO
			};
			v2f vert( appdata_full v )
			{
				v2f o;
				UNITY_SETUP_INSTANCE_ID( v );
				UNITY_INITIALIZE_OUTPUT( v2f, o );
				UNITY_INITIALIZE_VERTEX_OUTPUT_STEREO( o );
				UNITY_TRANSFER_INSTANCE_ID( v, o );
				Input customInputData;
				vertexDataFunc( v, customInputData );
				float3 worldPos = mul( unity_ObjectToWorld, v.vertex ).xyz;
				half3 worldNormal = UnityObjectToWorldNormal( v.normal );
				o.customPack1.xyzw = customInputData.uv_texcoord;
				o.customPack1.xyzw = v.texcoord;
				o.customPack2.x = customInputData.eyeDepth;
				o.worldPos = worldPos;
				TRANSFER_SHADOW_CASTER_NORMALOFFSET( o )
				o.screenPos = ComputeScreenPos( o.pos );
				o.color = v.color;
				return o;
			}
			half4 frag( v2f IN
			#if !defined( CAN_SKIP_VPOS )
			, UNITY_VPOS_TYPE vpos : VPOS
			#endif
			) : SV_Target
			{
				UNITY_SETUP_INSTANCE_ID( IN );
				Input surfIN;
				UNITY_INITIALIZE_OUTPUT( Input, surfIN );
				surfIN.uv_texcoord = IN.customPack1.xyzw;
				surfIN.eyeDepth = IN.customPack2.x;
				float3 worldPos = IN.worldPos;
				half3 worldViewDir = normalize( UnityWorldSpaceViewDir( worldPos ) );
				surfIN.screenPos = IN.screenPos;
				surfIN.vertexColor = IN.color;
				SurfaceOutput o;
				UNITY_INITIALIZE_OUTPUT( SurfaceOutput, o )
				surf( surfIN, o );
				#if defined( CAN_SKIP_VPOS )
				float2 vpos = IN.pos;
				#endif
				half alphaRef = tex3D( _DitherMaskLOD, float3( vpos.xy * 0.25, o.Alpha * 0.9375 ) ).a;
				clip( alphaRef - 0.01 );
				SHADOW_CASTER_FRAGMENT( IN )
			}
			ENDCG
		}
	}
	Fallback "Diffuse"
	CustomEditor "ASEMaterialInspector"
}