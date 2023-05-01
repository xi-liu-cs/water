Shader "custom/ray"
{
    Properties
    {
        main_tex("texture", 2D) = "white" {}
    }
    SubShader
    {
        Tags {"RenderType" = "Opaque"}
        LOD 100

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma multi_compile_fog
            #include "UnityCG.cginc"

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                float3 ray_o : TEXCOORD1, /* origin */
			    ray_d : TEXCOORD2; /* direction */
                float4 vertex : SV_POSITION;
            };

            sampler2D main_tex;
            float4 main_tex_ST;

            v2f vert(appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = v.uv;
                o.ray_o = _WorldSpaceCameraPos;
                o.ray_d = v.vertex.xyz - o.ray_o; /* o.ray_d = -ObjSpaceViewDir(v.vertex); o.ray_o = v.vertex.xyz - o.ray_d; */
                return o;
            }

            #define at(origin, direct, t)(origin + t * direct)

            /* (x - c_x) ^ 2 + (y - c_y) ^ 2 + (z - c_z) ^ 2 = r ^ 2
            center c = (c_x, c_y, c_z)
            point p = (x, y, z)
            p = c + r
            p - c = r
            (p - c) ^ 2 = r ^ 2
            ray p(t) = o + td
            (o + td - c) ^ 2 = r ^ 2
            (td + (o - c)) ^ 2 - r ^ 2 = 0
            (d ^ 2)t ^ 2 + (2d(o - c))t + (o - c) ^ 2 - r ^ 2 = 0
            quadratic
            a = d ^ 2
            b = 2d(o - c)
            c = (o - c) ^ 2 - r ^ 2
            float ray_sphere(float3 origin, float3 direct, float4 s)
            {
                origin -= s.xyz;
                float a = dot(direct, direct),
                b = 2. * dot(origin, direct),
                c = dot(origin, origin) - s.w * s.w,
                d = b * b - 4 * a * c;
                return d < 0. ? -1. : (-b - sqrt(d)) / (2. * a); 
            }
            b = 2h
            (-b +- sqrt{b ^ 2 - 4ac}) / (2a)
            = (-2h +- sqrt{(2h) ^ 2 - 4ac}) / (2a)
            = (-2h +- 2sqrt{h ^ 2 - ac}) / (2a)
            = (-h +- sqrt{h ^ 2 - ac}) / a */
            float ray_sphere(float3 origin, float3 direct, float4 s)
            {
                origin -= s.xyz;
                float a = dot(direct, direct),
                b = dot(origin, direct),
                c = dot(origin, origin) - s.w * s.w,
                d = b * b - a * c;
                return d < 0. ? -1. : (-b - sqrt(d)) / a; 
            }

            float3 shade_sphere(float3 p, float4 s, float3 light_direct)
            {
                float3 n = normalize(p - s.xyz),
                c = float3(1., 0., 0.);
                c += max(0., dot(n, light_direct));
                return c;
            }

            fixed4 frag(v2f i) : SV_Target
            {
                i.ray_d = normalize(i.ray_d);
                /* fixed4 col = tex2D(main_tex, i.uv); */
                fixed4 col = float4(0.1, 0.2, 0.6, 0.5);
                /* float3 light_direct = normalize(_WorldSpaceLightPos0.xyz); */
                float3 light_direct = normalize(float3(1., 1., 1.));
                float4 s = float4(0, 0, 0, .3);
                float t = ray_sphere(i.ray_o, i.ray_d, s);
                if(t > 0.)
                    col = float4(shade_sphere(at(i.ray_o, i.ray_d, t), s, light_direct), 1.);
                return col;
            }
            ENDCG
        }
    }
}