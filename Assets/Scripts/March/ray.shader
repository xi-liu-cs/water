Shader "custom/ray"
{
    Properties
    {
        main_tex("texture", 2D) = "white" {}
    }
    SubShader
    {
        Tags {"Queue" = "Transparent" "RenderType" = "Transparent"} // "RenderType" = "Opaque"
        LOD 200

        Pass
        {
            CGPROGRAM
            #include "UnityCG.cginc"
            #include "../Noise/noise.cginc"
            #pragma vertex vert
            #pragma fragment frag
            #pragma multi_compile_fog
            #pragma Standard addshadow fullforwardshadows alpha:fade

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
                o.vertex = UnityObjectToClipPos(v.vertex); /* o.vertex = mul(UNITY_MATRIX_MVP, float4(pos, 1.0)); model * view * projection */
                o.uv = v.uv;
                o.ray_d = -ObjSpaceViewDir(v.vertex);
				o.ray_o = v.vertex.xyz - o.ray_d; /* o.ray_o = _WorldSpaceCameraPos; o.ray_d = v.vertex.xyz - o.ray_o; */
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

            float turbulence(float3 p)
            {
                float t = 0., f = 1.;
                for(int i = 0; i < 10; ++i)
                {
                    t += abs(noise(f * p)) / f;
                    f *= 2.;
                }
                return t;
            }

            #define DRAG_MULT 0.048
            #define ITERATIONS_RAYMARCH 13
            #define ITERATIONS_NORMAL 48
            /* #define Mouse (mouse.xy / resolution.xy) */
            #define Mouse (float2(0.5, 0.5)) /* if((0 < mouseX) && (mouseX < canvas_x) && (0 < mouseY) && (mouseY < canvas_y) custom_shader.setUniform("mouse", [mouseX, mouseY]);
            else custom_shader.setUniform("mouse", [canvas_x_over_2, canvas_y_over_2]); */
            #define Resolution (_ScreenParams.xy)
            #define Time _Time.y

            float2 wavedx(float2 position, float2 direction, float speed, float frequency, float timeshift)
            {
                float x = dot(direction, position) * frequency + timeshift * speed;
                float wave = exp(sin(x) - 1.0);
                float dx = wave * cos(x);
                return float2(wave, -dx);
            }

            float getwaves(float2 position, int iterations)
            {
                float iter = 0.0;
                float phase = 6.0;
                float speed = 2.0;
                float weight = 1.0;
                float w = 0.0;
                float ws = 0.0;
                for(int i = 0; i < iterations; i++)
                {
                    float2 p = float2(sin(iter), cos(iter));
                    float2 res = wavedx(position, p, speed, phase, Time);
                    position += p * res.y * weight * DRAG_MULT;
                    w += res.x * weight;
                    iter += 12.0;
                    ws += weight;
                    weight = lerp(weight, 0.0, 0.2);
                    phase *= 1.18;
                    speed *= 1.07;
                }
                return w / ws;
            }

            float raymarchwater(float3 camera, float3 start, float3 end, float depth)
            {
                float3 pos = start;
                float h = 0.0;
                float hupper = depth;
                float hlower = 0.0;
                float2 zer = float2(0.0, 0.0);
                float3 dir = normalize(end - start);
                float eps = 0.01;
                for(int i=0;i<318;i++)
                {
                    h = getwaves(pos.xz * 0.1, ITERATIONS_RAYMARCH) * depth - depth;
                    float dist_pos = distance(pos, camera);
                    if(h + eps*dist_pos > pos.y)
                    {
                        return dist_pos;
                    }
                    pos += dir * (pos.y - h);
                    //eps *= 1.01;
                }
                return -1.0;
            }

            float3 wave_normal(float2 pos, float e, float depth)
            {
                float2 ex = float2(e, 0);
                float H = getwaves(pos.xy * 0.1, ITERATIONS_NORMAL) * depth;
                float3 a = float3(pos.x, H, pos.y);
                return (cross(normalize(a-float3(pos.x - e, getwaves((pos.xy - ex.xy)*0.1, ITERATIONS_NORMAL) * depth, pos.y)), 
                                    normalize(a-float3(pos.x, getwaves((pos.xy + ex.yx )* 0.1, ITERATIONS_NORMAL) * depth, pos.y + e))));
            }
            float3x3 rotmat(float3 axis, float angle)
            {
                float s = sin(angle);
                float c = cos(angle);
                float oc = 1.0 - c;
                return float3x3(oc * axis.x * axis.x + c, oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s, 
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s, 
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c);
            }

            float3 getRay(float2 uv)
            {
                uv = (uv * 2.0 - 1.0) * float2(Resolution.x / Resolution.y, 1.0);
                float3 proj = normalize(float3(uv.x, uv.y, 1.0) + float3(uv.x, uv.y, -1.0) * pow(length(uv), 2.0) * 0.05);	
                if(Resolution.x < 400.0) return proj;
                float3 ray = mul(rotmat(float3(0.0, -1.0, 0.0), 3.0 * (Mouse.x * 2.0 - 1.0)) * rotmat(float3(1.0, 0.0, 0.0), 1.5 * (Mouse.y * 2.0 - 1.0)), proj);
                return ray;
            }

            float intersectPlane(float3 origin, float3 direction, float3 p, float3 normal)
            { 
                return clamp(dot(p - origin, normal) / dot(direction, normal), -1.0, 9991999.0); 
            }

            float3 extra_cheap_atmosphere(float3 raydir, float3 sundir)
            {
                sundir.y = max(sundir.y, -0.07);
                float special_trick = 1.0 / (raydir.y * 1.0 + 0.1);
                float special_trick2 = 1.0 / (sundir.y * 11.0 + 1.0);
                float raysundt = pow(abs(dot(sundir, raydir)), 2.0);
                float sundt = pow(max(0.0, dot(sundir, raydir)), 8.0);
                float mymie = sundt * special_trick * 0.2;
                float3 suncolor = lerp(float3(1.0, 1.0, 1.0), max(float3(0.0, 0.0, 0.0), float3(1.0, 1.0, 1.0) - float3(5.5, 13.0, 22.4) / 22.4), special_trick2);
                float3 bluesky = float3(5.5, 13.0, 22.4) / 22.4 * suncolor;
                float3 bluesky2 = max(float3(0.0, 0.0, 0.0), bluesky - float3(5.5, 13.0, 22.4) * 0.002 * (special_trick + -6.0 * sundir.y * sundir.y));
                bluesky2 *= special_trick * (0.24 + raysundt * 0.24);
                return bluesky2 * (1.0 + 1.0 * pow(1.0 - raydir.y, 3.0)) + mymie * suncolor;
            }

            float3 getatm(float3 ray)
            {
                return extra_cheap_atmosphere(ray, normalize(float3(1.0, 1.0, 1.0))) * 0.5;
            }

            float sun(float3 ray)
            {
                float3 sd = normalize(float3(1.0, 1.0, 1.0));   
                return pow(max(0.0, dot(ray, sd)), 528.0) * 110.0;
            }

            float3 aces_tonemap(float3 color)
            {	
                float3x3 m1 = float3x3(
                    0.59719, 0.07600, 0.02840,
                    0.35458, 0.90834, 0.13383,
                    0.04823, 0.01566, 0.83777
                );
                float3x3 m2 = float3x3(
                    1.60475, -0.10208, -0.00327,
                    -0.53108,  1.10813, -0.07276,
                    -0.07367, -0.00605,  1.07602
                );
                float3 v = mul(m1, color);    
                float3 a = v * (v + 0.0245786) - 0.000090537;
                float3 b = v * (0.983729 * v + 0.4329510) + 0.238081;
                float c = 1.0 / 2.2;
                return pow(clamp(mul(m2, (a / b)), 0.0, 1.0), float3(c, c, c));	
            }

            fixed4 frag(v2f i) : SV_Target
            {
                i.ray_d = normalize(i.ray_d);
                /* fixed4 col = tex2D(main_tex, i.uv); */
                fixed4 col = float4(0.1, 0.2, 0.6, 0.5);
                /* float3 light_direct = normalize(_WorldSpaceLightPos0.xyz); */
                float3 light_direct = normalize(float3(1., 1., 1.));
                float4 s = float4(0, 0, 0, 2);
                float t = ray_sphere(i.ray_o, i.ray_d, s);
                float n = turbulence(i.vertex.xyz);
                col += float4(sin(_Time.y * n), cos(_Time.y * n), sin(_Time.y * 1.3 * n), 0.) - n;
                if(t >= 0.)
                    col = float4(shade_sphere(at(i.ray_o, i.ray_d, t), s, light_direct), 1.);
                /* if(i.vertex.y < 0.5 * _ScreenParams.y) col = float4(0, 1, 0, 1); */
                float3 color_offset = float3(clamp(-cos(Time), -0.2, 0.2), 0., clamp(-sin(Time), -0.2, 0.2));
                float2 uv = i.uv; /* gl_FragCoord.xy / resolution.xy */
                float waterdepth = 2.1;
                float3 wfloor = float3(0.0, -waterdepth, 0.0);
                float3 wceil = float3(0.0, 0.0, 0.0);
                float3 orig = i.ray_o, /* float3 orig = float3(0.0, 2.0, 0.0); */
                ray = getRay(uv); /* ray = i.ray_d; */
                float hihit = intersectPlane(orig, ray, wceil, float3(0.0, 1.0, 0.0));
                if(ray.y >= -0.01){
                    float3 C = getatm(ray) * 2.0 + sun(ray);
                    //tonemapping
                    return float4(aces_tonemap(C) + color_offset, 1.);   
                }
                float lohit = intersectPlane(orig, ray, wfloor, float3(0.0, 1.0, 0.0));
                float3 hipos = orig + ray * hihit;
                float3 lopos = orig + ray * lohit;
                float dist = raymarchwater(orig, hipos, lopos, waterdepth);
                float3 pos = orig + ray * dist;

                float3 N = wave_normal(pos.xz, 0.001, waterdepth);
                float2 velocity = N.xz * (1.0 - N.y);
                N = lerp(float3(0.0, 1.0, 0.0), N, 1.0 / (dist * dist * 0.01 + 1.0));
                float3 R = reflect(ray, N);
                float fresnel = (0.04 + (1.0-0.04)*(pow(1.0 - max(0.0, dot(-N, ray)), 5.0)));
                
                float3 C = fresnel * getatm(R) * 2.0 + fresnel * sun(R);
                //tonemapping
                return float4(aces_tonemap(C) + color_offset, 1.);
            }
            ENDCG
        }
    }
}