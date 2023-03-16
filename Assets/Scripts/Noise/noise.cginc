float noise(float3 a)
{
    float r = 0.;
    for(int i = 0; i < 16; ++i)
    {
        float3 D, p = a + fmod(float3(i, i / 4, i / 8), float3(4.0, 2.0, 2.0)) +
            1.7 * sin(float3(i, 5 * i, 8 * i)), C = floor(p), P = p - C - .5, A = abs(P);
        C += fmod(C.x + C.y + C.z, 2.) * step(max(A.yzx, A.zxy), A) * sign(P);
        D = 34. * sin(987. * float(i) + 876. * C + 76. * C.yzx + 765. * C.zxy);
        P = p - C - .5;
        r += sin(6.3 * dot(P, frac(D) - .5)) * pow(max(0., 1. - 2. * dot(P, P)), 4.);
    }
    return .5 * sin(r);
}