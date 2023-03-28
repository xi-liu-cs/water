__global__ void mainImage(uchar4 * fragColor, float iTime)
{
	int width = 1024;
	int height = 1024;
	unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int i = x + width * y;
	float2 iResolution = make_float2((float)width, (float)height);
	float2 fragCoord = make_float2((float)x, (float)y);
	float2 uv = make_float2(fragCoord.x / iResolution.x, fragCoord.y / iResolution.y);
	float4 color = make_float4(uv.x, uv.y, 0.0, 1.0);
	fragColor[i] = make_uchar4(color.x * 255, color.y * 255, color.z * 255, 255);
}