#include "Packages/com.unity.render-pipelines.high-definition/Runtime/RenderPipeline/PathTracing/Shaders/PathTracingIntersection.hlsl"
#include "Packages/com.unity.render-pipelines.high-definition/Runtime/RenderPipeline/PathTracing/Shaders/PathTracingMaterial.hlsl"
#include "Packages/com.unity.render-pipelines.high-definition/Runtime/RenderPipeline/PathTracing/Shaders/PathTracingBSDF.hlsl"

// Disney Reference
// --------------------------------------------------------------------------------------

float RoughnessToLongitudinalVariance(float roughness)
{
    float beta = roughness;
    float v = (0.726 * beta) + (0.812 * beta * beta) * (3.7 * pow(beta, 20.0));
    return v * v;
}

float RoughnessToLogisticalScale(float roughness)
{
    float beta = roughness;
    return (0.265 * beta) + (1.194 * beta * beta) + (5.372 * pow(beta, 22.0));
}

float Attenuation()
{
    return 0;
}

float HyperbolicCosecant(float x)
{
    return rcp(sinh(x));
}

// Precompute the factorials (should really precompute the squared value).
static const float FACTORIAL[11] = { 1.0,
                                     1.0,
                                     2.0,
                                     6.0,
                                     24.0,
                                     120.0,
                                     720.0,
                                     5040.0,
                                     40320.0,
                                     362880.0,
                                     3628800.0 };

// Modified Bessel Function of the First Kind
float BesselI(float x)
{
    float b = 0;

    for (int i = 0; i <= 10; ++i)
    {
        const float f = FACTORIAL[i];
        b += pow(x, 2.0 * i) / (pow(4, i) * f * f);

        // TODO: Confirm this is correct + cheaper.
        // b += pow(0.25 * x * x, i) / (f * f);
    }

    return b;
}

float AzimuthalDirection(uint p, float etaPrime, float h)
{
    float gammaI = asin(h);
    float gammaT = asin(h / etaPrime);
    return (2 * p * gammaT) - (2 * gammaI) + (p * PI);
}

// Ref: A Practical and Controllable Hair and Fur Model for Production Path Tracing
// https://www.desmos.com/calculator/cmy0eig6ln
float LogisticAzimuthalAngularDistribution(float s, float phi)
{
    return 0;
}

// Ref: An Energy-Conserving Hair Reflectance Model
// https://www.desmos.com/calculator/qamyefnrki
float LongitudinalScattering(float variance, float thetaI, float thetaR)
{
    float a = HyperbolicCosecant(1.0 / variance) / (2 * variance);
    float b = exp((sin(-thetaI) * sin(thetaR)) / variance);
    float c = BesselI((cos(-thetaI) * cos(thetaR)) / variance);
    return a * b * c;
}

float AzimuthalScattering(uint p)
{
    return 0;
}

// --------------------------------------------------------------------------------------

void ProcessBSDFData(PathIntersection pathIntersection, BuiltinData builtinData, inout BSDFData bsdfData)
{
    // TODO
}

bool CreateMaterialData(PathIntersection pathIntersection, BuiltinData builtinData, BSDFData bsdfData, inout float3 shadingPosition, inout float theSample, out MaterialData mtlData)
{
    // Kajiya not supported.
    if (HasFlag(bsdfData.materialFeatures, MATERIALFEATUREFLAGS_HAIR_KAJIYA_KAY))
        return false;

    // Alter values in the material's bsdfData struct, to better suit path tracing
    mtlData.bsdfData = bsdfData;
    ProcessBSDFData(pathIntersection, builtinData, mtlData.bsdfData);

    mtlData.bsdfWeight = 0.0;
    mtlData.V = -WorldRayDirection();

    if (!IsAbove(mtlData))
        return false;

    return true;
}

void EvaluateMaterial(MaterialData mtlData, float3 sampleDir, out MaterialResult result)
{
    Init(result);

    float3 T = mtlData.bsdfData.hairStrandDirectionWS;
    float3 N = mtlData.bsdfData.normalWS;
    float3 V = mtlData.V;
    float3 L = sampleDir;

    float TdotL = dot(T, L);
    float TdotV = dot(T, V);
    float NdotL = dot(N, L);
    float NdotV = dot(N, V);

    // TEMP: Deal in pure angles for now, optimize later.
    float thetaI = asin(TdotL);
    float thetaR = asin(TdotV);
    float phiI   = acos(NdotL);
    float phiR   = acos(NdotV);
    float phi    = phiR - phiI;

    // TODO: Move to ConvertSurfaceDataToBSDFData
    float v = RoughnessToLongitudinalVariance(mtlData.bsdfData.roughnessT);
    float s = RoughnessToLogisticalScale(mtlData.bsdfData.roughnessB);

    const float shiftR  = mtlData.bsdfData.cuticleAngle;
    const float shiftTT = -shiftR / 2.0;
    const float shiftTRT = 3 * -shiftR / 2.0;

    result.specValue  = LongitudinalScattering(v, thetaI, thetaR - shiftR) * (0.25 * abs(cos(phi * 0.5))); // * AzimuthalScattering(0, s); // R
    // result.specValue += LongitudinalScattering(0.0, 0.0, 0.0) * AzimuthalScattering(); // TT
    // result.specValue += LongitudinalScattering(0.0, 0.0, 0.0) * AzimuthalScattering(); // TRT
    // result.specValue += LongitudinalScattering(0.0, 0.0, 0.0) * AzimuthalScattering(); // TRRT

    // TODO: PDF
}

bool SampleMaterial(MaterialData mtlData, float3 inputSample, out float3 sampleDir, out MaterialResult result)
{
    Init(result);

    // TODO

    return false;
}

float AdjustPathRoughness(MaterialData mtlData, MaterialResult mtlResult, bool isSampleBelow, float pathRoughness)
{
    // TODO

    return pathRoughness;
}

float3 ApplyAbsorption(MaterialData mtlData, SurfaceData surfaceData, float dist, bool isSampleBelow, float3 value)
{
    // TODO

    return value;
}
