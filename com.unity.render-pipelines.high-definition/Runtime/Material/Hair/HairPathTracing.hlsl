#include "Packages/com.unity.render-pipelines.high-definition/Runtime/RenderPipeline/PathTracing/Shaders/PathTracingIntersection.hlsl"
#include "Packages/com.unity.render-pipelines.high-definition/Runtime/RenderPipeline/PathTracing/Shaders/PathTracingMaterial.hlsl"
#include "Packages/com.unity.render-pipelines.high-definition/Runtime/RenderPipeline/PathTracing/Shaders/PathTracingBSDF.hlsl"

// Disney Reference
// --------------------------------------------------------------------------------------

float HyperbolicCosecant(float x)
{
    return 1 / sinh(x);
}

// Modified Bessel Function of the First Kind
float BesselI(float x)
{
    // Compute directly from the series for N = 10.
    const float tenFactSquared = 1.3168189e+13;

    float bessel = 0;

    for (int n = 0; n < 10; ++n)
    {
        bessel += pow(0.25 * x * x, 10.0) / tenFactSquared;
    }

    return bessel;
}

// Ref: An Energy-Conserving Hair Reflectance Model
// https://www.desmos.com/calculator/val4bnkeqy
float LongitudinalScattering(float variance, float thetaI, float thetaR)
{
    float a = HyperbolicCosecant(1.0 / variance) / 2 * variance;
    float b = exp((sin(-thetaI) * sin(thetaR)) / variance);
    float c = BesselI((cos(-thetaI) * cos(thetaR)) / variance);
    return a * b * c;
}

// Ref: An Energy-Conserving Hair Reflectance Model
// https://www.desmos.com/calculator/piy9dqbjot
void AzimuthalScattering(uint p)
{

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

    result.specValue = LongitudinalScattering(0.02, 0.5, 0.5);

    // BRDF::EvaluateLambert(mtlData, sampleDir, result.diffValue, result.diffPdf);
}

bool SampleMaterial(MaterialData mtlData, float3 inputSample, out float3 sampleDir, out MaterialResult result)
{
    Init(result);

    // TODO

    float3 value;
    float pdf;

    if (!BRDF::SampleLambert(mtlData, inputSample, sampleDir, value, pdf))
        return false;

    result.diffValue = value;
    result.diffPdf = pdf;

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
