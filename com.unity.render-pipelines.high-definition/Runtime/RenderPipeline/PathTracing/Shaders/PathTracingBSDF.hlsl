#ifndef UNITY_PATH_TRACING_BSDF_INCLUDED
#define UNITY_PATH_TRACING_BSDF_INCLUDED

#include "Packages/com.unity.render-pipelines.high-definition/Runtime/RenderPipeline/PathTracing/Shaders/PathTracingSampling.hlsl"
#include "Packages/com.unity.render-pipelines.high-definition/Runtime/RenderPipeline/Raytracing/Shaders/SubSurface.hlsl"

#define DELTA_PDF 1000000.0
#define MIN_GGX_ROUGHNESS 0.00001
#define MAX_GGX_ROUGHNESS 0.99999

float Lambda_AnisoGGX(float roughnessX,
                      float roughnessY,
                      float3 V)
{
    return 0.5 * (sqrt(1.0 + (Sq(roughnessX * V.x) + Sq(roughnessY * V.y)) / Sq(V.z)) - 1.0);
}

float G_AnisoGGX(float roughnessX,
                 float roughnessY,
                 float3 V)
{
    return rcp(1.0 + Lambda_AnisoGGX(roughnessX, roughnessY, V));
}

float D_AnisoGGX(float roughnessX,
                 float roughnessY,
                 float3 H)
{
    return rcp(PI * roughnessX * roughnessY * Sq(Sq(H.x / roughnessX) + Sq(H.y / roughnessY) + Sq(H.z)));
}

namespace BRDF
{

bool SampleGGX(MaterialData mtlData,
               float roughness,
               float3 fresnel0,
               float3 inputSample,
           out float3 outgoingDir,
           out float3 value,
           out float pdf,
           out float3 fresnel)
{
    roughness = clamp(roughness, MIN_GGX_ROUGHNESS, MAX_GGX_ROUGHNESS);

    float NdotL, NdotH, VdotH;
    float3x3 localToWorld = GetLocalFrame(GetSpecularNormal(mtlData));
    SampleGGXDir(inputSample.xy, mtlData.V, localToWorld, roughness, outgoingDir, NdotL, NdotH, VdotH);

    if (NdotL < 0.001 || !IsAbove(mtlData, outgoingDir))
        return false;

    float D = D_GGX(NdotH, roughness);
    pdf = D * NdotH / (4.0 * VdotH);

    if (pdf < 0.001)
        return false;

    float NdotV = dot(GetSpecularNormal(mtlData), mtlData.V);
    float Vg = V_SmithJointGGX(NdotL, NdotV, roughness);
    fresnel = F_Schlick(fresnel0, VdotH);

    value = fresnel * D * Vg * NdotL;

    return true;
}

bool SampleGGX(MaterialData mtlData,
               float roughness,
               float3 fresnel0,
               float3 inputSample,
           out float3 outgoingDir,
           out float3 value,
           out float pdf)
{
    float3 dummyFresnel;
    return SampleGGX(mtlData, roughness, fresnel0, inputSample, outgoingDir, value, pdf, dummyFresnel);
}

void EvaluateGGX(MaterialData mtlData,
                 float roughness,
                 float3 fresnel0,
                 float3 outgoingDir,
             out float3 value,
             out float pdf,
             out float3 fresnel)
{
    float NdotV = dot(GetSpecularNormal(mtlData), mtlData.V);
    if (NdotV < 0.001)
    {
        value = 0.0;
        pdf = 0.0;
        return;
    }
    float NdotL = dot(GetSpecularNormal(mtlData), outgoingDir);

    roughness = clamp(roughness, MIN_GGX_ROUGHNESS, MAX_GGX_ROUGHNESS);

    float3 H = normalize(mtlData.V + outgoingDir);
    float NdotH = dot(GetSpecularNormal(mtlData), H);
    float VdotH = dot(mtlData.V, H);
    float D = D_GGX(NdotH, roughness);
    pdf = D * NdotH / (4.0 * VdotH);

    float Vg = V_SmithJointGGX(NdotL, NdotV, roughness);
    fresnel = F_Schlick(fresnel0, VdotH);

    value = fresnel * D * Vg * NdotL;
}

void EvaluateGGX(MaterialData mtlData,
                float roughness,
                float3 fresnel0,
                float3 outgoingDir,
            out float3 value,
            out float pdf)
{
    float3 dummyFresnel;
    EvaluateGGX(mtlData, roughness, fresnel0, outgoingDir, value, pdf, dummyFresnel);
}

bool SampleAnisoGGX(MaterialData mtlData,
                    float3 fresnel0,
                    float3 inputSample,
                out float3 outgoingDir,
                out float3 value,
                out float pdf,
                out float3 fresnel)
{
    float roughnessX = clamp(mtlData.bsdfData.roughnessT, MIN_GGX_ROUGHNESS, MAX_GGX_ROUGHNESS);
    float roughnessY = clamp(mtlData.bsdfData.roughnessB, MIN_GGX_ROUGHNESS, MAX_GGX_ROUGHNESS);

    float VdotH;
    float3 localV, localH;
    float3x3 localToWorld = GetSpecularTangentFrame(mtlData);
    SampleAnisoGGXVisibleNormal(inputSample.xy, mtlData.V, localToWorld, roughnessX, roughnessY, localV, localH, VdotH);

    // Compute the reflection direction
    float3 localL = 2.0 * VdotH * localH - localV;
    outgoingDir = mul(localL, localToWorld);

    if (localL.z < 0.001 || !IsAbove(mtlData, outgoingDir))
        return false;

    float pdfNoGV = D_AnisoGGX(roughnessX, roughnessY, localH) / (4.0 * localV.z);
    float lambdaVPlusOne = Lambda_AnisoGGX(roughnessX, roughnessY, localV) + 1.0;
    pdf = pdfNoGV / lambdaVPlusOne;

    if (pdf < 0.001)
        return false;

    float lambdaL = Lambda_AnisoGGX(roughnessX, roughnessY, localL);
    fresnel = F_Schlick(fresnel0, VdotH);
    value = fresnel * pdfNoGV / (lambdaVPlusOne + lambdaL);

    return true;
}

bool SampleAnisoGGX(MaterialData mtlData,
                    float3 fresnel0,
                    float3 inputSample,
                out float3 outgoingDir,
                out float3 value,
                out float pdf)
{
    float3 dummyFresnel;
    return SampleAnisoGGX(mtlData, fresnel0, inputSample, outgoingDir, value, pdf, dummyFresnel);
}

void EvaluateAnisoGGX(MaterialData mtlData,
                      float3 fresnel0,
                      float3 outgoingDir,
                  out float3 value,
                  out float pdf,
                  out float3 fresnel)
{
    float NdotV = dot(GetSpecularNormal(mtlData), mtlData.V);
    if (NdotV < 0.001)
    {
        value = 0.0;
        pdf = 0.0;
        return;
    }

    float roughnessX = clamp(mtlData.bsdfData.roughnessT, MIN_GGX_ROUGHNESS, MAX_GGX_ROUGHNESS);
    float roughnessY = clamp(mtlData.bsdfData.roughnessB, MIN_GGX_ROUGHNESS, MAX_GGX_ROUGHNESS);

    float3x3 worldToLocal = transpose(GetSpecularTangentFrame(mtlData));
    float3 localV = mul(mtlData.V, worldToLocal);
    float3 localL = mul(outgoingDir, worldToLocal);
    float3 localH = normalize(localV + localL);
    float VdotH = dot(localV, localH);

    float pdfNoGV = D_AnisoGGX(roughnessX, roughnessY, localH) / (4.0 * localV.z);
    float lambdaVPlusOne = Lambda_AnisoGGX(roughnessX, roughnessY, localV) + 1.0;
    float lambdaL = Lambda_AnisoGGX(roughnessX, roughnessY, localL);

    fresnel = F_Schlick(fresnel0, VdotH);
    value = fresnel * pdfNoGV / (lambdaVPlusOne + lambdaL);
    pdf = pdfNoGV / lambdaVPlusOne;
}

void EvaluateAnisoGGX(MaterialData mtlData,
                      float3 fresnel0,
                      float3 outgoingDir,
                  out float3 value,
                  out float pdf)
{
    float3 dummyFresnel;
    EvaluateAnisoGGX(mtlData, fresnel0, outgoingDir, value, pdf, dummyFresnel);
}

bool SampleDelta(MaterialData mtlData,
             out float3 outgoingDir,
             out float3 value,
             out float pdf)
{
    if (IsAbove(mtlData))
    {
        outgoingDir = reflect(-mtlData.V, GetSpecularNormal(mtlData));
        float NdotV = dot(GetSpecularNormal(mtlData), mtlData.V);
        value = F_Schlick(mtlData.bsdfData.fresnel0, NdotV);
    }
    else // Below
    {
        outgoingDir = -reflect(mtlData.V, GetSpecularNormal(mtlData));
        float NdotV = -dot(GetSpecularNormal(mtlData), mtlData.V);
        value = F_FresnelDielectric(1.0 / mtlData.ior, NdotV);
    }

    value *= DELTA_PDF;
    pdf = DELTA_PDF;

    return any(outgoingDir);
}

bool SampleLambert(MaterialData mtlData,
                   float3 inputSample,
               out float3 outgoingDir,
               out float3 value,
               out float pdf)
{
    outgoingDir = SampleHemisphereCosine(inputSample.x, inputSample.y, GetDiffuseNormal(mtlData));

    if (!IsAbove(mtlData, outgoingDir))
        return false;

    pdf = dot(GetDiffuseNormal(mtlData), outgoingDir) * INV_PI;

    if (pdf < 0.001)
        return false;

    value = mtlData.bsdfData.diffuseColor * pdf;

    return true;
}

void EvaluateLambert(MaterialData mtlData,
                     float3 outgoingDir,
                 out float3 value,
                 out float pdf)
{
    pdf = saturate(dot(GetDiffuseNormal(mtlData), outgoingDir)) * INV_PI;
    value = mtlData.bsdfData.diffuseColor * pdf;
}

bool SampleBurley(MaterialData mtlData,
                  float3 inputSample,
              out float3 outgoingDir,
              out float3 value,
              out float pdf)
{
    outgoingDir = SampleHemisphereCosine(inputSample.x, inputSample.y, GetDiffuseNormal(mtlData));

    if (!IsAbove(mtlData, outgoingDir))
        return false;

    float NdotL = dot(GetDiffuseNormal(mtlData), outgoingDir);
    pdf = NdotL * INV_PI;

    if (pdf < 0.001)
        return false;

    float NdotV = saturate(dot(GetDiffuseNormal(mtlData), mtlData.V));
    float LdotV = saturate(dot(outgoingDir, mtlData.V));
    value = mtlData.bsdfData.diffuseColor * DisneyDiffuseNoPI(NdotV, NdotL, LdotV, mtlData.bsdfData.perceptualRoughness) * pdf;

    return true;
}

void EvaluateBurley(MaterialData mtlData,
                    float3 outgoingDir,
                out float3 value,
                out float pdf)
{
    float NdotL = saturate(dot(GetDiffuseNormal(mtlData), outgoingDir));
    float NdotV = saturate(dot(GetDiffuseNormal(mtlData), mtlData.V));
    float LdotV = saturate(dot(outgoingDir, mtlData.V));

    pdf = NdotL * INV_PI;
    value = mtlData.bsdfData.diffuseColor * DisneyDiffuseNoPI(NdotV, NdotL, LdotV, mtlData.bsdfData.perceptualRoughness) * pdf;
}

bool SampleDiffuse(MaterialData mtlData,
                   float3 inputSample,
               out float3 outgoingDir,
               out float3 value,
               out float pdf)
{
#ifdef USE_DIFFUSE_LAMBERT_BRDF
    return SampleLambert(mtlData, inputSample, outgoingDir, value, pdf);
#else
    return SampleBurley(mtlData, inputSample, outgoingDir, value, pdf);
#endif
}

void EvaluateDiffuse(MaterialData mtlData,
                     float3 outgoingDir,
                 out float3 value,
                 out float pdf)
{
#ifdef USE_DIFFUSE_LAMBERT_BRDF
    EvaluateLambert(mtlData, outgoingDir, value, pdf);
#else
    EvaluateBurley(mtlData, outgoingDir, value, pdf);
#endif
}

void EvaluateSheen(MaterialData mtlData,
                   float3 outgoingDir,
               out float3 value,
               out float pdf)
{
    // We use cosine-weighted sampling for this lobe
    float NdotL = saturate(dot(GetSpecularNormal(mtlData), outgoingDir));
    pdf = NdotL * INV_PI;

    float NdotV = dot(GetSpecularNormal(mtlData), mtlData.V);
    if (NdotV < 0.001)
    {
        value = 0.0;
        return;
    }

    float roughness = clamp(mtlData.bsdfData.roughnessT, MIN_GGX_ROUGHNESS, MAX_GGX_ROUGHNESS);

    float3 H = normalize(mtlData.V + outgoingDir);
    float NdotH = dot(GetSpecularNormal(mtlData), H);

    float D = D_Charlie(NdotH, roughness);

    // We use this visibility term to match the raster implementation (Fabric.hlsl)
    float Vg = V_Ashikhmin(NdotL, NdotV);
    //float Vg = V_Charlie(NdotL, NdotV, roughness);

    value = mtlData.bsdfData.fresnel0 * D * Vg * NdotL;
}

} // namespace BRDF

namespace BTDF
{

bool SampleGGX(MaterialData mtlData,
               float3 inputSample,
           out float3 outgoingDir,
           out float3 value,
           out float pdf)
{
    float roughness = clamp(mtlData.bsdfData.roughnessT, MIN_GGX_ROUGHNESS, MAX_GGX_ROUGHNESS);

    float NdotL, NdotH, VdotH;
    float3x3 localToWorld = GetLocalFrame(GetSpecularNormal(mtlData));
    SampleGGXDir(inputSample.xy, mtlData.V, localToWorld, roughness, outgoingDir, NdotL, NdotH, VdotH);

    // FIXME: won't be necessary after new version of SampleGGXDir()
    float3 H = normalize(mtlData.V + outgoingDir);
    outgoingDir = refract(-mtlData.V, H, 1.0 / mtlData.ior);
    NdotL = dot(GetSpecularNormal(mtlData), outgoingDir);

    if (NdotL > -0.001 || !IsBelow(mtlData, outgoingDir))
        return false;

    float NdotV = dot(GetSpecularNormal(mtlData), mtlData.V);
    float LdotH = dot(outgoingDir, H);

    float3 F = F_Schlick(mtlData.bsdfData.fresnel0, VdotH);
    float  D = D_GGX(NdotH, roughness);
    float Vg = V_SmithJointGGX(-NdotL, NdotV, roughness);

    // Compute the Jacobian
    float jacobian = max(abs(VdotH + mtlData.ior * LdotH), 0.001);
    jacobian = Sq(mtlData.ior) * abs(LdotH) / Sq(jacobian);

    pdf = D * NdotH * jacobian;
    value = abs(4.0 * (1.0 - F) * D * Vg * NdotL * VdotH * jacobian);

    return (pdf > 0.001);
}

bool SampleAnisoGGX(MaterialData mtlData,
                    float3 inputSample,
                out float3 outgoingDir,
                out float3 value,
                out float pdf)
{
    float roughnessX = clamp(mtlData.bsdfData.roughnessT, MIN_GGX_ROUGHNESS, MAX_GGX_ROUGHNESS);
    float roughnessY = clamp(mtlData.bsdfData.roughnessB, MIN_GGX_ROUGHNESS, MAX_GGX_ROUGHNESS);

    float VdotH;
    float3 localV, localH;
    float3x3 localToWorld = GetSpecularTangentFrame(mtlData);
    SampleAnisoGGXVisibleNormal(inputSample.xy, mtlData.V, localToWorld, roughnessX, roughnessY, localV, localH, VdotH);

    // Compute refraction direction instead of reflection
    float3 localL = refract(-localV, localH, 1.0 / mtlData.ior);
    outgoingDir = mul(localL, localToWorld);

    if (localL.z > -0.001 || !IsBelow(mtlData, outgoingDir))
        return false;

    // Compute the Jacobian
    float LdotH = dot(localL, localH);
    float jacobian = max(abs(VdotH + mtlData.ior * LdotH), 0.001);
    jacobian = Sq(mtlData.ior) * abs(LdotH) / Sq(jacobian);

    float3 F = F_Schlick(mtlData.bsdfData.fresnel0, VdotH);
    float  D = D_AnisoGGX(roughnessX, roughnessY, localH);

    float pdfNoGV = D * VdotH * jacobian / localV.z;
    float lambdaVPlusOne = Lambda_AnisoGGX(roughnessX, roughnessY, localV) + 1.0;
    float lambdaL = Lambda_AnisoGGX(roughnessX, roughnessY, localL);

    pdf = pdfNoGV / lambdaVPlusOne;
    value = abs((1.0 - F) * pdfNoGV / (lambdaVPlusOne + lambdaL));

    return (pdf > 0.001);
}

bool SampleDelta(MaterialData mtlData,
             out float3 outgoingDir,
             out float3 value,
             out float pdf)
{
    if (IsAbove(mtlData))
    {
        outgoingDir = refract(-mtlData.V, GetSpecularNormal(mtlData), 1.0 / mtlData.ior);
        float NdotV = dot(GetSpecularNormal(mtlData), mtlData.V);
        value = 1.0 - F_Schlick(mtlData.bsdfData.fresnel0, NdotV);
    }
    else // Below
    {
        outgoingDir = -refract(mtlData.V, GetSpecularNormal(mtlData), mtlData.ior);
        float NdotV = -dot(GetSpecularNormal(mtlData), mtlData.V);
        value = 1.0 - F_FresnelDielectric(1.0 / mtlData.ior, NdotV);
    }

    value *= DELTA_PDF;
    pdf = DELTA_PDF;

    return any(outgoingDir);
}

bool SampleLambert(MaterialData mtlData,
                   float3 inputSample,
               out float3 outgoingDir,
               out float3 value,
               out float pdf)
{
    bool retVal = BRDF::SampleLambert(mtlData, inputSample, outgoingDir, value, pdf);
    outgoingDir = -outgoingDir;
    return retVal;
}

void EvaluateLambert(MaterialData mtlData,
                     float3 outgoingDir,
                 out float3 value,
                 out float pdf)
{
    BRDF::EvaluateLambert(mtlData, -outgoingDir, value, pdf);
}

} // namespace BTDF

namespace SSS
{

#define MAX_WALK_STEPS 16
#define DIM_OFFSET 42

struct Result
{
    float3 throughput;
    float3 exitPosition;
    float3 exitNormal;
};

bool RandomWalk(float3 position, float3 normal, float3 diffuseColor, float3 meanFreePath, uint2 pixelCoord, out Result result)
{
    // Remap from our user-friendly parameters to and sigmaS and sigmaT
    float3 sigmaS, sigmaT;
    RemapSubSurfaceScatteringParameters(diffuseColor, meanFreePath, sigmaS, sigmaT);

    // Initialize the intersection structure
    PathIntersection intersection;
    intersection.remainingDepth = _RaytracingMaxRecursion + 1;

    // Initialize the walk parameters
    RayDesc rayDesc;
    rayDesc.Origin = position - normal * _RaytracingRayBias;
    rayDesc.TMin = 0.0;

    bool hit;
    uint walkIdx = 0;

    result.throughput = 1.0;

    do // Start our random walk
    {
        // Samples for direction, distance and channel selection
        float dirSample0 = GetSample(pixelCoord, _RaytracingSampleIndex, DIM_OFFSET + 4 * walkIdx + 0);
        float dirSample1 = GetSample(pixelCoord, _RaytracingSampleIndex, DIM_OFFSET + 4 * walkIdx + 1);
        float distSample = GetSample(pixelCoord, _RaytracingSampleIndex, DIM_OFFSET + 4 * walkIdx + 2);
        float channelSample = GetSample(pixelCoord, _RaytracingSampleIndex, DIM_OFFSET + 4 * walkIdx + 3);

        // Compute the per-channel weight
        float3 weights = result.throughput * SafeDivide(sigmaS, sigmaT);

        // Normalize our weights
        float wSum = weights.x + weights.y + weights.z;
        float3 channelWeights = SafeDivide(weights, wSum);

        // Evaluate what channel we should be using for this sample
        uint channelIdx = GetChannel(channelSample, channelWeights);

        // Evaluate the length of our steps
        rayDesc.TMax = -log(1.0 - distSample) / sigmaT[channelIdx];

        // Sample our next path segment direction
        rayDesc.Direction = walkIdx ?
            SampleSphereUniform(dirSample0, dirSample1) : SampleHemisphereCosine(dirSample0, dirSample1, -normal);

        // Initialize the intersection data
        intersection.t = -1.0;

        // Do the next step
        TraceRay(_RaytracingAccelerationStructure, RAY_FLAG_FORCE_OPAQUE | RAY_FLAG_CULL_FRONT_FACING_TRIANGLES,
                 RAYTRACINGRENDERERFLAG_PATH_TRACING, 0, 1, 1, rayDesc, intersection);

        // Check if we hit something
        hit = intersection.t > 0.0;

        // How much did the ray travel?
        float t = hit ? intersection.t : rayDesc.TMax;

        // Evaluate the transmittance for the current segment
        float3 transmittance = exp(-t * sigmaT);

        // Evaluate the pdf for the current segment
        float pdf = dot((hit ? transmittance : sigmaT * transmittance), channelWeights);

        // Contribute to the throughput
        result.throughput *= SafeDivide(hit ? transmittance : sigmaS * transmittance, pdf);

        // Compute the next path position
        rayDesc.Origin += rayDesc.Direction * t;

        // increment the path depth
        walkIdx++;
    }
    while (!hit && walkIdx < MAX_WALK_STEPS);

    // Set the exit intersection position and normal
    if (!hit)
    {
        result.exitPosition = position;
        result.exitNormal = normal;
        result.throughput = diffuseColor;

        // By not returning false here, we default to a diffuse BRDF when an intersection is not found;
        // this is physically wrong, but may prove more convenient for a user, as results will look
        // like diffuse instead of getting slightly darker when the mean free path becomes shorter.
        //return false;
    }
    else
    {
        result.exitPosition = rayDesc.Origin;
        result.exitNormal = intersection.value;
    }

    return true;
}

} // namespace SSS

#endif // UNITY_PATH_TRACING_BSDF_INCLUDED
