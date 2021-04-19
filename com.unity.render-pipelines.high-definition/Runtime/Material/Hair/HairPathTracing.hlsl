#include "Packages/com.unity.render-pipelines.high-definition/Runtime/RenderPipeline/PathTracing/Shaders/PathTracingIntersection.hlsl"
#include "Packages/com.unity.render-pipelines.high-definition/Runtime/RenderPipeline/PathTracing/Shaders/PathTracingMaterial.hlsl"
#include "Packages/com.unity.render-pipelines.high-definition/Runtime/RenderPipeline/PathTracing/Shaders/PathTracingBSDF.hlsl"

// Disney Reference
// --------------------------------------------------------------------------------------

struct ReferenceInputs
{
    float thetaI;
    float thetaR;
    float thetaD;
    float thetaT;

    float phiI;
    float phiR;
    float phi;

    float LdotV;

    float h;

    float eta;
    float etaP;

    float shifts[3];
    float variances[3];
    float logisticScale;

    float3 absorption;
    float3 absorptionP;
};

float FresnelHair(float u)
{
    // TODO: SurfaceData.IOR (Exposed for animal fur, instead of hardcoded 1.55 for human hair).
    return F_Schlick(DEFAULT_HAIR_SPECULAR_VALUE, u);
}

float ModifiedIOR(float ior, float thetaD)
{
    float sinThetaD = sin(thetaD);
    float num = (ior * ior) - (sinThetaD * sinThetaD);
    return sqrt(num) / cos(thetaD);
}

float HyperbolicCosecant(float x)
{
    return rcp(sinh(x));
}

// Plot: https://www.desmos.com/calculator/4dnfmn9xal
float RoughnessToLongitudinalVariance(float roughness)
{
    float beta = roughness;
    float v = (0.726 * beta) + (0.812 * beta * beta) + (3.7 * pow(beta, 20.0));
    return v * v;
}

float RoughnessToLogisticalScale(float roughness)
{
    float beta = roughness;
    return (0.265 * beta) + (1.194 * beta * beta) + (5.372 * pow(beta, 22.0));
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

    UNITY_UNROLL
    for (int i = 0; i <= 10; ++i)
    {
        const float f = FACTORIAL[i];
        b += pow(x, 2.0 * i) / (pow(4, i) * f * f);
    }

    return b;
}

// Remap the azimuthal direction to the normalized logistic function on -PI to PI.
float RemapLogisticAngle(float a)
{
    if (a < -PI)
        a += TWO_PI;

    if (a > +PI)
        a -= TWO_PI;

    return a;
}

// Ref: Light Scattering from Human Hair Fibers
float AzimuthalDirection(uint p, float etaPrime, float h)
{
    float gammaI = asin(h);
    float gammaT = asin(h / etaPrime);
    float omega = (2 * p * gammaT) - (2 * gammaI) + p * PI;

    // Remap to the logistic function.
    return RemapLogisticAngle(omega);
}

float3 Attenuation(uint p, float h, float LdotV, float thetaD, float etaPrime, float3 absorption)
{
    float3 A;

    if (p == 0)
    {
        // Attenuation term for R is a special case.
        A = FresnelHair(0.5 * acos(LdotV));
    }
    else
    {
        float f = FresnelHair(acos(cos(thetaD) * cos(asin(h))));
        float gammaT = asin(h / etaPrime);
        float3 T = exp(-2 * absorption * (1 + cos(2 * gammaT)));

        // NOTE: This can't be used this due to NaNs via pow(f, 0)..
        // A = pow(1 - f, 2.0) * pow(f, p - 1) * pow(T, p);

        if (p == 1)
            A = pow(1 - f, 2.0) * T;
        else
            A = pow(1 - f, 2.0) * f * (T * T);
    }

    return A;
}

// Ref: [A Practical and Controllable Hair and Fur Model for Production Path Tracing]
// Plot: https://www.desmos.com/calculator/cmy0eig6ln
float LogisticAzimuthalAngularDistribution(float s, float phi)
{
    const float a = -PI;
    const float b = +PI;

    const float scalePeakTerm = sqrt(PI / 8.0);
    s *= scalePeakTerm;

    float normalizeTerm = rcp(rcp(1 + exp(a / s)) - rcp(1 + exp(b / s)));

    float distributionN = exp(-phi / s);

    float distributionD = 1 + exp(-phi / s);
    distributionD = s * distributionD * distributionD;

    return normalizeTerm * (distributionN / distributionD);
}

// Ref: [An Energy-Conserving Hair Reflectance Model]
// Plot: https://www.desmos.com/calculator/jmf1ofgfdv
float LongitudinalScattering(uint p, ReferenceInputs inputs)
{
    const float v      = max(0.0001, inputs.variances[p]);
    const float thetaI = inputs.thetaI;
    const float thetaR = inputs.thetaR - radians(inputs.shifts[p]);

    float M;

#if 1
    if (v < 0.1)
    {
        // Ref: [https://publons.com/review/414383/]
        // Small variances (< ~0.1) produce numerical issues due to limited floating precision.
        float a = (cos(-thetaI) * cos(thetaR)) / v;
        float b = (sin(-thetaI) * sin(thetaR)) / v;

        // The log of the bessel function may also be problematic for larger inputs (> ~12)...
        float lnI0;
        if (a > 12)
        {
            // ...in which case it's approximated.
            lnI0 = a + 0.5 * (-log(TWO_PI) + log(rcp(a)) + rcp(8 * a));
        }
        else
        {
            lnI0 = log(BesselI(a));
        }

        M = exp(lnI0 + b - rcp(v) + 0.6931 + log(rcp(2 * v)));
    }
    else
    {
        M  = HyperbolicCosecant(rcp(v)) / (2 * v);
        M *= exp((sin(-thetaI) * sin(thetaR)) / v);
        M *= BesselI((cos(-thetaI) * cos(thetaR)) / v);
    }
#else
    // TODO: Marschner Gaussian
    M = 1;
#endif

    return M;
}

float3 AzimuthalScattering(uint p, ReferenceInputs inputs)
{
    float3 N;

#if 1
    // Disney Integration of N(phi, h) (Near-Field).
    float3 A = Attenuation(p, inputs.h, inputs.LdotV, inputs.thetaD, inputs.etaP, inputs.absorptionP);

    float azimuth = AzimuthalDirection(p, inputs.etaP, inputs.h);
    float D = LogisticAzimuthalAngularDistribution(inputs.logisticScale, inputs.phi - azimuth);

    N = A * D;
#else
    // TODO: D'Eon's integration over fiber width and Gaussian Detector (Far-Field).
    N = 1;
#endif

    return N;
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

    // Construct a local frame with respect to strand and outgoing direction
    float3 U = mtlData.bsdfData.hairStrandDirectionWS;
    float3 V = normalize(cross(U, mtlData.V));
    float3 W = normalize(cross(U, V));

    // Transform to the local frame
    float3x3 frame = float3x3(W, V, U);
    float3 I = mul(frame, sampleDir);
    float3 R = mul(frame, mtlData.V);

    ReferenceInputs inputs;
    ZERO_INITIALIZE(ReferenceInputs, inputs);

    // Model Reference Inputs.
    // Notation Ref: Light Scattering from Human Hair Fibers
    {
        // Longitudinal
        inputs.thetaI = HALF_PI - acos(I.z);
        inputs.thetaR = HALF_PI - acos(R.z);
        inputs.thetaD = (inputs.thetaR - inputs.thetaI) * 0.5;

        // Azimuthal
        float phiI = atan2(I.y, I.x);
        float phiR = atan2(R.y, R.x);
        inputs.phi = phiR - phiI;

        // TODO: Move to ConvertSurfaceDataToBSDFData
        // Only the first lobe is attenuated by primary reflection roughness.
        float variance = RoughnessToLongitudinalVariance(mtlData.bsdfData.roughnessT);
        inputs.variances[0] = variance * mtlData.bsdfData.roughnessPrimaryReflection;
        inputs.variances[1] = 0.5 * variance;
        inputs.variances[2] = 2.0 * variance;

        inputs.shifts[0] = -mtlData.bsdfData.cuticleAngle;
        inputs.shifts[1] = -inputs.shifts[0] / 2.0;
        inputs.shifts[2] = 3 * -inputs.shifts[0] / 2.0;

        inputs.eta    = mtlData.bsdfData.ior;
        inputs.etaP   = ModifiedIOR(inputs.eta, inputs.thetaD);

        inputs.LdotV  = dot(sampleDir, mtlData.V);

        // Evaluation of h in the normal plane, given by gammaI = asin(h), where gammaI is the incident angle.
        // Since we are using a near-field method, we can use the true h value (rather than integrating over the whole fiber width).
        inputs.h = sin(acos(dot(cross(mtlData.bsdfData.normalWS, U), frame[1])));

        inputs.logisticScale = RoughnessToLogisticalScale(mtlData.bsdfData.roughnessB);

        float thetaT = asin(sin(inputs.thetaR / inputs.eta));
        inputs.absorptionP = mtlData.bsdfData.transmittance / cos(thetaT);
    }

    float3 S = 0;

    // Factored lobe representation. Sigma Sp(thetai, thetao, phi) = Mp(thetai, thetao) * Np(phi).
    for (uint p = 0; p < 3; p++)
    {
        // TEMP: Lobe selection
        // if (p == 0) continue;
        // if (p == 1) continue;
        // if (p == 2) continue;

        S += LongitudinalScattering(p, inputs) * AzimuthalScattering(p, inputs);
    }

    result.specValue = S;

    // TODO: Importance Sample
    result.specPdf = INV_FOUR_PI;
}

bool SampleMaterial(MaterialData mtlData, float3 inputSample, out float3 sampleDir, out MaterialResult result)
{
    Init(result);

    sampleDir = SampleSphereUniform(inputSample.x, inputSample.y);
    EvaluateMaterial(mtlData, sampleDir, result);

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
