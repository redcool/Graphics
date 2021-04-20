using System;
using System.Linq;
using UnityEditor.ShaderGraph;

namespace UnityEditor.Rendering.Universal.ShaderGraph
{
    internal class UniversalSpriteSubTargetShared
    {
        #region Passes
        internal static class SpriteSharedPasses
        {
            public static readonly PassDescriptor SceneSelection = new PassDescriptor()
            {
                // Definition
                displayName = "SceneSelectionPass",
                referenceName = "SHADERPASS_DEPTHONLY",
                lightMode = "SceneSelectionPass",
                useInPreview = true,

                // Template
                passTemplatePath = GenerationUtils.GetDefaultTemplatePath("PassMesh.template"),
                sharedTemplateDirectories = UniversalTarget.kSharedTemplateDirectories,

                // Port Mask
                validVertexBlocks = CoreBlockMasks.Vertex,
                validPixelBlocks = CoreBlockMasks.FragmentAlphaOnly,

                // Fields
                structs = CoreStructCollections.Default,
                fieldDependencies = CoreFieldDependencies.Default,

                // Conditional State
                renderStates = CoreRenderStates.SceneSelection,
                pragmas = CorePragmas.Instanced,
                defines = CoreDefines.SceneSelection,
                includes = SpriteSharedIncludes.DepthOnly,

                // Custom Interpolator Support
                customInterpolators = CoreCustomInterpDescriptors.Common
            };

            public static readonly PassDescriptor ScenePicking = new PassDescriptor()
            {
                // Definition
                displayName = "ScenePickingPass",
                referenceName = "SHADERPASS_DEPTHONLY",
                lightMode = "Picking",
                useInPreview = true,

                // Template
                passTemplatePath = GenerationUtils.GetDefaultTemplatePath("PassMesh.template"),
                sharedTemplateDirectories = UniversalTarget.kSharedTemplateDirectories,

                // Port Mask
                validVertexBlocks = CoreBlockMasks.Vertex,
                validPixelBlocks = CoreBlockMasks.FragmentAlphaOnly,

                // Fields
                structs = CoreStructCollections.Default,
                fieldDependencies = CoreFieldDependencies.Default,

                // Conditional State
                renderStates = CoreRenderStates.ScenePicking,
                pragmas = CorePragmas.Instanced,
                defines = CoreDefines.ScenePicking,
                includes = SpriteSharedIncludes.DepthOnly,

                // Custom Interpolator Support
                customInterpolators = CoreCustomInterpDescriptors.Common
            };
        }
        #endregion

        #region Includes
        static class SpriteSharedIncludes
        {
            const string kSpriteDepthOnlyPass = "Packages/com.unity.render-pipelines.universal/Editor/2D/ShaderGraph/Includes/SpriteDepthOnlyPass.hlsl";

            public static IncludeCollection DepthOnly = new IncludeCollection
            {
                // Pre-graph
                { CoreIncludes.CorePregraph },
                { CoreIncludes.ShaderGraphPregraph },

                // Post-graph
                { CoreIncludes.CorePostgraph },
                { kSpriteDepthOnlyPass, IncludeLocation.Postgraph },
            };
        }
        #endregion
    }
}
