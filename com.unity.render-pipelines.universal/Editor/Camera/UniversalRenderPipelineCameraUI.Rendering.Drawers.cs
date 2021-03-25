using System.Linq;
using NUnit.Framework.Constraints;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.Universal;

namespace UnityEditor.Rendering.Universal
{
    using CED = CoreEditorDrawer<UniversalRenderPipelineSerializedCamera>;

    static partial class UniversalRenderPipelineCameraUI
    {
        public partial class Rendering
        {
            static bool s_PostProcessingWarningShown = false;

            static readonly CED.IDrawer PostProcessingWarningInit = CED.Group(
                (serialized, owner) => s_PostProcessingWarningShown = false
            );

            private static readonly CED.IDrawer PostProcessingWarningDrawer = CED.Conditional(
                (serialized, owner) => IsAnyRendererHasPostProcessingEnabled(serialized, UniversalRenderPipeline.asset) && serialized.renderPostProcessing.boolValue,
                (serialized, owner) =>
                {
                    EditorGUILayout.HelpBox(Styles.disabledPostprocessing, MessageType.Warning);
                    s_PostProcessingWarningShown = true;
                });

            private static readonly CED.IDrawer PostProcessingAAWarningDrawer = CED.Conditional(
                (serialized, owner) => !s_PostProcessingWarningShown && IsAnyRendererHasPostProcessingEnabled(serialized, UniversalRenderPipeline.asset) && (AntialiasingMode)serialized.antialiasing.intValue != AntialiasingMode.None,
                (serialized, owner) =>
                {
                    EditorGUILayout.HelpBox(Styles.disabledPostprocessing, MessageType.Warning);
                    s_PostProcessingWarningShown = true;
                });

            private static readonly CED.IDrawer PostProcessingStopNaNsWarningDrawer = CED.Conditional(
                (serialized, owner) => !s_PostProcessingWarningShown && IsAnyRendererHasPostProcessingEnabled(serialized, UniversalRenderPipeline.asset) && serialized.stopNaNs.boolValue,
                (serialized, owner) =>
                {
                    EditorGUILayout.HelpBox(Styles.disabledPostprocessing, MessageType.Warning);
                    s_PostProcessingWarningShown = true;
                });

            private static readonly CED.IDrawer PostProcessingDitheringWarningDrawer = CED.Conditional(
                (serialized, owner) => !s_PostProcessingWarningShown && IsAnyRendererHasPostProcessingEnabled(serialized, UniversalRenderPipeline.asset) && serialized.dithering.boolValue,
                (serialized, owner) =>
                {
                    EditorGUILayout.HelpBox(Styles.disabledPostprocessing, MessageType.Warning);
                    s_PostProcessingWarningShown = true;
                });

            static readonly CED.IDrawer BaseCameraRenderTypeDrawer = CED.Conditional(
                (serialized, owner) => (CameraRenderType) serialized.cameraType.intValue == CameraRenderType.Base,
                CED.Group(
                    Drawer_Rendering_RenderPostProcessing
                ),
                PostProcessingWarningDrawer,
                CED.Group(
                    Drawer_Rendering_Antialiasing
                ),
                PostProcessingAAWarningDrawer,
                CED.Conditional(
                    (serialized, owner) => !serialized.antialiasing.hasMultipleDifferentValues,
                    CED.Group(
                        GroupOption.Indent,
                        CED.Conditional(
                            (serialized, owner) => (AntialiasingMode) serialized.antialiasing.intValue ==
                                                   AntialiasingMode.SubpixelMorphologicalAntiAliasing,
                            CED.Group(
                                Drawer_Rendering_SMAAQuality
                            )
                        )
                    )
                ),
                CED.Group(
                    CameraUI.Rendering.Drawer_Rendering_StopNaNs
                    ),
                PostProcessingStopNaNsWarningDrawer,
                CED.Conditional(
                    (serialized, owner) => serialized.stopNaNs.boolValue && CoreEditorUtils.buildTargets.Contains(GraphicsDeviceType.OpenGLES2),
                    (serialized, owner) => EditorGUILayout.HelpBox("Stop NaNs has no effect on GLES2 platforms.", MessageType.Warning)
                ),
                CED.Group(
                    CameraUI.Rendering.Drawer_Rendering_Dithering
                ),
                PostProcessingDitheringWarningDrawer,
                CED.Group(
                    Drawer_Rendering_RenderShadows,
                    Drawer_Rendering_Priority,
                    Drawer_Rendering_OpaqueTexture,
                    Drawer_Rendering_DepthTexture
                )
            );

            static readonly CED.IDrawer OverlayCameraRenderTypeDrawer = CED.Conditional(
                (serialized, owner) => (CameraRenderType) serialized.cameraType.intValue == CameraRenderType.Overlay,
                CED.Group(
                    Drawer_Rendering_RenderPostProcessing
                ),
                PostProcessingWarningDrawer,
                CED.Group(
                    Drawer_Rendering_ClearDepth
                )
            );

            public static readonly CED.IDrawer Drawer = CED.FoldoutGroup(
                CameraUI.Rendering.Styles.header,
                CameraUI.Expandable.Rendering,
                k_ExpandedState,
                FoldoutOption.Indent,
                PostProcessingWarningInit,
                CED.Group(
                    Drawer_Rendering_Renderer
                    ),
                BaseCameraRenderTypeDrawer,
                OverlayCameraRenderTypeDrawer,
                CED.Group(
                    Drawer_Rendering_RenderShadows,
                    CameraUI.Rendering.Drawer_Rendering_CullingMask,
                    CameraUI.Rendering.Drawer_Rendering_OcclusionCulling
                )
            );

            static void Drawer_Rendering_Renderer(UniversalRenderPipelineSerializedCamera p, Editor owner)
            {
                var rpAsset = UniversalRenderPipeline.asset;

                int selectedRendererOption = p.renderer.intValue;
                EditorGUI.BeginChangeCheck();

                Rect controlRect = EditorGUILayout.GetControlRect(true);
                EditorGUI.BeginProperty(controlRect, Styles.rendererType, p.renderer);

                EditorGUI.showMixedValue = p.renderer.hasMultipleDifferentValues;
                int selectedRenderer = EditorGUI.IntPopup(controlRect, Styles.rendererType, selectedRendererOption, rpAsset.rendererDisplayList, rpAsset.rendererIndexList);
                EditorGUI.EndProperty();
                if (!rpAsset.ValidateRendererDataList())
                {
                    EditorGUILayout.HelpBox(Styles.noRendererError, MessageType.Error);
                }
                else if (!rpAsset.ValidateRendererData(selectedRendererOption))
                {
                    EditorGUILayout.HelpBox(Styles.missingRendererWarning, MessageType.Warning);
                    var rect = EditorGUI.IndentedRect(EditorGUILayout.GetControlRect());
                    if (GUI.Button(rect, "Select Render Pipeline Asset"))
                    {
                        Selection.activeObject = AssetDatabase.LoadAssetAtPath<UniversalRenderPipelineAsset>(AssetDatabase.GetAssetPath(UniversalRenderPipeline.asset));
                    }
                    GUILayout.Space(5);
                }

                if (EditorGUI.EndChangeCheck())
                    p.renderer.intValue = selectedRenderer;
            }

            static bool IsAnyRendererHasPostProcessingEnabled(UniversalRenderPipelineSerializedCamera p, UniversalRenderPipelineAsset rpAsset)
            {
                int selectedRendererOption = p.renderer.intValue;

                if (selectedRendererOption < -1 || selectedRendererOption > rpAsset.m_RendererDataList.Length || p.renderer.hasMultipleDifferentValues)
                    return false;

                var rendererData = selectedRendererOption == -1 ? rpAsset.m_RendererData : rpAsset.m_RendererDataList[selectedRendererOption];

                var forwardRendererData = rendererData as UniversalRendererData;
                if (forwardRendererData != null && forwardRendererData.postProcessData == null)
                    return true;

                var renderer2DData = rendererData as UnityEngine.Experimental.Rendering.Universal.Renderer2DData;
                return renderer2DData != null && renderer2DData.postProcessData == null;
            }

            static void Drawer_Rendering_Antialiasing(UniversalRenderPipelineSerializedCamera p, Editor owner)
            {
                Rect antiAliasingRect = EditorGUILayout.GetControlRect();
                EditorGUI.BeginProperty(antiAliasingRect, Styles.antialiasing, p.antialiasing);
                {
                    EditorGUI.BeginChangeCheck();
                    int selectedValue = (int)(AntialiasingMode)EditorGUI.EnumPopup(antiAliasingRect, Styles.antialiasing, (AntialiasingMode)p.antialiasing.intValue);
                    if (EditorGUI.EndChangeCheck())
                        p.antialiasing.intValue = selectedValue;
                }
                EditorGUI.EndProperty();
            }

            static void Drawer_Rendering_ClearDepth(UniversalRenderPipelineSerializedCamera p, Editor owner)
            {
                EditorGUILayout.PropertyField(p.clearDepth, Styles.clearDepth);
            }

            static void Drawer_Rendering_RenderShadows(UniversalRenderPipelineSerializedCamera p, Editor owner)
            {
                EditorGUILayout.PropertyField(p.renderShadows, Styles.renderingShadows);
            }

            static void Drawer_Rendering_SMAAQuality(UniversalRenderPipelineSerializedCamera p, Editor owner)
            {
                EditorGUILayout.PropertyField(p.antialiasingQuality, Styles.antialiasingQuality);

                if (CoreEditorUtils.buildTargets.Contains(GraphicsDeviceType.OpenGLES2))
                    EditorGUILayout.HelpBox("Sub-pixel Morphological Anti-Aliasing isn't supported on GLES2 platforms.", MessageType.Warning);
            }

            static void Drawer_Rendering_RenderPostProcessing(UniversalRenderPipelineSerializedCamera p, Editor owner)
            {
                EditorGUILayout.PropertyField(p.renderPostProcessing, Styles.renderPostProcessing);
            }

            static void Drawer_Rendering_Priority(UniversalRenderPipelineSerializedCamera p, Editor owner)
            {
                EditorGUILayout.PropertyField(p.baseCameraSettings.depth, Styles.priority);
            }

            static void Drawer_Rendering_DepthTexture(UniversalRenderPipelineSerializedCamera p, Editor owner)
            {
                EditorGUILayout.PropertyField(p.renderDepth, Styles.requireDepthTexture);
            }

            static void Drawer_Rendering_OpaqueTexture(UniversalRenderPipelineSerializedCamera p, Editor owner)
            {
                EditorGUILayout.PropertyField(p.renderOpaque, Styles.requireOpaqueTexture);
            }
        }
    }
}
