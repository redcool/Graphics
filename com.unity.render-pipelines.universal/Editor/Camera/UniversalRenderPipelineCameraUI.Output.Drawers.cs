using System;
using UnityEngine;
using UnityEngine.Rendering.Universal;

namespace UnityEditor.Rendering.Universal
{
    using CED = CoreEditorDrawer<UniversalRenderPipelineSerializedCamera>;

    static partial class UniversalRenderPipelineCameraUI
    {
        public partial class Output
        {
            public static readonly CED.IDrawer Drawer = CED.Conditional(
                    (serialized, owner) => (CameraRenderType)serialized.cameraType.intValue == CameraRenderType.Base,
                    CED.FoldoutGroup(
                        CameraUI.Output.Styles.header,
                        CameraUI.Expandable.Output,
                        k_ExpandedState,
                        FoldoutOption.Indent,
                        CED.Group(
                            Drawer_Output_TargetTexture
                        ),
                        CED.Conditional(
                            (serialized, owner) => serialized.serializedObject.targetObject is Camera camera && camera.targetTexture == null,
                            CED.Group(
                                Drawer_Output_HDR,
                                Drawer_Output_MSAA,
                                Drawer_Output_AllowDynamicResolution,
                                Drawer_Output_MultiDisplay
                            )
                        ),
#if ENABLE_VR && ENABLE_XR_MODULE
                        CED.Group(Drawer_Output_XRRendering),
#endif
                        CED.Group(
                            Drawer_Output_NormalizedViewPort
                        )
                    )
            );

            static void Drawer_Output_MultiDisplay(UniversalRenderPipelineSerializedCamera p, Editor owner)
            {
                using (var checkScope = new EditorGUI.ChangeCheckScope())
                {
                    p.baseCameraSettings.DrawMultiDisplay();
                    if (checkScope.changed)
                    {
                        UpdateStackCamerasOutput(p, camera =>
                        {
                            bool isChanged = false;
                            // Force same target display
                            int targetDisplay = p.baseCameraSettings.targetDisplay.intValue;
                            if (camera.targetDisplay != targetDisplay)
                            {
                                camera.targetDisplay = targetDisplay;
                                isChanged = true;
                            }

                            // Force same target display
                            StereoTargetEyeMask stereoTargetEye = (StereoTargetEyeMask)p.baseCameraSettings.targetEye.intValue;
                            if (camera.stereoTargetEye != stereoTargetEye)
                            {
                                camera.stereoTargetEye = stereoTargetEye;
                                isChanged = true;
                            }

                            return isChanged;
                        });
                    }
                }
            }

            static void Drawer_Output_AllowDynamicResolution(UniversalRenderPipelineSerializedCamera p, Editor owner)
            {
                using (var checkScope = new EditorGUI.ChangeCheckScope())
                {
                    CameraUI.Output.Drawer_Output_AllowDynamicResolution(p, owner);
                    if (checkScope.changed)
                    {
                        UpdateStackCamerasOutput(p, camera =>
                        {
                            bool allowDynamicResolution = p.allowDynamicResolution.boolValue;

                            Debug.Log($"allowDynamicResolution {allowDynamicResolution} - {camera.name} : {camera.allowDynamicResolution}");

                            if (camera.allowDynamicResolution == p.allowDynamicResolution.boolValue)
                                return false;

                            Debug.Log("Changed");

                            EditorUtility.SetDirty(camera);

                            camera.allowDynamicResolution = allowDynamicResolution;
                            return true;

                            //return false;
                        });
                    }
                }
            }

            static void Drawer_Output_NormalizedViewPort(UniversalRenderPipelineSerializedCamera p, Editor owner)
            {
                using (var checkScope = new EditorGUI.ChangeCheckScope())
                {
                    CameraUI.Output.Drawer_Output_NormalizedViewPort(p, owner);
                    if (checkScope.changed)
                    {
                        UpdateStackCamerasOutput(p, camera =>
                        {
                            Rect rect = p.baseCameraSettings.normalizedViewPortRect.rectValue;
                            if (camera.rect != rect)
                            {
                                camera.rect = p.baseCameraSettings.normalizedViewPortRect.rectValue;
                                return true;
                            }

                            return false;
                        });
                    }
                }
            }

            static void UpdateStackCamerasOutput(UniversalRenderPipelineSerializedCamera p, Func<Camera, bool> updateOutputProperty)
            {
                int cameraCount = p.cameras.arraySize;
                for (int i = 0; i < cameraCount; ++i)
                {
                    SerializedProperty cameraProperty = p.cameras.GetArrayElementAtIndex(i);
                    Camera overlayCamera = cameraProperty.objectReferenceValue as Camera;
                    if (overlayCamera != null)
                    {
                        Undo.RecordObject(overlayCamera, Styles.inspectorOverlayCameraText);
                        if (updateOutputProperty(overlayCamera))
                            EditorUtility.SetDirty(overlayCamera);
                    }
                }
            }

            static void Drawer_Output_TargetTexture(UniversalRenderPipelineSerializedCamera p, Editor owner)
            {
                var rpAsset = UniversalRenderPipeline.asset;
                using (var checkScope = new EditorGUI.ChangeCheckScope())
                {
                    EditorGUILayout.PropertyField(p.baseCameraSettings.targetTexture, Styles.targetTextureLabel);

                    var texture = p.baseCameraSettings.targetTexture.objectReferenceValue as RenderTexture;
                    if (!p.baseCameraSettings.targetTexture.hasMultipleDifferentValues && rpAsset != null)
                    {
                        int pipelineSamplesCount = rpAsset.msaaSampleCount;

                        if (texture && texture.antiAliasing > pipelineSamplesCount)
                        {
                            string pipelineMSAACaps = (pipelineSamplesCount > 1) ? $"is set to support {pipelineSamplesCount}x" : "has MSAA disabled";
                            EditorGUILayout.HelpBox(
                                $"Camera target texture requires {texture.antiAliasing}x MSAA. Universal pipeline {pipelineMSAACaps}.",
                                MessageType.Warning, true);
                        }
                    }

                    if (checkScope.changed)
                    {
                        UpdateStackCamerasOutput(p, camera =>
                        {
                            if (camera.targetTexture == texture)
                                return false;

                            camera.targetTexture = texture;
                            return true;

                        });
                    }
                }
            }

#if ENABLE_VR && ENABLE_XR_MODULE
            static void Drawer_Output_XRRendering(UniversalRenderPipelineSerializedCamera p, Editor owner)
            {
                Rect controlRect = EditorGUILayout.GetControlRect(true);
                EditorGUI.BeginProperty(controlRect, Styles.xrTargetEye, p.allowXRRendering);
                {
                    using (var checkScope = new EditorGUI.ChangeCheckScope())
                    {
                        int selectedValue = !p.allowXRRendering.boolValue ? 0 : 1;
                        bool allowXRRendering = EditorGUI.IntPopup(controlRect, Styles.xrTargetEye, selectedValue, Styles.xrTargetEyeOptions, Styles.xrTargetEyeValues) == 1;
                        if (checkScope.changed)
                            p.allowXRRendering.boolValue = allowXRRendering;
                    }
                }
                EditorGUI.EndProperty();
            }
#endif

            static void Drawer_Output_HDR(UniversalRenderPipelineSerializedCamera p, Editor owner)
            {
                Rect controlRect = EditorGUILayout.GetControlRect(true);
                EditorGUI.BeginProperty(controlRect, Styles.allowHDR, p.baseCameraSettings.HDR);
                {
                    using (var checkScope = new EditorGUI.ChangeCheckScope())
                    {
                        int selectedValue = !p.baseCameraSettings.HDR.boolValue ? 0 : 1;
                        var allowHDR = EditorGUI.IntPopup(controlRect, Styles.allowHDR, selectedValue, Styles.displayedCameraOptions, Styles.cameraOptions) == 1;
                        if (checkScope.changed)
                        {
                            p.baseCameraSettings.HDR.boolValue = allowHDR;
                            UpdateStackCamerasOutput(p, camera =>
                            {
                                if (camera.allowHDR == allowHDR)
                                    return false;

                                camera.allowHDR = allowHDR;
                                return true;
                            });
                        }
                    }
                }
                EditorGUI.EndProperty();
            }

            static void Drawer_Output_MSAA(UniversalRenderPipelineSerializedCamera p, Editor owner)
            {
                Rect controlRect = EditorGUILayout.GetControlRect(true);
                EditorGUI.BeginProperty(controlRect, Styles.allowMSAA, p.baseCameraSettings.allowMSAA);
                {
                    using (var checkScope = new EditorGUI.ChangeCheckScope())
                    {
                        int selectedValue = !p.baseCameraSettings.allowMSAA.boolValue ? 0 : 1;
                        var allowMSAA = EditorGUI.IntPopup(controlRect, Styles.allowMSAA,
                            selectedValue, Styles.displayedCameraOptions, Styles.cameraOptions) == 1;
                        if (checkScope.changed)
                        {
                            p.baseCameraSettings.allowMSAA.boolValue = allowMSAA;
                            UpdateStackCamerasOutput(p, camera =>
                            {
                                if (camera.allowMSAA == allowMSAA)
                                    return false;

                                camera.allowMSAA = allowMSAA;
                                return true;
                            });
                        }
                    }
                }
                EditorGUI.EndProperty();
            }
        }
    }
}
