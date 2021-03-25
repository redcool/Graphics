using UnityEngine;
using UnityEngine.Rendering.Universal;

namespace UnityEditor.Rendering.Universal
{
    using CED = CoreEditorDrawer<UniversalRenderPipelineSerializedCamera>;

    static partial class UniversalRenderPipelineCameraUI
    {
        static readonly ExpandedState<CameraUI.Expandable, Camera> k_ExpandedState = new ExpandedState<CameraUI.Expandable, Camera>(CameraUI.Expandable.Projection, "URP");

        public static readonly CED.IDrawer SectionProjectionSettings = CED.FoldoutGroup(
            CameraUI.Styles.projectionSettingsHeaderContent,
            CameraUI.Expandable.Projection,
            k_ExpandedState,
            FoldoutOption.Indent,
            CED.Group(
                CameraUI.Drawer_Projection
            ),
            PhysicalCamera.Drawer,
            CED.Group(
                CameraUI.Drawer_FieldClippingPlanes
            )
        );

        public static readonly CED.IDrawer[] Inspector =
        {
            CED.Group(
                Drawer_CameraType
            ),
            SectionProjectionSettings,
            Rendering.Drawer,
            Environment.Drawer,
            Output.Drawer
        };

        static void Drawer_CameraType(UniversalRenderPipelineSerializedCamera p, Editor owner)
        {
            int selectedRenderer = p.renderer.intValue;
            ScriptableRenderer scriptableRenderer = UniversalRenderPipeline.asset.GetRenderer(selectedRenderer);
            bool isDeferred = scriptableRenderer is UniversalRenderer renderer ? renderer.renderingMode == RenderingMode.Deferred : false;

            EditorGUI.BeginChangeCheck();

            CameraRenderType originalCamType = (CameraRenderType)p.cameraType.intValue;
            CameraRenderType camType = (originalCamType != CameraRenderType.Base && isDeferred) ? CameraRenderType.Base : originalCamType;

            camType = (CameraRenderType)EditorGUILayout.EnumPopup(
                Styles.cameraType,
                camType,
                e =>
                {
                    return isDeferred ? (CameraRenderType)e != CameraRenderType.Overlay : true;
                },
                false
            );

            if (EditorGUI.EndChangeCheck() || camType != originalCamType)
            {
                p.cameraType.intValue = (int)camType;
            }
        }
    }
}
