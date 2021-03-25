namespace UnityEditor.Rendering
{
    public static partial class CameraUI
    {
        public partial class Output
        {
            public static void Drawer_Output_AllowDynamicResolution(ISerializedCamera p, Editor owner)
            {
                EditorGUILayout.PropertyField(p.allowDynamicResolution, Styles.allowDynamicResolution);
                p.baseCameraSettings.allowDynamicResolution.boolValue = p.allowDynamicResolution.boolValue;
            }

            public static void Drawer_Output_NormalizedViewPort(ISerializedCamera p, Editor owner)
            {
                EditorGUILayout.PropertyField(p.baseCameraSettings.normalizedViewPortRect, Styles.viewport);
            }
        }
    }
}
