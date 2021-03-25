using System;
using System.Linq;
using System.Runtime.CompilerServices;
using UnityEngine;

namespace UnityEditor.Rendering
{
    public static partial class CameraUI
    {
        public partial class Rendering
        {
            public static void Drawer_Rendering_StopNaNs(ISerializedCamera p, Editor owner)
            {
                EditorGUILayout.PropertyField(p.stopNaNs, Styles.stopNaNs);
            }

            public static void Drawer_Rendering_Dithering(ISerializedCamera p, Editor owner)
            {
                EditorGUILayout.PropertyField(p.dithering, Styles.dithering);
            }

            public static void Drawer_Rendering_CullingMask(ISerializedCamera p, Editor owner)
            {
                EditorGUILayout.PropertyField(p.baseCameraSettings.cullingMask, Styles.cullingMask);
            }

            public static void Drawer_Rendering_OcclusionCulling(ISerializedCamera p, Editor owner)
            {
                EditorGUILayout.PropertyField(p.baseCameraSettings.occlusionCulling, Styles.occlusionCulling);
            }
        }
    }
}
