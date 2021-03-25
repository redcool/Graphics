using System;
using System.Linq;
using System.Runtime.CompilerServices;
using UnityEngine;

namespace UnityEditor.Rendering
{
    public static partial class CameraUI
    {
        public partial class Environment
        {
            public static void Drawer_Environment_VolumeLayerMask(ISerializedCamera p, Editor owner)
            {
                EditorGUILayout.PropertyField(p.volumeLayerMask, Styles.volumeLayerMask);
            }
        }
    }
}
