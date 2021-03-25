using System.Linq;
using UnityEngine;

namespace UnityEditor.Rendering
{
    public static partial class CameraUI
    {
        public partial class Environment
        {
            public class Styles
            {
                public static readonly GUIContent header = EditorGUIUtility.TrTextContent("Environment", "These settings control what the camera background looks like.");

                public static readonly GUIContent volumeLayerMask = EditorGUIUtility.TrTextContent("Volume Mask", "This camera will only be affected by volumes in the selected scene-layers.");
            }
        }
    }
}
