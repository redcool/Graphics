using System.Linq;
using UnityEngine;

namespace UnityEditor.Rendering
{
    public static partial class CameraUI
    {
        public partial class Output
        {
            public class Styles
            {
                public static readonly GUIContent header = EditorGUIUtility.TrTextContent("Output", "These settings control how the camera output is formatted.");

#if ENABLE_MULTIPLE_DISPLAYS
                public static readonly GUIContent targetDisplay = EditorGUIUtility.TrTextContent("Target Display");
#endif

                public static readonly GUIContent viewport = EditorGUIUtility.TrTextContent("Viewport Rect", "Four values that indicate where on the screen HDRP draws this Camera view. Measured in Viewport Coordinates (values in the range of [0, 1]).");
                public static readonly GUIContent allowDynamicResolution = EditorGUIUtility.TrTextContent("Allow Dynamic Resolution", "Whether to support dynamic resolution.");
            }
        }
    }
}
