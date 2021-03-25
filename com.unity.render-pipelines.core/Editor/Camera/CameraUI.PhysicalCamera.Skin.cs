using System.Linq;
using UnityEngine;

namespace UnityEditor.Rendering
{
    public static partial class CameraUI
    {
        public partial class PhysicalCamera
        {
            public class Styles
            {
                // Camera Body
                public static readonly GUIContent cameraBody = EditorGUIUtility.TrTextContent("Camera Body");
                public static readonly GUIContent sensorType = EditorGUIUtility.TrTextContent("Sensor Type", "Common sensor sizes. Choose an item to set Sensor Size, or edit Sensor Size for your custom settings.");
                public static readonly string[] apertureFormatNames = CameraEditor.Settings.ApertureFormatNames.ToArray();
                public static readonly Vector2[] apertureFormatValues = CameraEditor.Settings.ApertureFormatValues.ToArray();
                public static readonly int customPresetIndex = apertureFormatNames.Length - 1;

                public static readonly GUIContent sensorSize = EditorGUIUtility.TrTextContent("Sensor Size", "The size of the camera sensor in millimeters.");
                public static readonly GUIContent gateFit = EditorGUIUtility.TrTextContent("Gate Fit", "Determines how the rendered area (resolution gate) fits into the sensor area (film gate).");

                // Lens
                public static readonly GUIContent lens = EditorGUIUtility.TrTextContent("Lens");
                public static readonly GUIContent focalLength = EditorGUIUtility.TrTextContent("Focal Length", "The simulated distance between the lens and the sensor of the physical camera. Larger values give a narrower field of view.");
                public static readonly GUIContent shift = EditorGUIUtility.TrTextContent("Shift", "Offset from the camera sensor. Use these properties to simulate a shift lens. Measured as a multiple of the sensor size.");
            }
        }
    }
}
