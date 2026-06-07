// PolyCubeMapBakerWindow.cs
// Thin EditorWindow over PolyCubeMapBaker.BakeFromMesh for arbitrary meshes.

using UnityEditor;
using UnityEngine;

public class PolyCubeMapBakerWindow : EditorWindow
{
    private Mesh _mesh;
    private Texture2D _sourceTexture;
    private int _resolution = 8;
    private bool _fillInterior = true;
    private int _squareletSizeIndex = 1; // -> 8
    private int _atlasSizeIndex = 1;     // -> 1024
    private string _outputBaseName = "";

    private static readonly int[] SquareletSizes = { 4, 8, 16 };
    private static readonly string[] SquareletLabels = { "4", "8", "16" };
    private static readonly int[] AtlasSizes = { 512, 1024, 2048 };
    private static readonly string[] AtlasLabels = { "512", "1024", "2048" };

    [MenuItem("PolyCubeMap/Bake/From Mesh…")]
    public static void ShowWindow()
    {
        var w = GetWindow<PolyCubeMapBakerWindow>(true, "PolyCubeMap Baker", true);
        w.minSize = new Vector2(360, 320);
        w.Show();
    }

    private void OnGUI()
    {
        EditorGUILayout.LabelField("Bake a PolyCubeMap from an arbitrary mesh.", EditorStyles.boldLabel);
        EditorGUILayout.Space();

        EditorGUI.BeginChangeCheck();
        var newMesh = (Mesh)EditorGUILayout.ObjectField("Mesh", _mesh, typeof(Mesh), false);
        if (EditorGUI.EndChangeCheck())
        {
            _mesh = newMesh;
            // Auto-populate the output base name from the mesh name.
            if (_mesh != null && string.IsNullOrWhiteSpace(_outputBaseName))
                _outputBaseName = _mesh.name + "_pcm";
        }

        _sourceTexture = (Texture2D)EditorGUILayout.ObjectField(
            new GUIContent("Source Texture", "Optional. If null, the baker uses an HSV checkerboard."),
            _sourceTexture, typeof(Texture2D), false);

        _resolution = EditorGUILayout.IntSlider(
            new GUIContent("Resolution", "Max cubes along the longest axis. Capped at 14 by the shader."),
            _resolution, 2, 14);

        _fillInterior = EditorGUILayout.Toggle(
            new GUIContent("Fill interior", "Flood-fills the inside. Only valid for closed meshes."),
            _fillInterior);

        _squareletSizeIndex = EditorGUILayout.Popup("Squarelet size", _squareletSizeIndex, SquareletLabels);
        _atlasSizeIndex = EditorGUILayout.Popup("Atlas size", _atlasSizeIndex, AtlasLabels);

        if (string.IsNullOrWhiteSpace(_outputBaseName) && _mesh != null)
            _outputBaseName = _mesh.name + "_pcm";
        _outputBaseName = EditorGUILayout.TextField("Output base name", _outputBaseName);

        EditorGUILayout.Space();
        EditorGUILayout.HelpBox(
            "Limits:\n" +
            "  - Resolution is capped at 14 (shader packs LUT as cellX + 16*cellZ).\n" +
            "  - 'Fill interior' assumes a closed/watertight mesh; open meshes will leak.\n" +
            "  - Source mesh must be Read/Write enabled in its importer.\n" +
            "  - Source texture (if provided) must be Read/Write enabled.\n" +
            "  - Cell types 6a/6b and the iterative MIPS optimization are not implemented.",
            MessageType.Info);

        EditorGUILayout.Space();
        using (new EditorGUI.DisabledScope(_mesh == null))
        {
            if (GUILayout.Button("Bake", GUILayout.Height(28)))
            {
                var sz = AtlasSizes[_atlasSizeIndex];
                var sq = SquareletSizes[_squareletSizeIndex];
                var go = PolyCubeMapBaker.BakeFromMesh(_mesh, _sourceTexture, _resolution,
                                                      _fillInterior, sq, sz, sz, _outputBaseName);
                if (go != null) EditorGUIUtility.PingObject(go);
            }
        }
    }
}
