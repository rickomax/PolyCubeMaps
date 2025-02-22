using System;
using UnityEngine;
using System.IO;
using System.Collections.Generic;
using System.Globalization;
using UnityEditor;

public class OBJParser : MonoBehaviour
{
    [MenuItem("OBJ/Import")]
    private static void LoadOBJFile()
    {
        var objFilePath = EditorUtility.OpenFilePanel("Select OBJ File", "", "obj");

        if (string.IsNullOrEmpty(objFilePath))
        {
            Debug.LogWarning("No file selected.");
            return;
        }

        if (!File.Exists(objFilePath))
        {
            Debug.LogError("OBJ file does not exist at the specified path.");
            return;
        }

        var mesh = ParseOBJFile(objFilePath);

        if (mesh != null)
        {
            SaveMeshAsAsset(mesh, Path.GetFileNameWithoutExtension(objFilePath));
        }
    }

    private static Mesh ParseOBJFile(string filePath)
    {
        var vertices = new List<Vector3>();
        var uv = new List<Vector3>();
        var triangles = new List<int>();

        var lines = File.ReadAllLines(filePath);

        foreach (var line in lines)
        {
            if (line.StartsWith("v ")) // Vertex position
            {
                var parts = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                var x = float.Parse(parts[1], CultureInfo.InvariantCulture);
                var y = float.Parse(parts[2], CultureInfo.InvariantCulture);
                var z = float.Parse(parts[3], CultureInfo.InvariantCulture);
                vertices.Add(new Vector3(x, y, z));
            }
            else if (line.StartsWith("vt ")) // Texture coordinate
            {
                var parts = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                var u = float.Parse(parts[1], CultureInfo.InvariantCulture);
                var v = float.Parse(parts[2], CultureInfo.InvariantCulture);
                var t = float.Parse(parts[3], CultureInfo.InvariantCulture);
                uv.Add(new Vector3(u, v, t));
            }
            else if (line.StartsWith("f ")) // Face
            {
                var parts = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                for (var i = 1; i <= 3; i++)
                {
                    var vertexData = parts[i].Split('/');
                    var vertexIndex = int.Parse(vertexData[0], CultureInfo.InvariantCulture) - 1;
                    triangles.Add(vertexIndex);
                }
            }
        }

        if (vertices.Count == 0)
        {
            Debug.LogError("No vertices found in the OBJ file.");
            return null;
        }

        var mesh = new Mesh();
        mesh.vertices = vertices.ToArray();
        mesh.SetUVs(0, uv);
        mesh.triangles = triangles.ToArray();
        mesh.RecalculateNormals();

        return mesh;
    }

    private static void SaveMeshAsAsset(Mesh mesh, string assetName)
    {
        var path = "Assets/" + assetName + ".asset";
        AssetDatabase.CreateAsset(mesh, path);
        AssetDatabase.SaveAssets();
        Debug.Log("Mesh saved as: " + path);
    }
}