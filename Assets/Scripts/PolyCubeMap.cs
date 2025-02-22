using System.Collections.Generic;
using UnityEngine;

public class PolyCubeMap : MonoBehaviour
{
    public Texture2D PolyCubemapTexture;
    public int SquareletSize = 16;

    public bool ShowParametrization;

    private void Start()
    {
        if (PolyCubemapTexture != null)
        {
            if (ShowParametrization)
            {
                var meshFilter = GetComponent<MeshFilter>();
                if (meshFilter != null)
                {
                    var uvs = new List<Vector3>();
                    meshFilter.mesh.GetUVs(0, uvs);
                    meshFilter.mesh.SetVertices(uvs);
                    meshFilter.mesh.UploadMeshData(false);
                }
            }
            var meshRenderer = GetComponent<MeshRenderer>();
            if (meshRenderer != null)
            {
                var d0x = 1.0f / PolyCubemapTexture.width;
                var d0y = 1.0f / PolyCubemapTexture.height;
                var d1x = 255.0f * SquareletSize / PolyCubemapTexture.width;
                var d1y = 255.0f * SquareletSize / PolyCubemapTexture.height;
                var d2x = 1.0f * SquareletSize / PolyCubemapTexture.width;
                var d2y = 1.0f * SquareletSize / PolyCubemapTexture.height;
                var material = meshRenderer.material;
                material.mainTexture = PolyCubemapTexture;
                material.SetVector("text_coord_normalizer", new Vector4(d0x, d0y, 0f, 0f));
                material.SetVector("text_coord_normalizer_on_TS_times_255", new Vector4(d1x, d1y, 0f, 0f));
                material.SetVector("text_coord_normalizer_on_TS_0", new Vector4(d2x, d2y, 1f, 0f));
                material.SetVector("text_coord_normalizer_on_TS_1", new Vector4(d2x, d2y, 1f, 1f));
                meshRenderer.material = material;
            }
        }
    }
}
