// Original shader:
// https://vcg.isti.cnr.it/polycubemaps/resources/polycubemap.fp

Shader "Custom/PolyCubeMapUnlit"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
        text_coord_normalizer ("text_coord_normalizer Size)", Vector) = (0, 0, 0, 0)
        text_coord_normalizer_on_TS_times_255 ("text_coord_normalizer_on_TS_times_255 Size)", Vector) = (0, 0, 0, 0)
        text_coord_normalizer_on_TS_0 ("text_coord_normalizer_on_TS_0 Size", Vector) = (0, 0, 0, 0)
        text_coord_normalizer_on_TS_1 ("text_coord_normalizer_on_TS_1 Size", Vector) = (0, 0, 0, 0)
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #include "UnityCG.cginc"

            struct appdata
            {
                float4 vertex : POSITION;
                float3 texcoord : TEXCOORD0;
            };

            struct v2f
            {
                float4 vertex : SV_POSITION;
                float3 texcoord : TEXCOORD0;
            };

            sampler2D _MainTex;
            float4 _MainTex_ST;

            // Environment parameters (equivalent to program.env[])
            float4 text_coord_normalizer;                // {1/TEXT_SIZE_X, 1/TEXT_SIZE_Y, 0, 0}
            float4 text_coord_normalizer_on_TS_times_255;// {255*HX, 255*HY, 0, 0}
            float4 text_coord_normalizer_on_TS_0;        // {HX, HY, 1, 0}
            float4 text_coord_normalizer_on_TS_1;        // {HX, HY, 1, 1}

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.texcoord = v.texcoord;
                return o;
            }

            float4 tex2DAdjust(sampler2D tex, float4 uv) {
                return tex2Dbias(tex, uv);
            }

            fixed4 frag (v2f i) : SV_Target
            {
                float4 tc = float4(i.texcoord, 0.0);

                // PARAM consts = {1.99, 0.995, 1, 1};
                float4 consts = {1.99, 0.995, 1, 1};

                // FLR map, tc;
                float4 map = floor(tc);

                // DP3 map.x, map, {1,0,16,0};
                map.x = dot(map, float4(1, 0, 16, 0));

                // MAD map, map, text_coord_normalizer, {0,0,0,-1000};
                map = map * text_coord_normalizer + float4(0, 0, 0, -1000);
                // TXB map, map, texture, 2D ;
                map = tex2DAdjust(_MainTex, map);

                // FRC sub, tc;
                float4 sub = frac(tc);

                // MAD sub, sub, consts.x, -consts.y; 
                sub = sub * consts.x  + (-consts.y); 

                // MUL decoded, map.z, {128.0,64.0,32.0, 8.0}; # /2 /4 /8 /32
                float4 decoded = map.z * float4(128.0,64.0,32.0, 8.0);
                // FRC decoded, decoded;
                decoded = frac(decoded);
                // SUB decoded, decoded, {0.49,0.49,0.49,0.234375};  
                decoded = decoded -  float4(0.49,0.49,0.49,0.234375);

                // SWZ tmp, sub, -z,-y,-x,0;
                float4 tmp = float4(-sub.z, -sub.y, -sub.x, 0.0);
                // CMP sub, decoded.z, sub, tmp;
                sub = (decoded.z < 0) ? sub : tmp;
                // SWZ tmp, sub, y,z,x,0;
                tmp = float4(sub.y, sub.z, sub.x,0);
                // CMP sub, decoded.w, sub, tmp;
                sub = (decoded.w < 0) ? sub : tmp;
                // SWZ tmp, sub, y,z,x,0;
                tmp = float4(sub.y, sub.z, sub.x, 0.0);
                // SUB decoded.w, decoded.w, {0,0,0,0.375};
                decoded.w = decoded.w  - 0.375; // Quickly compute 5th bit
                // CMP sub, decoded.w, sub, tmp;
                sub = (decoded.w < 0) ? sub : tmp;
                // MUL tmp, sub, {-1,-1,1,0};
                tmp = sub * float4(-1, -1, 1,0);
                // CMP sub, decoded.x, sub, tmp;
                sub = (decoded.x < 0) ? sub : tmp;
                // MUL tmp, sub, {1,-1,-1,0};
                tmp = sub * float4(1, -1, -1,0);
                // CMP sub, decoded.y, sub, tmp;
                sub = (decoded.y < 0) ? sub : tmp;

                //SUB decoded, {0.12375,0.24875,0.37375,0.49875}, map.z ; 
                decoded = float4(0.12375, 0.24875, 0.37375, 0.49875) - map.z;

                // MOV res, sub;
                float4 res = sub;

                float4 expa;
                // ADD expa.z, sub.x, sub.z;
                expa.z = sub.x + sub.z;
                // CMP expa.w, expa.z, -sub.z, sub.x;
                expa.w = (expa.z < 0) ? -sub.z : sub.x;
                // ADD expa.y, expa.w, {1,1,1,1};
                expa.y = expa.w + 1;
                // RCP expa.y, expa.y;
                expa.y = 1.0 / expa.y;
                // MUL expa.x, expa.z, expa.y;
                expa.x = expa.z * expa.y;
                // CMP res.x, decoded.x, expa.x, res.x;
                res.x = (decoded.x < 0) ? expa.x : res.x;

                float4 try3;
                // ADD try3.y, sub.y, -expa.w; 
                try3.y = sub.y  + (- expa.w);
                // MUL try3.y, try3.y, expa.y;
                try3.y = try3.y * expa.y;
                // CMP res.y, decoded.y, try3.y, res.y;
                res.y = (decoded.y < 0) ? try3.y : res.y;

                float4 try5;
                // ADD try5.y, sub.y, expa.w;
                try5.y = sub.y + expa.w;
                // SUB try5.z, {1,1,1,1}, expa.w;
                try5.z = 1.0 - expa.w;
                // RCP try5.z, try5.z;
                try5.z = 1.0 / try5.z;
                // MUL try5.y, try5.y, try5.z;
                try5.y = try5.y*  try5.z;
                //CMP res.y, decoded.z, try5.y, res.y;
                res.y = (decoded.z < 0) ? try5.y : res.y;


                float4 expb;
                float4 tma;
                float4 tryA;
                float4 try1;
                //SWZ expb, sub, z,-x,y,1;
                expb = float4(sub.z, -sub.x, sub.y, 1.0);
                // ADD expb, expb, sub.y; 
                expb = expb + float4(sub.y, sub.y, sub.y, sub.y);
                // RCP expb.w, expb.w;
                expb.w = 1.0 / expb.w;
                // CMP tma, expa.z, {0,-1,0,0}, {1, 0,0,0};
                tma = (expa.z < 0) ? float4(0, -1, 0, 0) : float4(1, 0, 0, 0);
                // MUL tma, expa.x, tma;
                tma = expa.x * tma;
                // MAD tma, tma, expb.z, -expb.z;
                tma = tma * expb.z - expb.z;
                //ADD tryA, expb, tma;   
                tryA = expb + tma;
                // MUL tryA, tryA, {-1,1,1,1}; 
                tryA = tryA * float4(-1, 1, 1, 1);
                // CMP try, decoded.z, tryA, expb;
                try1 = (decoded.z < 0) ? tryA : expb;
                // CMP tmp, try.y, {3,0,0,0}, {1,-1,0,0};
                tmp = (try1.y < 0) ?  float4(3, 0, 0, 0) : float4(1, -1, 0, 0);

                //MAD try, try, expb.w, tmp;
                try1 = try1 * expb.w + tmp;
                // CMP try, decoded.y, try, res;
                try1 = (decoded.y < 0) ?  try1 : res;
                // CMP res, res.y, res, try;
                res = (res.y < 0) ? res : try1;

                // MAD res, res, text_coord_normalizer_on_TS_0, text_coord_normalizer_on_TS_1;
                res = res * text_coord_normalizer_on_TS_0 + text_coord_normalizer_on_TS_1;
                // MAD res, map, text_coord_normalizer_on_TS_times_255, res;
                res = map * text_coord_normalizer_on_TS_times_255 + res;

                // TXB result.color, res, texture, 2D ;
                float4 color = tex2DAdjust(_MainTex, res);
                return color;
            }
            ENDCG
        }
    }
}
