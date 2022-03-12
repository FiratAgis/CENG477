#version 430

layout (location = 0) in vec3 VertexPosition;
layout (location = 1) in vec3 VertexNormal;
layout (location = 2) in vec2 VertexTex;

uniform vec3 lightPosition; //parsed
uniform vec3 cameraPosition; //parsed

uniform mat4 ProjectionMatrix; //parsed
uniform mat4 ViewMatrix; //parsed
uniform mat4 NormalMatrix; //parsed
uniform mat4 MVP; //parsed

uniform sampler2D TexColor;
uniform sampler2D TexGrey;
uniform float textureOffset; //parsed
uniform float orbitDegree; //parsed

uniform float heightFactor; //parsed
uniform float imageWidth; //parsed
uniform float imageHeight; //parsed


out vec3 Position;
out vec3 Normal;
out vec2 TexCoord;



out vec3 LightVector;
out vec3 CameraVector;



void main()
{
    vec3 n  = normalize((MVP * vec4(VertexNormal, 1.0f)).xyz - (MVP * vec4(VertexPosition, 1.0f)).xyz);

    vec4 texGrey = texture(TexGrey, VertexTex);
    float height = texGrey.r * (heightFactor);
   
    vec3 dummyVP = VertexPosition + height * VertexNormal ;

    LightVector = normalize(lightPosition - (MVP * vec4(dummyVP, 1.0f)).xyz);
    CameraVector = normalize(cameraPosition - (MVP * vec4(dummyVP, 1.0f)).xyz);

    Normal = n;
    TexCoord = VertexTex;
    Position = (MVP * vec4(dummyVP, 1.0f)).xyz ;
    gl_Position = MVP * vec4(dummyVP, 1.0f); // this is a placeholder. It does not correctly set the position

}