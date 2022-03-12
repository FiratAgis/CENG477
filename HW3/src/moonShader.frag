#version 430

in vec3 Position;
in vec3 Normal;
in vec2 TexCoord;
in vec3 LightVector;
in vec3 CameraVector;

uniform vec3 lightPosition;
uniform sampler2D TexColor;
uniform sampler2D MoonTexColor;
uniform sampler2D TexGrey;
uniform float textureOffset;

out vec4 color;

vec3 ambientReflectenceCoefficient = vec3(0.5f, 0.5f, 0.5f);
vec3 ambientLightColor = vec3(0.6f, 0.6f, 0.6f);
vec3 specularReflectenceCoefficient= vec3(1.0f, 1.0f, 1.0f);
vec3 specularLightColor = vec3(1.0f, 1.0f, 1.0f);
float SpecularExponent = 10;
vec3 diffuseReflectenceCoefficient;
vec3 diffuseLightColor = vec3(1.0f);


void main()
{
    vec2 textureCoordinate = TexCoord;
    vec4 texColor = texture(MoonTexColor, textureCoordinate);

    diffuseReflectenceCoefficient = texColor.xyz;
    vec3 ambient = (ambientReflectenceCoefficient * ambientLightColor).xyz;    
    vec3 diffuse = max(dot(Normal, LightVector), 0.0f) * (diffuseReflectenceCoefficient * diffuseLightColor);
    vec3 spec = pow(max(dot(reflect(-LightVector, Normal), CameraVector), 0.0f), SpecularExponent) * (specularReflectenceCoefficient * specularLightColor).xyz;

    color = vec4(clamp((ambient + diffuse + spec) * texColor.xyz, 0.0f, 1.0f), 1.0);
}
