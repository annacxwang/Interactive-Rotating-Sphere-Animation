/* 
File Name: "vshader53.glsl":
Vertex shader:
  - Per vertex shading for a single point light source;
    distance attenuation is Yet To Be Completed.
  - Entire shading computation is done in the Eye Frame.
*/

#version 150  // YJC: Comment/un-comment this line to resolve compilation errors
              //      due to different settings of the default GLSL version


in vec3 v;
in vec4 vColor;

out vec4 color;
out float y;
uniform mat4 ModelView;
uniform mat4 Projection;

uniform float t;

void main()
{
    vec4 vPosition;
    float a = -0.00000049;
    vPosition = vec4( 0.001 * v.x * t,0.1 +  0.001 * v.y * t + 0.5 * a * t * t,0.001 * v.z * t,1.0);
    y = vPosition.y;
    gl_Position = Projection * ModelView * vPosition;
    color = vColor;
}
