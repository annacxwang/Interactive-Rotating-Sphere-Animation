/* 
File Name: "fshader53.glsl":
           Fragment Shader
*/

#version 150  // YJC: Comment/un-comment this line to resolve compilation errors
              //      due to different settings of the default GLSL version

in  vec4 color;
in float y;
out vec4 fColor;



void main() 
{ 
    
    if (y<0.1) discard;
    fColor = color;


} 

