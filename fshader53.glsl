/* 
File Name: "fshader53.glsl":
           Fragment Shader
*/

#version 150  // YJC: Comment/un-comment this line to resolve compilation errors
              //      due to different settings of the default GLSL version

in  vec4 color;
in float z;
in vec2 texCoord;
in float texFloat;
in vec2 latticeCoord;
out vec4 fColor;

uniform int fogType;
uniform sampler1D texture_1D;
uniform sampler2D texture_2D;
uniform int texType; // 0: no tex; 1: 1d texture; 2 : 2d texture
uniform int object;
uniform int lattice;

uniform int if_wireframe;


void main() 
{ 
    fColor = color;

    float tex_f;

    if (if_wireframe==0 || object ==1){ //no texture for wire-frame sphere
        if(texType == 1){
            fColor = fColor * texture( texture_1D, texFloat );  
        }
        else if (texType ==2){
            vec4 tex_color = texture( texture_2D, texCoord );
            if (object ==2 && tex_color[0] == 0)  tex_color = vec4(0.9, 0.1, 0.1, 1.0);
            fColor = fColor * tex_color;  

        }
    }

    if(lattice ==1 && object > 1 ){ //lattice for shadow and sphere
        float s = latticeCoord[0];
        float t = latticeCoord[1];
        if ((fract (4*s) < 0.35) && (fract (4*t) < 0.35)) discard;
    }

    float f;
    float f_start = 0.0;
    float f_end = 18.0;
    float f_d = 0.09;
    vec4 fog_color = vec4(0.7,0.7,0.7,0.5);
    float temp_a = fColor.a;
    if(fogType ==1){
        f = (f_end - z)/(f_end - f_start);
    }
    else if (fogType ==2){
        f = 1/exp(f_d*z);
    }
    else if (fogType ==3){
        f = 1/exp(f_d*z*f_d*z);
    }
    if (fogType != 0){
        clamp(f,0.0,1.0);
        fColor = mix(fog_color,fColor,f);
        fColor.a = temp_a; /* preserves the alpha value before and after fog effect */
    }


} 

