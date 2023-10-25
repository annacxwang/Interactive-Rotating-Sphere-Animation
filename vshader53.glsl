/* 
File Name: "vshader53.glsl":
Vertex shader:
  - Per vertex shading for a single point light source;
    distance attenuation is Yet To Be Completed.
  - Entire shading computation is done in the Eye Frame.
*/

#version 150  // YJC: Comment/un-comment this line to resolve compilation errors
              //      due to different settings of the default GLSL version

in  vec4 vPosition;
in  vec3 vNormal;
in vec4 vColor;

in vec2 vTexCoord;

out vec4 color;
out float z;
out vec2 texCoord;
out float texFloat;
out vec2 latticeCoord;

//vec4 AmbientProduct, DiffuseProduct, SpecularProduct;
uniform vec4 MaterialAmbient, MaterialDiffuse, MaterialSpecular;
uniform vec4 DirectionalLightAmbient, DirectionalLightDiffuse, DirectionalLightSpecular;
uniform vec4 PositionalLightAmbient, PositionalLightDiffuse, PositionalLightSpecular;
uniform mat4 ModelView;
uniform mat4 Projection;
uniform mat3 Normal_Matrix;
uniform vec3 LightDirection; 
uniform vec3 Towards;
uniform vec4 LightPosition;   // Must be in Eye Frame
uniform float Shininess;

uniform int lighting; //0: no lighting; 1: with lighting
uniform int lightType; //0: directional light; 1: point source; 2: spot light

uniform float ConstAtt;  // Constant Attenuation
uniform float LinearAtt; // Linear Attenuation
uniform float QuadAtt;   // Quadratic Attenuation
uniform float Exponent;
uniform float Cutoff;

uniform int texType; // 0: no tex; 1: 1d texture; 2 : 2d texture

uniform int if_object; // 1: object frame, 0: eye frame
uniform int object;  //0:axes, 1: floor, 2:sphere,3:shade
uniform int vert;
uniform int slanted;

uniform int lattice;
uniform int upright;
uniform int tilted; 

void main()
{
    // Transform vertex  position into eye coordinates
    vec3 pos = (ModelView * vPosition).xyz;
	vec3 eye = vec3(0.0,0.0,0.0); 
    z = length(eye - pos);
    if(lighting == 0){
        color = vColor;
    }
    else{
    vec4 GlobalAmbient= vec4(1.0,1.0,1.0,1.0);
    vec4 DirectionalAmbientProduct, DirectionalDiffuseProduct, DirectionalSpecularProduct, GlobalAmbientProduct;
    vec4 PositionalAmbientProduct,PositionalDiffuseProduct,PositionalSpecularProduct;
    // Compute products
    DirectionalAmbientProduct = DirectionalLightAmbient * MaterialAmbient;
    DirectionalDiffuseProduct = DirectionalLightDiffuse * MaterialDiffuse;
    DirectionalSpecularProduct = DirectionalLightSpecular * MaterialSpecular;
    GlobalAmbientProduct = GlobalAmbient *MaterialAmbient;

    PositionalAmbientProduct = PositionalLightAmbient * MaterialAmbient;
    PositionalDiffuseProduct = PositionalLightDiffuse * MaterialDiffuse;
    PositionalSpecularProduct = PositionalLightSpecular * MaterialSpecular;

    vec3 L_Directional = normalize(-LightDirection);
    vec3 L_Positional = normalize( LightPosition.xyz - pos );

    vec3 E = normalize( -pos );
    vec3 H_Directional = normalize( L_Directional + E );
    vec3 H_Positional = normalize( L_Positional + E );

    // Transform vertex normal into eye coordinates
      // vec3 N = normalize( ModelView*vec4(vNormal, 0.0) ).xyz;
    vec3 N = normalize(Normal_Matrix * vNormal);

// YJC Note: N must use the one pointing *toward* the viewer
   if ( dot(N, E) < 0 ) {N = -N;}
   
     /* Compute directional illumination equation terms */
    vec4 ambient = DirectionalAmbientProduct;

    float d = max( dot(L_Directional, N), 0.0 );
    vec4  diffuse = DirectionalDiffuseProduct *d;

    float s = pow( max(dot(N, H_Directional), 0.0), Shininess );
    vec4  specular =  DirectionalSpecularProduct *s;
    
    if( dot(L_Directional, N) < 0.0 ) {
	    specular = vec4(0.0, 0.0, 0.0, 1.0); //perpendicular
    } 


    float attenuation;

    attenuation = 1.0;
    color = GlobalAmbientProduct + attenuation * (ambient + diffuse + specular);

    if(lightType !=0) { /*positional light*/
         ambient = PositionalAmbientProduct;
         d = max( dot(L_Positional, N), 0.0 );
         diffuse = d * PositionalDiffuseProduct;
         s = pow( max(dot(N, H_Directional), 0.0), Shininess );
         specular = s * PositionalSpecularProduct;
        if( dot(L_Positional, N) < 0.0 ) {
	        specular = vec4(0.0, 0.0, 0.0, 1.0); //perpendicular
        } 
        float dist = length(LightPosition.xyz-pos);
        attenuation=1/(ConstAtt + LinearAtt*dist + QuadAtt*dist*dist);
        if (lightType ==2){ /* spot light*/
            vec3 L_f = normalize(Towards.xyz - LightPosition.xyz);
            float spot = dot(L_f,-L_Positional);
            if (spot >= cos(Cutoff)) attenuation = attenuation*pow(spot,Exponent);
            else attenuation = 0;
        }
        color = color +  attenuation * (ambient + diffuse + specular);
    }
}
    gl_Position = Projection * ModelView * vPosition;
    if (object == 1) {
        texCoord = vTexCoord;}
    else{
    /* Computes texture coordinates */
        vec3 pos_in_frame;
        if(if_object == 1){
            pos_in_frame = vPosition.xyz;
        }
        else{
            pos_in_frame = pos;
        }
        if(texType == 1){
            if (vert == 1){
                texFloat = pos_in_frame[0] * 2.5; // s = 2.5x
            }
            else if (slanted ==1){
                texFloat = (pos_in_frame[0] + pos_in_frame[1] + pos_in_frame[2]) * 1.5; //s = 1.5(x + y + z)
            }
        }
        else if (texType == 2){
            if (vert == 1){
                texCoord[0] = (pos_in_frame[0] + 1)  * 0.75; // s = 0.75(x + 1) 
                texCoord[1] = (pos_in_frame[1] + 1)  * 0.75; // t = 0.75(y + 1) 
            }
            else if (slanted ==1){
                texCoord[0] = (pos_in_frame[0] + pos_in_frame[1] + pos_in_frame[2]) * 0.45; //s = 0.45(x + y + z) 
                texCoord[1] = (pos_in_frame[0] - pos_in_frame[1] + pos_in_frame[2]) * 0.45; //t = 0.45(x - y + z) 
            }
            
        }
        /* Computes lattice coordinates*/
        if(lattice == 1){
            pos_in_frame = vPosition.xyz; //using object space coordinates
            if(upright ==1){
                latticeCoord[0] = (pos_in_frame[0] + 1)  * 0.5; latticeCoord[1] = (pos_in_frame[1] + 1)  * 0.5;
            //s = 0.5(x + 1) and t = 0.5(y + 1)
            }
            else if (tilted ==1){
                latticeCoord[0] = (pos_in_frame[0] + pos_in_frame[1] + pos_in_frame[2]) * 0.3;
                latticeCoord[1] = (pos_in_frame[0] - pos_in_frame[1] + pos_in_frame[2]) * 0.3;
            // s = 0.3(x + y + z) and t = 0.3(x − y + z),
            }
        }
    }

}
