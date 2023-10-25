
#include "Angel-yjc.h"
#include <string>
#include <fstream>

#define PI 3.14159265358979323846

typedef Angel::vec4  color4;
typedef Angel::vec3  point3;

GLuint Angel::InitShader(const char* vShaderFile, const char* fShaderFile);

GLuint program;       /* shader program object id */
GLuint program2;
GLuint sphere_buffer;   /* vertex buffer object id for sphere */
GLuint floor_buffer;  /* vertex buffer object id for floor */
GLuint axes_buffer;  /* vertex buffer object id for axes */
GLuint shadow_buffer;  /* vertex buffer object id for axes */

GLuint sphere_flat_buffer;
GLuint sphere_smooth_buffer;
GLuint floor_shading_buffer;

GLuint floor_tex_buffer;
GLuint particle_buffer;
// Projection transformation parameters
GLfloat  fovy = 45.0;  // Field-of-view in Y direction angle (in degrees)
GLfloat  aspect;       // Viewport aspect ratio
GLfloat  zNear = 0.5, zFar = 13.0;

GLfloat angle = 0.0; // rotation angle in degrees
GLfloat radius = 1.0;
GLfloat theta = 0.5f; //rotation speed in degrees
vec4 init_eye(7.0, 3.0, -10.0, 1.0); // initial viewer position
//vec4 init_eye(3.0, 2.0, 0.0, 1.0); // initial viewer position
vec4 eye = init_eye;               // current viewer position

int animationFlag = 0; // 1: animation; 0: non-animation. 
int beginFlag = 0; 
int sphereFlag = 1; //0:wire-frame sphere; 1:solid sphere
int shadowFlag = 1;
int lightFlag = 1;
int shadingFlag = 1; //0:smooth shading; 1:flat shading
int lightType = 1; //0: directional light; 1: point source; 2: spot light

int fogFlag = 0; //0: no fog ; 1: linear ; 2: exponential ;3 : exponential square
int blendFlag = 1; // 0:no blending; 1: blending shadow
int floorTexFlag = 2;// 0:no texture; 2: 2d texture
int sphereTexFlag = 2;  // 0:no texture; 1: stripe texture; 2: 2d texture

int objFlag = 1; // 1: object space ; 0: eye space
int verticalFlag = 1;
int slantedFlag = 0;

int latticeFlag = 0;
int uprightFlag = 1;
int tiltedFlag = 0;

int firework = 1;
float time_elapsed =0.0f;
float time_start = 0.0f;
int reset_time = 0;
float T_max = 9999.0f; // to be changed

point3* sphere_points;
color4* sphere_colors;
vec3* sphere_flat_normals;
vec3* sphere_smooth_normals;

color4* sphere_shadow_colors;
int sphere_NumVertices;

const int floor_NumVertices = 6; //(1 face)*(2 triangles/face)*(3 vertices/triangle)
point3 floor_points[floor_NumVertices]; // positions for all vertices
color4 floor_colors[floor_NumVertices]; // colors for all vertices
vec3 floor_normals[floor_NumVertices];

const int NumParticles = 300; //NumParticles indicates # of particles
vec3 particle_velocity[NumParticles];
color4 particle_color[NumParticles];

const int axes_NumVertices = 6; //(2 vertices)*(3 axes)
point3 axes_points[axes_NumVertices]; // positions for all vertices
color4 axes_colors[axes_NumVertices]; // colors for all vertices

point3 turningPts[3] = {point3(-4.0,1.0,4.0),point3(-1.0,1.0,-4.0),point3(3.0,1.0,5.0)};
int seg = 0;
point3 origin = turningPts[seg]; 
point3 center = origin;//set starting point A

point3 transVectors[3];
point3 rotationVectors[3];
point3 currRV;

mat4 accumulatedRM = identity();

point3 L(-14.0, 12.0, -3.0);

mat4 shadow_proj(L.y, 0.0f, 0.0f, 0.0f,
    -L.x, 0.0f, -L.z, -1.0f,
    0.0f, 0.0f, L.y, 0.0f,
    0.0f, 0.0f, 0.0f, L.y);
// q = (−f(x, y, z), 0, −g(x, y, z), −h(x, y, z)) s.t. w >0

// Vertices of the floor
point3 floor_vertices[4] = {
    point3(5, 0, 8), //a
    point3(5,  0,  -4), //b
    point3(-5, 0,  -4), // c
    point3(-5, 0,  8) //d
};

point3 axes_vertices[6] = {
    //x-axis
    point3(0, 0.02f, 0),
    point3(100,  0.02f, 0),
    //point3(200,  0.02f, 0),
    //y-axis
    point3(0, 0,  0),
    point3(0, 100,  0),
   // point3(0, 200,  0),
    //z-axis
    point3(0, 0.02f,  0),
    point3(0, 0.02f,  100),
    //point3(0, 0.02f,  200)
};


// RGBA colors
color4 vertex_colors[8] = {
    color4( 0.0, 0.0, 0.0,0.0),  // black
    color4( 1.0, 0.0, 0.0,0.0),  // red
    color4( 1.0, 1.0, 0.0,0.0),  // yellow
    color4( 0.0, 1.0, 0.0,0.0),  // green
    color4( 0.0, 0.0, 1.0,0.0),  // blue
    color4( 1.0, 0.0, 1.0,0.0),  // magenta
    color4( 1.0, 1.0, 1.0,0.0),  // white
    color4( 0.0, 1.0, 1.0,0.0)   // cyan
};

vec2 floor_texCoord[6] = {
  vec2(5.0, 6.0),  // for a
  vec2(5.0, 0.0),  // for b
  vec2(0.0, 0.0),  // for c

  vec2(0.0, 0.0),  // for c
  vec2(0.0, 6.0),  // for d
  vec2(5.0, 6.0),  // for a 
};

void particles() {
    for (int i = 0; i < NumParticles; i++) {
        particle_velocity[i] = vec3(2.0 * ((rand() % 256) / 256.0 - 0.5), 
            1.2 * 2.0 * ((rand() % 256) / 256.0),// [0,2.4]
            2.0 * ((rand() % 256) / 256.0 - 0.5)); //[-1,1]
        particle_color[i] = color4(rand() % 256 / 256.0 , rand() % 256 / 256.0
            , rand() % 256 / 256.0,1.0);

    }
}
//-------------------------------
// 
// generate 2 triangles: 6 vertices and 6 colors

void floor()
{

    floor_colors[0] = vertex_colors[3]; floor_points[0] = floor_vertices[0]; floor_normals[0] = vec3(0.0f,1.0f,0.0f);
    floor_colors[1] = vertex_colors[3]; floor_points[1] = floor_vertices[1]; floor_normals[1] = vec3(0.0f, 1.0f, 0.0f);
    floor_colors[2] = vertex_colors[3]; floor_points[2] = floor_vertices[2]; floor_normals[2] = vec3(0.0f, 1.0f, 0.0f);

    floor_colors[3] = vertex_colors[3]; floor_points[3] = floor_vertices[0]; floor_normals[3] = vec3(0.0f, 1.0f, 0.0f);
    floor_colors[4] = vertex_colors[3]; floor_points[4] = floor_vertices[3]; floor_normals[4] = vec3(0.0f, 1.0f, 0.0f);
    floor_colors[5] = vertex_colors[3]; floor_points[5] = floor_vertices[2]; floor_normals[5] = vec3(0.0f, 1.0f, 0.0f);
}
// generate 3 axes : 6 vertices and 6 colors

void axes() {
    //x axis

    axes_colors[0] = vertex_colors[1]; axes_points[0] = axes_vertices[0];
    axes_colors[1] = vertex_colors[1]; axes_points[1] = axes_vertices[1];
    //axes_colors[2] = vertex_colors[1]; axes_points[2] = axes_vertices[2];

    //y axis
    
    axes_colors[2] = vertex_colors[5]; axes_points[2] = axes_vertices[2];
    axes_colors[3] = vertex_colors[5]; axes_points[3] = axes_vertices[3];
    //axes_colors[5] = vertex_colors[5]; axes_points[5] = axes_vertices[5];
    // z axis
    axes_colors[4] = vertex_colors[4]; axes_points[4] = axes_vertices[4];
    axes_colors[5] = vertex_colors[4]; axes_points[5] = axes_vertices[5];
    //axes_colors[8] = vertex_colors[4]; axes_points[8] = axes_vertices[8];
    
}

//function to compute distance btw 2 points

float distance(point3 p1, point3 p2) {
    float x = p1.x - p2.x;
    float y = p1.y - p2.y;
    float z = p1.z - p2.z;
    return sqrt(x*x+y*y+z*z);
}

// direction vector from p1 to p2
point3 computeTranslationVector(point3 p1, point3 p2) {
    point3 trans;
    trans.x = p2.x - p1.x;
    trans.y = p2.y - p1.y;
    trans.z = p2.z - p1.z;

    float norm = sqrt(trans.x*trans.x + trans.y *trans.y + trans.z *trans.z);
    trans.x /= norm;
    trans.y /= norm;
    trans.z /= norm;

    return trans;
}

//compute cross product of 2 vectors
point3 ComputeXProduct(point3 u, point3 v) {
    point3 n;
    n.x = u.y * v.z - u.z * v.y;
    n.y = u.z * v.x - u.x * v.z;
    n.z = u.x * v.y - u.y * v.x;
    return n;
}

//computes rotaionvector by the cross product of y-axis and translation vectors
void computeRotationVector() {

    for (int i = 0; i < 3; i++) {
        transVectors[i] = computeTranslationVector(turningPts[i],turningPts[(i+1)%3]);
    }
    for (int j = 0; j < 3; j++) {
        rotationVectors[j] = ComputeXProduct(point3(0,1.0,0),transVectors[j]);
        //rotationVectors[j] = cross((point3(0, 1.0, 0), transVectors[j]));
    }
}

void computeFlatNormals(int triangleNo) {
    
    for (int i = 0; i < triangleNo; i++) {
        /*Compute and store flat normals*/
        vec3 u = sphere_points[3*i + 1] - sphere_points[3*i];
        vec3 v = sphere_points[3*i + 2] - sphere_points[3*i];

        vec3 normal = normalize(cross(u, v));

        sphere_flat_normals[3*i] = normal;
        sphere_flat_normals[3*i + 1] = normal;
        sphere_flat_normals[3*i + 2] = normal;


    }

}

//load data from input file and store in matrix
void load() {
    std::cout << "Please enter name of input file(sphere.8/128/256/1024):";

    std::string fname, myline;
    std::cin >> fname;
    std::ifstream myfile;
    myfile.open(fname);

    if (myfile.is_open()) {
        int numTriangles, vertice;
        myfile >> numTriangles; //number of triangles
        sphere_NumVertices = 3 * numTriangles;//number of vertices
        sphere_points = new point3[sphere_NumVertices];
        sphere_colors = new color4[sphere_NumVertices];
        sphere_flat_normals = new vec3[sphere_NumVertices];
        sphere_smooth_normals = new vec3[sphere_NumVertices];
        //sphere_shadow_points = new point3[sphere_NumVertices];
        sphere_shadow_colors = new color4[sphere_NumVertices];
        float curPoint[3];
        int count = 0;
        
        //loop to load point locations
        for (int i = 0; i < numTriangles; i++) {
            myfile >> vertice;
            //loop to load each point
            for (int j = 0; j < vertice; j++) {
                for (int m = 0; m < 3; m++) {
                    myfile >> curPoint[m];
                }
                sphere_points[count] = point3(curPoint[0], curPoint[1], curPoint[2]);

                /*Compute and store smooth normals*/
                sphere_smooth_normals[count] = normalize(sphere_points[count]);

                count++;

            }

        }
       
        // loop to load all vertices colors 
        for (int i = 0; i < count; i++) {
            sphere_colors[i] = color4(1.0f, 0.84f, 0.0f, 0.0f);
            sphere_shadow_colors[i] = color4(0.25f, 0.25f, 0.25f, 0.65f);
        }
         computeFlatNormals(numTriangles);


    }
    else {
        std::cerr << "Couldn't open file " << fname << " please start over\n";
        exit(0);
    }
}

    /*----- Shader Lighting Parameters -----*/
color4 global_ambient_light(1.0, 1.0, 1.0, 1.0);
/*directional parameters*/


/// a directional (distant) light source with black ambient color (0.0, 0.0, 0.0, 1.0), diffuse color
//(0.8, 0.8, 0.8, 1.0), specular color(0.2, 0.2, 0.2, 1.0), and direction(0.1, 0.0, −1.0, 0.0) in the eye
//coordinate system
color4 directional_light_ambient(0.0f, 0.0f, 0.0f, 1.0f);
    color4 directional_light_diffuse(0.8f, 0.8f, 0.8f, 1.0f);
    color4 directional_light_specular(0.2f, 0.2f, 0.2f, 1.0f);
    vec3 light_direction(0.1f, 0.0f, -1.0f); //in eye frame
/*positional parameters*/
//with white diffuse and specular color(1.0, 1.0, 1.0, 1.0), black ambient color(0.0, 0.0, 0.0, 1.0), with position at L = (−14.0, 12.0, −3.0, 1.0)
//(in the usual world coordinate system).

    color4 positional_light_ambient(0.0f, 0.0f, 0.0f, 1.0f);
    color4 positional_light_diffuse(1.0f, 1.0f, 1.0f, 1.0f);
    color4 positional_light_specular(1.0f, 1.0f, 1.0f, 1.0f);
    float const_att = 2.0f;
    float linear_att = 0.01f;
    float quad_att = 0.001f;
    vec4 light_position(-14.0f, 12.0f, -3.0f,1.0f);
    // In World frame.
    // Needs to transform it to Eye Frame
    // before sending it to the shader(s).
/*spot light parameters*/
    // Spot Light, whose direction is from its position toward the point (−6.0, 0.0, −4.5, 1.0) (againin the usual world coordinate system), 
    //with the exponent value 15.0 and the cutoff angle 20.0
    vec4 towards(-6.0, 0.0, -4.5,1.0);
    float exponent = 15.0;
    float cutoff = 20.0f *PI / 180.0f;


   /* a green diffuse color(0.0, 1.0, 0.0, 1.0),
        with ambient color(0.2, 0.2, 0.2, 1.0) and specular color(0.0, 0.0, 0.0, 1.0), and give your sphere a
        golden yellow diffuseand specular color(1.0, 0.84, 0.0, 1.0), with ambient color(0.2, 0.2, 0.2, 1.0)
        and a shininess coefficient of 125.0*/
    color4 floor_ambient(0.2f, 0.2f, 0.2f, 1.0f);
    color4 floor_diffuse(0.0f, 1.0f, 0.0f, 1.0f);
    color4 floor_specular(0.0f, 0.0f, 0.0f, 1.0f);

    color4 sphere_ambient(0.2f, 0.2f, 0.2f, 1.0f);
    color4 sphere_diffuse(1.0f, 0.84f, 0.0f, 1.0f);
    color4 sphere_specular(1.0f, 0.84f, 0.0f, 1.0f);
    float  sphere_shininess = 125.0f;

    void SetUp_Lighting_Uniform_Vars(mat4 mv);

    void SetUp_Material(int mode);

    int Index = 0;


    /* global definitions for constants and global image arrays */

#define ImageWidth  32
#define ImageHeight 32
    GLubyte Image[ImageHeight][ImageWidth][4];

#define	stripeImageWidth 32
    GLubyte stripeImage[4 * stripeImageWidth];
    static GLuint texCheck;
    static GLuint texStripe;

    /*************************************************************
    void image_set_up(void):
      generate checkerboard and stripe images.

    * Inside init(), call this function and set up texture objects
      for texture mapping.
      (init() is called from main() before calling glutMainLoop().)
    ***************************************************************/
    void image_set_up(void)
    {
        int i, j, c;

        /* --- Generate checkerboard image to the image array ---*/
        for (i = 0; i < ImageHeight; i++)
            for (j = 0; j < ImageWidth; j++)
            {
                c = (((i & 0x8) == 0) ^ ((j & 0x8) == 0));

                if (c == 1) /* white */
                {
                    c = 255;
                    Image[i][j][0] = (GLubyte)c;
                    Image[i][j][1] = (GLubyte)c;
                    Image[i][j][2] = (GLubyte)c;
                }
                else  /* green */
                {
                    Image[i][j][0] = (GLubyte)0;
                    Image[i][j][1] = (GLubyte)150;
                    Image[i][j][2] = (GLubyte)0;
                }

                Image[i][j][3] = (GLubyte)255;
            }

        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

        /*--- Generate 1D stripe image to array stripeImage[] ---*/
        for (j = 0; j < stripeImageWidth; j++) {
            /* When j <= 4, the color is (255, 0, 0),   i.e., red stripe/line.
               When j > 4,  the color is (255, 255, 0), i.e., yellow remaining texture
             */
            stripeImage[4 * j] = (GLubyte)255;
            stripeImage[4 * j + 1] = (GLubyte)((j > 4) ? 255 : 0);
            stripeImage[4 * j + 2] = (GLubyte)0;
            stripeImage[4 * j + 3] = (GLubyte)255;
        }

        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        /*----------- End 1D stripe image ----------------*/

        /*--- texture mapping set-up is to be done in
              init() (set up texture objects),
              display() (activate the texture object to be used, etc.)
              and in shaders.
         ---*/

    } /* end function */


//----------------------------------------------------------------------------
// OpenGL initialization
void init()
{   
    /*set up texture objects for check board*/
    image_set_up();
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    /*--- Create and Initialize a check texture object ---*/
    glGenTextures(1, &texCheck);      // Generate texture obj name(s)

    glActiveTexture(GL_TEXTURE0);  // Set the active texture unit to be 0 
    glBindTexture(GL_TEXTURE_2D, texCheck); // Bind the texture to this texture unit
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, ImageWidth, ImageHeight,
        0, GL_RGBA, GL_UNSIGNED_BYTE, Image);

    /*--- Create and Initialize a stripe texture object ---*/
    glGenTextures(1, &texStripe);      // Generate texture obj name(s)

    glActiveTexture(GL_TEXTURE1);  // Set the active texture unit to be 0 
    glBindTexture(GL_TEXTURE_1D, texStripe); // Bind the texture to this texture unit
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    //glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, stripeImageWidth,
        0, GL_RGBA, GL_UNSIGNED_BYTE, stripeImage);


    glutIdleFunc(NULL);

    particles();
    glGenBuffers(1, &particle_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, particle_buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vec3) * NumParticles + sizeof(color4) * NumParticles ,NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0,sizeof(vec3) * NumParticles, particle_velocity);
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(vec3) * NumParticles, sizeof(color4) * NumParticles,
        particle_color);

    floor();     
 // Create and initialize a vertex buffer object for floor, to be used in display()
    glGenBuffers(1, &floor_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, floor_buffer);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(point3) * floor_NumVertices + sizeof(color4) * floor_NumVertices + sizeof(vec2) * floor_NumVertices,
        NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0,
        sizeof(point3) * floor_NumVertices, floor_points);
    glBufferSubData(GL_ARRAY_BUFFER,
        sizeof(point3) * floor_NumVertices,
        sizeof(color4) * floor_NumVertices,
        floor_colors);
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3) * floor_NumVertices + sizeof(color4) * floor_NumVertices,
        sizeof(vec2) * floor_NumVertices, floor_texCoord);

    // Create and initialize a vertex buffer object for floor w/ shding
    glGenBuffers(1, &floor_shading_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, floor_shading_buffer);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(point3) * floor_NumVertices + sizeof(vec3) * floor_NumVertices + sizeof(vec2) * floor_NumVertices,
        NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0,
        sizeof(point3) * floor_NumVertices, floor_points);
    glBufferSubData(GL_ARRAY_BUFFER,
        sizeof(point3) * floor_NumVertices,
        sizeof(vec3) * floor_NumVertices,
        floor_normals);
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(point3) * floor_NumVertices + sizeof(vec3) * floor_NumVertices,
        sizeof(vec2) * floor_NumVertices, floor_texCoord);

    axes();
    // Create and initialize a vertex buffer object for axes

    glGenBuffers(1, &axes_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, axes_buffer);

    glBufferData(GL_ARRAY_BUFFER,
        sizeof(point3) * axes_NumVertices + sizeof(color4) * axes_NumVertices,
        NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0,
        sizeof(point3) * axes_NumVertices, axes_points);
    glBufferSubData(GL_ARRAY_BUFFER,
        sizeof(point3) * axes_NumVertices,
        sizeof(color4) * axes_NumVertices,
        axes_colors);


    load();
    computeRotationVector();
    currRV = rotationVectors[seg];
    //std::cout << transVectors[0].x << transVectors[0].y << transVectors[0].z << std::endl;
    // Create and initialize a vertex buffer object for sphere
    glGenBuffers(1, &sphere_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, sphere_buffer);

    glBufferData(GL_ARRAY_BUFFER,
        sizeof(point3) * sphere_NumVertices + sizeof(color4) * sphere_NumVertices,
        NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0,
        sizeof(point3) * sphere_NumVertices, sphere_points);
    glBufferSubData(GL_ARRAY_BUFFER,
        sizeof(point3) * sphere_NumVertices,
        sizeof(color4) * sphere_NumVertices,
        sphere_colors);

    // Create and initialize a vertex buffer object for shadow
    glGenBuffers(1, &shadow_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, shadow_buffer);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(point3) * sphere_NumVertices + sizeof(color4) * sphere_NumVertices,
        NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0,
        sizeof(point3) * sphere_NumVertices, sphere_points);
    glBufferSubData(GL_ARRAY_BUFFER,
        sizeof(point3) * sphere_NumVertices,
        sizeof(color4) * sphere_NumVertices,
        sphere_shadow_colors);


    // Create and initialize a vertex buffer object for sphere w/ flat shding
    glGenBuffers(1, &sphere_flat_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, sphere_flat_buffer);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(point3) * sphere_NumVertices + sizeof(vec3) * sphere_NumVertices,
        NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0,
        sizeof(point3) * sphere_NumVertices, sphere_points);
    glBufferSubData(GL_ARRAY_BUFFER,
        sizeof(point3) * sphere_NumVertices,
        sizeof(vec3) * sphere_NumVertices,
        sphere_flat_normals);

    // Create and initialize a vertex buffer object for sphere w/ smooth shding, to be used in display()
    glGenBuffers(1, &sphere_smooth_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, sphere_smooth_buffer);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(point3) * sphere_NumVertices + sizeof(vec3) * sphere_NumVertices,
        NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0,
        sizeof(point3) * sphere_NumVertices, sphere_points);
    glBufferSubData(GL_ARRAY_BUFFER,
        sizeof(point3) * sphere_NumVertices,
        sizeof(vec3) * sphere_NumVertices,
        sphere_smooth_normals);


 // Load shaders and create a shader program (to be used in display())
    program = InitShader("vshader53.glsl", "fshader53.glsl");
    program2 = InitShader("vshaderP.glsl", "fshaderP.glsl");
    
    glEnable( GL_DEPTH_TEST ); /* Always enable z-buffer testing*/
    glClearColor(0.529f, 0.807f, 0.92f, 0.0f);    /* sky blue background*/
    glLineWidth(1.0);
    glPointSize(3.0);



}


//----------------------------------------------------------------------
// SetUp_Lighting_Uniform_Vars(mat4 mv):
// Set up lighting parameters that are uniform variables in shader.
//
// Note: "LightPosition" in shader must be in the Eye Frame.
//       So we use parameter "mv", the model-view matrix, to transform
//       light_position to the Eye Frame.
//----------------------------------------------------------------------

void SetUp_Lighting_Uniform_Vars(mat4 mv) {

    glUniform1i(glGetUniformLocation(program, "lightType"), lightType);

    glUniform1i(glGetUniformLocation(program, "fogType"), fogFlag);

        glUniform4fv(glGetUniformLocation(program, "DirectionalLightAmbient"),
            1, directional_light_ambient);
        glUniform4fv(glGetUniformLocation(program, "DirectionalLightDiffuse"),
            1, directional_light_diffuse);
        glUniform4fv(glGetUniformLocation(program, "DirectionalLightSpecular"),
            1, directional_light_specular);
        glUniform3fv(glGetUniformLocation(program, "LightDirection"),
            1, light_direction);

       // std::cout << "directional light" << std::endl;
    
        glUniform4fv(glGetUniformLocation(program, "PositionalLightAmbient"),
            1, positional_light_ambient);
        glUniform4fv(glGetUniformLocation(program, "PositionalLightDiffuse"),
            1, positional_light_diffuse);
        glUniform4fv(glGetUniformLocation(program, "PositionalLightSpecular"),
            1, positional_light_specular);
        // The Light Position in Eye Frame
        vec4 light_position_eyeFrame = mv * light_position;
        glUniform3fv(glGetUniformLocation(program, "LightPosition"),
            1, light_position_eyeFrame);

        glUniform1f(glGetUniformLocation(program, "ConstAtt"),
            const_att);
        glUniform1f(glGetUniformLocation(program, "LinearAtt"),
            linear_att);
        glUniform1f(glGetUniformLocation(program, "QuadAtt"),
            quad_att);
        /*spot light*/

        vec4 towards_eyeFrame = mv * towards;
        glUniform3fv(glGetUniformLocation(program, "Towards"),
            1, towards_eyeFrame);
        glUniform1f(glGetUniformLocation(program, "Exponent"),
            exponent);
        glUniform1f(glGetUniformLocation(program, "Cutoff"),
            cutoff);


}
void SetUp_Material(int mode) /*mode  0: floor; 1: sphere*/
{
    glUniform1f(glGetUniformLocation(program, "Shininess"),
        sphere_shininess);
    
    if (mode == 0) {
        glUniform4fv(glGetUniformLocation(program, "MaterialAmbient"),
            1, floor_ambient);
        glUniform4fv(glGetUniformLocation(program, "MaterialDiffuse"),
            1, floor_diffuse);
        glUniform4fv(glGetUniformLocation(program, "MaterialSpecular"),
            1, floor_specular);
       // std::cout << "color set for floor" << std::endl;
    }
    else if (mode ==1){
        glUniform4fv(glGetUniformLocation(program, "MaterialAmbient"),
            1, sphere_ambient);
        glUniform4fv(glGetUniformLocation(program, "MaterialDiffuse"),
            1, sphere_diffuse);
        glUniform4fv(glGetUniformLocation(program, "MaterialSpecular"),
            1, sphere_specular);
        

      //  std::cout << "color set for sphere" << std::endl;
    }


}




//----------------------------------------------------------------------------
// drawObj(buffer, num_vertices):
//   draw the object that is associated with the vertex buffer object "buffer"
//   and has "num_vertices" vertices.
//
void drawObj(GLuint buffer, int num_vertices, GLenum mode, int light = lightFlag, int tex = 0, int objectID = 0)
{
    //--- Activate the vertex buffer object to be drawn ---//
    glBindBuffer(GL_ARRAY_BUFFER, buffer);

    /*----- Set up vertex attribute arrays for each vertex attribute -----*/
    GLuint vPosition = glGetAttribLocation(program, "vPosition");
    glEnableVertexAttribArray(vPosition);
    glVertexAttribPointer(vPosition, 3, GL_FLOAT, GL_FALSE, 0,
			  BUFFER_OFFSET(0) );
    glUniform1i(glGetUniformLocation(program, "lighting"),light);


    glUniform1i(glGetUniformLocation(program, "texType"), tex);
    //std::cout << tex << std::endl;
    glUniform1i(glGetUniformLocation(program, "object"), objectID);

    GLuint vTexCoord;

    if (light ==0) {
        GLuint vColor = glGetAttribLocation(program, "vColor");
        glEnableVertexAttribArray(vColor);
        glVertexAttribPointer(vColor, 4, GL_FLOAT, GL_FALSE, 0,
            BUFFER_OFFSET(sizeof(point3) * num_vertices));
        // the offset is the (total) size of the previous vertex attribute array(s)
        /* Draw a sequence of geometric objs (triangles/lines) from the vertex buffer
       (using the attributes specified in each enabled vertex attribute array) */
        if ((tex != 0) && (objectID ==1)  ) {
            vTexCoord = glGetAttribLocation(program, "vTexCoord");
            glEnableVertexAttribArray(vTexCoord);
            glVertexAttribPointer(vTexCoord, 2, GL_FLOAT, GL_FALSE, 0,
                BUFFER_OFFSET((sizeof(point3)+ sizeof(color4)) * num_vertices));
        }

        glDrawArrays(mode, 0, num_vertices);
        glDisableVertexAttribArray(vColor);

    }
    else {

        GLuint vNormal = glGetAttribLocation(program, "vNormal");
        glEnableVertexAttribArray(vNormal);
        glVertexAttribPointer(vNormal, 3, GL_FLOAT, GL_FALSE, 0,
            BUFFER_OFFSET(sizeof(point3) * num_vertices));
        // the offset is the (total) size of the previous vertex attribute array(s)

        if ((tex != 0) && (objectID == 1)) {
            vTexCoord = glGetAttribLocation(program, "vTexCoord");
            glEnableVertexAttribArray(vTexCoord);
            glVertexAttribPointer(vTexCoord, 2, GL_FLOAT, GL_FALSE, 0,
                BUFFER_OFFSET((sizeof(point3) + sizeof(point3)) * num_vertices));
        }

        /* Draw a sequence of geometric objs (triangles) from the vertex buffer
           (using the attributes specified in each enabled vertex attribute array) */
        glDrawArrays(mode, 0, num_vertices);
        /*--- Disable each vertex attribute array being enabled ---*/
        glDisableVertexAttribArray(vNormal);
        

    }
    if ((tex != 0) && (objectID == 1)) {
        glDisableVertexAttribArray(vTexCoord);
    }

    /*--- Disable each vertex attribute array being enabled ---*/
    glDisableVertexAttribArray(vPosition);
        
}

void drawObj2(GLuint buffer) {
    //--- Activate the vertex buffer object to be drawn ---//
    glBindBuffer(GL_ARRAY_BUFFER, buffer);
    /*----- Set up vertex attribute arrays for each vertex attribute -----*/
    GLuint velocity = glGetAttribLocation(program2, "v");
    glEnableVertexAttribArray(velocity);
    glVertexAttribPointer(velocity, 3, GL_FLOAT, GL_FALSE, 0,
        BUFFER_OFFSET(0));
    GLuint vColor = glGetAttribLocation(program2, "vColor");
    glEnableVertexAttribArray(vColor);
    glVertexAttribPointer(vColor, 4, GL_FLOAT, GL_FALSE, 0,
        BUFFER_OFFSET(sizeof(vec3) * NumParticles));

    glDrawArrays(GL_POINTS,0,NumParticles); //%^&*

    /*--- Disable each vertex attribute array being enabled ---*/
    glDisableVertexAttribArray(velocity);
    glDisableVertexAttribArray(vColor);
};

void drawShadow(mat4 view,GLuint model_view) {
    if (shadowFlag == 1 && eye.y >= 0) { //draw the shadow

        mat4 mv = view * shadow_proj * Translate(center.x, center.y, center.z) * accumulatedRM;
        glUniformMatrix4fv(model_view, 1, GL_TRUE, mv);

        if (sphereFlag == 1) { // Filled shadow
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }
        else {              // Wireframe shadow
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        }
        drawObj(shadow_buffer, sphere_NumVertices, GL_TRIANGLES, 0, 0, 3);  // draw the shadow
    }
}
//----------------------------------------------------------------------------
void display( void )
{
  GLuint  model_view;  // model-view matrix uniform shader variable location
  GLuint  projection;  // projection matrix uniform shader variable location

    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glUseProgram(program); // Use the shader program

    model_view = glGetUniformLocation(program, "ModelView" );
    projection = glGetUniformLocation(program, "Projection" );

    // Set the value of the fragment shader texture sampler variable
  //   ("texture_2D") to the appropriate texture unit. In this case,
  //   0, for GL_TEXTURE0 which was previously set in init() by calling
  //   glActiveTexture( GL_TEXTURE0 ).
    glUniform1i(glGetUniformLocation(program, "texture_2D"), 0);
    glUniform1i(glGetUniformLocation(program, "texture_1D"), 1);

    glUniform1i(glGetUniformLocation(program, "if_object"), objFlag);
    glUniform1i(glGetUniformLocation(program, "vert"), verticalFlag);
    glUniform1i(glGetUniformLocation(program, "slanted"), slantedFlag);

    glUniform1i(glGetUniformLocation(program, "lattice"), latticeFlag);
    glUniform1i(glGetUniformLocation(program, "upright"), uprightFlag);
    glUniform1i(glGetUniformLocation(program, "tilted"), tiltedFlag);


    glClearColor(0.529f, 0.807f, 0.92f, 0.0f);    /* sky blue background*/

/*---  Set up and pass on Projection matrix to the shader ---*/
    mat4  p = Perspective(fovy, aspect, zNear, zFar);
    glUniformMatrix4fv(projection, 1, GL_TRUE, p); // GL_TRUE: matrix is row-major

/*---  Set up and pass on Model-View matrix to the shader ---*/
    // eye is a global variable of vec4 set to init_eye and updated by keyboard()
    //at = VRP + VPN
    vec4    at(0.0, 0.0, 0.0, 1.0);
    vec4    up(0.0, 1.0, 0.0, 0.0);
    mat4  mv = LookAt(eye, at, up);

    glUniformMatrix4fv(model_view, 1, GL_TRUE, mv); // GL_TRUE: matrix is row-major
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    drawObj(axes_buffer, axes_NumVertices,GL_LINES,0,0,0); // draw the axes
    
    if (lightFlag == 1) {
        SetUp_Lighting_Uniform_Vars(mv);
    }

/*----- Set Up the Model-View matrix for the sphere -----*/
    accumulatedRM = Rotate(theta,currRV.x,currRV.y,currRV.z)*accumulatedRM;
    mv = mv* Translate(center.x, center.y, center.z)*accumulatedRM;
    glUniformMatrix4fv(model_view, 1, GL_TRUE, mv); // GL_TRUE: matrix is row-major
    if (sphereFlag == 1) { // Filled sphere
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        if (lightFlag == 0) {
            drawObj(sphere_buffer, sphere_NumVertices, GL_TRIANGLES, 0, sphereTexFlag,2);  // draw the sphere without light
        }
        else {
            SetUp_Material(1);
            mat3 normal_matrix = NormalMatrix(mv, 1);
            glUniformMatrix3fv(glGetUniformLocation(program, "Normal_Matrix"),
                1, GL_TRUE, normal_matrix);
            if (shadingFlag == 1) { /* Flat shading*/
                drawObj(sphere_flat_buffer, sphere_NumVertices, GL_TRIANGLES, 1, sphereTexFlag,2);  
            }
            else { /*smooth Shading*/
                drawObj(sphere_smooth_buffer, sphere_NumVertices, GL_TRIANGLES, 1, sphereTexFlag,2); 
            }
           
        }
    }
    else {              // Wireframe sphere
        
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        drawObj(sphere_buffer, sphere_NumVertices, GL_TRIANGLES,0, 0,2);  // draw the Wireframe sphere
        
    }
    /*----- Set Up the Model-View matrix for the floor and shadow -----*/

    

    glDepthMask(GL_FALSE); /*Disable writing to Z-buffer*/

    mv = LookAt(eye, at, up);

    glUniformMatrix4fv(model_view, 1, GL_TRUE, mv); // GL_TRUE: matrix is row-major
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    
    /*draw the floor to frame buffer*/
    if (lightFlag == 1) {
        SetUp_Material(0);
        mat3 normal_matrix = NormalMatrix(mv, 1);
        glUniformMatrix3fv(glGetUniformLocation(program, "Normal_Matrix"),
            1, GL_TRUE, normal_matrix);

        drawObj(floor_shading_buffer, floor_NumVertices, GL_TRIANGLES, 1,floorTexFlag,1);  // draw the floor with light
      //  std::cout << "shading floor" << std::endl;

    }
    else {
        drawObj(floor_buffer, floor_NumVertices, GL_TRIANGLES, 0, floorTexFlag,1);  
    }
    mv = LookAt(eye, at, up);
    if (blendFlag == 0) {
        glDepthMask(GL_TRUE); 
        drawShadow(mv,model_view); //Draw shadow to both buffers
    }
    else {
        //Enable shadow blending when selected
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        drawShadow(mv, model_view); //Draw shadow to frame buffer only
        glDisable(GL_BLEND);
        glDepthMask(GL_TRUE); //Enable writing to Z-buffer
    }
      
    glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE); /* Disable writing to frame buffer*/
    mv = LookAt(eye, at, up);
    glUniformMatrix4fv(model_view, 1, GL_TRUE, mv); // GL_TRUE: matrix is row-major
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    /*draw the floor to z - buffer*/
    if (lightFlag == 1) {
        drawObj(floor_shading_buffer, floor_NumVertices, GL_TRIANGLES, 1, floorTexFlag,1);  // draw the floor with light
    }
    else {
        drawObj(floor_buffer, floor_NumVertices, GL_TRIANGLES, 0, floorTexFlag,1);  // draw the floor to frame buffer
    }

    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE); /*Enable writing to frame buffer*/


    glUseProgram(program2); // Use the shader program for particles
    model_view = glGetUniformLocation(program2, "ModelView");
    projection = glGetUniformLocation(program2, "Projection");

    glUniformMatrix4fv(projection, 1, GL_TRUE, p); // GL_TRUE: matrix is row-major
    mv = LookAt(eye, at, up);
    glUniformMatrix4fv(model_view, 1, GL_TRUE, mv); // GL_TRUE: matrix is row-major
    
    
    if (beginFlag == 0) {
        time_elapsed = 0.0f;
    }
    else {
        if (reset_time == 1) {
            time_start = glutGet(GLUT_ELAPSED_TIME);
            reset_time = 0;
        }
        time_elapsed = glutGet(GLUT_ELAPSED_TIME) - time_start;
        if (time_elapsed > T_max) {
            time_elapsed = 0.0f;
            time_start = glutGet(GLUT_ELAPSED_TIME);
        }
    }
    glUniform1f(glGetUniformLocation(program2, "t"), time_elapsed);

    if (firework == 1) { 
        glPolygonMode(GL_FRONT_AND_BACK, GL_POINT); 
        drawObj2(particle_buffer); 
        //std::cout << "firework" << std::endl;
    }
    glutSwapBuffers();
   // std::cout << "display finish" << std::endl;
}
//---------------------------------------------------------------------------
void idle (void)
{
      angle += theta;  //YJC: change this value to adjust the sphere rotation speed, in degrees

      GLfloat travelled = theta / 180.0f * PI * radius;
      center.x += transVectors[seg].x * travelled;
      center.y += transVectors[seg].y * travelled;
      center.z += transVectors[seg].z * travelled;

      if (distance(center,origin) > distance(turningPts[seg+1],origin)) { //proceed with next segment
          seg = (seg + 1) % 3;
          origin = turningPts[seg];
          currRV = rotationVectors[seg];
      }
    glutPostRedisplay();
}
//----------------------------------------------------------------------------
void keyboard(unsigned char key, int x, int y)
{
    switch(key) {
	case 033: // Escape Key
	case 'q': case 'Q':
	    exit( EXIT_SUCCESS );
	    break;

        case 'X': eye[0] += 1.0; break;
	case 'x': eye[0] -= 1.0; break;
        case 'Y': eye[1] += 1.0; break;
	case 'y': eye[1] -= 1.0; break;
        case 'Z': eye[2] += 1.0; break;
	case 'z': eye[2] -= 1.0; break;

    case 'B': case 'b': //start animation
            if (beginFlag == 0) reset_time = 1;
            beginFlag = 1;
            animationFlag = 1;
            glutIdleFunc(idle);
            break;
    case 'V': case 'v': // vertical texture
        verticalFlag = 1;
        slantedFlag = 0;
        break;
    case 'S': case 's': //slanted texture
        verticalFlag = 0;
        slantedFlag = 1;
        break;
    case 'O': case 'o': //toggles btw eye/object frames
        objFlag = 1 - objFlag;
        break;
    case 'U': case 'u': //upright lattice
        uprightFlag = 1;
        tiltedFlag = 0;
        break;
    case 'T': case 't': //tilted lattice
        uprightFlag = 0;
        tiltedFlag = 1;
        break;
    case 'L': case 'l': //toggle to enable/disable lattice
        latticeFlag = 1 - latticeFlag;
        break;

	case ' ':  // reset to initial viewer/eye position
	    eye = init_eye;
	    break;
    }

    //std::cout << eye[0]<<eye[1]<<eye[2] << std::endl;
    glutPostRedisplay();
}

void mouse(int button, int state, int x, int y) {
    if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN && beginFlag) {
        animationFlag = 1 - animationFlag;}
    if (animationFlag == 1) glutIdleFunc(idle);
    else                    glutIdleFunc(NULL);

}

// functions to handle menu action
void mainMenu(int id) {

    switch (id) {
    case 1:
        eye = init_eye;
        break;
    case 2:
        sphereFlag = 0;
        break;
    case 3:
        exit(EXIT_SUCCESS);
        break;
    }
    glutPostRedisplay();
}

void shadowMenu(int id) {
    shadowFlag = id % 2;
    glutPostRedisplay();
}

void lightingMenu(int id) {
    lightFlag = id % 2;
    glutPostRedisplay();
}
/*
void sphereMenu(int id) {
    sphereFlag = id % 2;
    glutPostRedisplay();
}
*/
void shadingMenu(int id) {

    sphereFlag = 1; //draws a solid sphere 
    shadingFlag = id % 2;
    glutPostRedisplay();
}

void sourceMenu(int id) {
    lightType = id;
    glutPostRedisplay();
}

void fogMenu(int id) {
    fogFlag = id % 4;
    glutPostRedisplay();
}

void blendMenu(int id) {
    blendFlag = id % 2;
    glutPostRedisplay();
}

void floorTexMenu(int id) {
    floorTexFlag = id;
    glutPostRedisplay();
}

void sphereTexMenu(int id) {
    sphereTexFlag = id % 3;
    glutPostRedisplay();
}

void fireworkMenu(int id) {
    if (firework != id % 2) { // status change occurs
        if (firework == 0) { // fresh start
            reset_time = 1;
            firework = 1;
        }
        else {
            firework = 0;
        }
    }
}

void createMenu() {
    
    int shadowMenuID = glutCreateMenu(shadowMenu);
    glutAddMenuEntry(" Yes ",1);
    glutAddMenuEntry(" No ",2);
    glutSetMenuFont(shadowMenuID, GLUT_BITMAP_HELVETICA_18);
    int lightingMenuID = glutCreateMenu(lightingMenu);
    glutAddMenuEntry(" Yes ", 1);
    glutAddMenuEntry(" No ", 2);
    glutSetMenuFont(lightingMenuID, GLUT_BITMAP_HELVETICA_18);
    /*
    int sphereMenuID = glutCreateMenu(sphereMenu);
    glutAddMenuEntry(" Yes ", 2);
    glutAddMenuEntry(" No ", 1);
    */
    int shadingMenuID = glutCreateMenu(shadingMenu);
    glutAddMenuEntry(" Flat Shading ", 1);
    glutAddMenuEntry(" Smooth Shading ", 2);
    glutSetMenuFont(shadingMenuID, GLUT_BITMAP_HELVETICA_18);
    int sourceMenuID = glutCreateMenu(sourceMenu);
    glutAddMenuEntry(" Spot Light ", 2);
    glutAddMenuEntry(" Point Source ", 1);
    glutSetMenuFont(sourceMenuID, GLUT_BITMAP_HELVETICA_18);

    int fogMenuID = glutCreateMenu(fogMenu);
    glutAddMenuEntry(" No Fog ", 4);
    glutAddMenuEntry(" Linear ", 1);
    glutAddMenuEntry(" Exponential ", 2);
    glutAddMenuEntry(" Exponential Square ", 3);
    glutSetMenuFont(fogMenuID, GLUT_BITMAP_HELVETICA_18);

    int blendMenuID = glutCreateMenu(blendMenu);
    glutAddMenuEntry(" No ", 2);
    glutAddMenuEntry(" Yes ", 1);
    glutSetMenuFont(blendMenuID, GLUT_BITMAP_HELVETICA_18);

    int floorTexMenuID = glutCreateMenu(floorTexMenu);
    glutAddMenuEntry(" No ", 0);
    glutAddMenuEntry(" Yes ", 2);
    glutSetMenuFont(floorTexMenuID, GLUT_BITMAP_HELVETICA_18);

    int sphereTexMenuID = glutCreateMenu(sphereTexMenu);
    glutAddMenuEntry(" No ", 3);
    glutAddMenuEntry(" Yes - Contour Lines ", 1);
    glutAddMenuEntry(" Yes - Checkerboard ", 2);
    glutSetMenuFont(sphereTexMenuID, GLUT_BITMAP_HELVETICA_18);

    int fireworkMenuID = glutCreateMenu(fireworkMenu);
    glutAddMenuEntry(" Yes ", 1);
    glutAddMenuEntry(" No ", 2);
    glutSetMenuFont(fireworkMenuID, GLUT_BITMAP_HELVETICA_18);
    
    int menuID = glutCreateMenu(mainMenu);
    glutAddMenuEntry(" Default View Point ", 1);
    glutAddSubMenu(" Shadow ",shadowMenuID);
    glutAddSubMenu(" Enable Lighting ",lightingMenuID);
    glutAddMenuEntry(" Wire Frame Sphere ", 2);
   // glutAddSubMenu(" Wire Frame Sphere ",sphereMenuID);
    glutAddSubMenu(" Shading ",shadingMenuID);
    glutAddSubMenu(" Light Source ", sourceMenuID);
    glutAddSubMenu(" Fog Options ", fogMenuID);
    glutAddSubMenu(" Blending Shadow ", blendMenuID);
    glutAddSubMenu(" Texture Mapped Ground ", floorTexMenuID);
    glutAddSubMenu(" Texture Mapped Sphere ", sphereTexMenuID);
    glutAddSubMenu(" Firework ", fireworkMenuID);
    glutAddMenuEntry(" Quit ", 3);
    glutSetMenuFont(menuID, GLUT_BITMAP_HELVETICA_18);
    glutAttachMenu(GLUT_LEFT_BUTTON);

}

//----------------------------------------------------------------------------
void reshape(int width, int height)
{
    glViewport(0, 0, width, height);
    aspect = (GLfloat) width  / (GLfloat) height;
    glutPostRedisplay();
}
//----------------------------------------------------------------------------
int main( int argc, char **argv )
{
    glutInit(&argc, argv);
#ifdef __APPLE__ // Enable core profile of OpenGL 3.2 on macOS.
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_3_2_CORE_PROFILE);
#else
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
#endif
    glutInitWindowSize(512, 512);
    glutCreateWindow("Rolling Sphere");
    createMenu();
#ifdef __APPLE__ // on macOS
    // Core profile requires to create a Vertex Array Object (VAO).
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
#else           // on Linux or Windows, we still need glew
    /* Call glewInit() and error checking */
    int err = glewInit();
    if (GLEW_OK != err)
    { 
        printf("Error: glewInit failed: %s\n", (char*) glewGetErrorString(err)); 
        exit(1);
    }
#endif

    // Get info of GPU and supported OpenGL version
    printf("Renderer: %s\n", glGetString(GL_RENDERER));
    printf("OpenGL version supported %s\n", glGetString(GL_VERSION));

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(NULL);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    init();
    glutMainLoop();
    return 0;
}
