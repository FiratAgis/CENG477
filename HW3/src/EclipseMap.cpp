#include "EclipseMap.h"

using namespace std;

struct vertex {
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec2 texture;

    vertex() {}

    vertex(const glm::vec3 &position, const glm::vec3 &normal, const glm::vec2 &texture) : position(position),
                                                                                           normal(normal),
                                                                                           texture(texture) {}
};

struct triangle {
    int vertex1;
    int vertex2;
    int vertex3;

    triangle() {}

    triangle(const int &vertex1, const int &vertex2, const int &vertex3) : vertex1(vertex1), vertex2(vertex2),
                                                                           vertex3(vertex3) {}
};

//Added
glm::vec3 getStepCoordinates(float size, int horizontal_step, int horizontal_split_count, int vertical_step, int vertical_split_count)
{
    float alpha, beta;
    float x, y, z;
    alpha = 2 * PI * ((float) horizontal_step / (float) horizontal_split_count);
    beta = PI * ((float) vertical_step / (float) vertical_split_count);
    z = size * cos(beta);
    y = size * sin(beta) * sin(alpha);
    x = size * sin(beta) * cos(alpha);
    return glm::vec3(x, y, z);
}

glm::vec3 getStepCoordinates(glm::vec3 center, float size, int horizontal_step, int horizontal_split_count, int vertical_step, int vertical_split_count)
{
    glm::vec3 stepPos = getStepCoordinates(size, horizontal_step, horizontal_split_count, vertical_step, vertical_split_count);
    return stepPos + center;
}

glm::vec3 getNormalOfStepCoordinates(glm::vec3 stepPos)
{
    return glm::normalize(stepPos);
}

glm::vec3 getNormalOfStepCoordinates(glm::vec3 stepPos, glm::vec3 center)
{
    return glm::normalize(stepPos-center);
}

void setupVertices(vector<float> *vertices, vector<unsigned int> *indices, float size, int horizontal_split_count, int vertical_split_count)
{
    for(int i = 0; i <= horizontal_split_count; i++)
    {
        for(int j = 0; j <= vertical_split_count; j++)
        {
            glm::vec3 stepPos = getStepCoordinates(size, i, horizontal_split_count, j, vertical_split_count);
            vertex temp = vertex(stepPos, getNormalOfStepCoordinates(stepPos), glm::vec2((float) i / horizontal_split_count, (float) j / vertical_split_count));
            vertices->push_back(temp.position.x);
            vertices->push_back(temp.position.y);
            vertices->push_back(temp.position.z);
            vertices->push_back(temp.normal.x);
            vertices->push_back(temp.normal.y);
            vertices->push_back(temp.normal.z);
            vertices->push_back(temp.texture.x);
            vertices->push_back(temp.texture.y);
        }
    }

    for(int i = 0; i < vertical_split_count; i++)
    {
        for(int j = 0; j < horizontal_split_count; j++)
        {
            indices->push_back(j + i * (horizontal_split_count + 1));
            indices->push_back(j + (i + 1) * (horizontal_split_count + 1));
            indices->push_back(j + i * (horizontal_split_count + 1) + 1);

            indices->push_back(j + i * (horizontal_split_count + 1) + 1);
            indices->push_back(j + (i + 1) * (horizontal_split_count + 1));
            indices->push_back(j + (i + 1) * (horizontal_split_count + 1) + 1);

        }
    }
}

void EclipseMap::bindBuffer(unsigned int vbo, unsigned int ebo, unsigned int normal, unsigned int text)
{
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GL_FLOAT), 0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, normal);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GL_FLOAT), 0);
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, text);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(GL_FLOAT), 0);
    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
}


void EclipseMap::cameraPositionSet() {
    cameraPosition = cameraPosition + cameraSpeed *  cameraDirection ;
    center = cameraPosition + cameraDirection * near;
    view = glm::lookAt(cameraPosition, center, cameraUp);    
    
    normal = glm::inverseTranspose(view);
    
    projection = glm::perspective(projectionAngle, aspectRatio, near, far);

    mvpMoon = projection * view * modelMoon;
    mvpWorld = projection * view * modelWorld;
}


//to parse cam and matrix to shader
void EclipseMap::cameraSetGL(GLuint ShaderID) {
    GLint cam_loc = glGetUniformLocation(ShaderID, "cameraPosition");
    glUniform3fv(cam_loc, 1, &cameraPosition[0]);
    
    GLint projection_loc = glGetUniformLocation(ShaderID, "ProjectionMatrix");
    glUniformMatrix4fv(projection_loc, 1, GL_FALSE, &projection[0][0]);

    GLint view_loc = glGetUniformLocation(ShaderID, "ViewMatrix");
    glUniformMatrix4fv(view_loc, 1, GL_FALSE, &view[0][0]);

    GLint normal_loc = glGetUniformLocation(ShaderID, "NormalMatrix");
    glUniformMatrix4fv(normal_loc, 1, GL_FALSE, &normal[0][0]);
}

void EclipseMap::setOthers(GLuint ShaderID, bool is_moon) {
    GLint light_loc = glGetUniformLocation(ShaderID, "lightPosition");
    glUniform3fv(light_loc, 1, &lightPos[0]);


    GLint offset_loc = glGetUniformLocation(ShaderID, "textureOffset");
    glUniform1f(offset_loc, textureOffset);
    if (is_moon) {
        GLint width_loc = glGetUniformLocation(ShaderID, "imageWidth");
        glUniform1f(width_loc, moonImageWidth);

        GLint imageHeight_loc = glGetUniformLocation(ShaderID, "imageHeight");
        glUniform1f(imageHeight_loc, moonImageHeight);

        GLint orbit_loc = glGetUniformLocation(ShaderID, "orbitDegree");
        glUniform1f(orbit_loc, orbitDegree);

        GLint mvp_loc = glGetUniformLocation(ShaderID, "MVP");
        glUniformMatrix4fv(mvp_loc, 1, GL_FALSE, &mvpMoon[0][0]);

    }
    else {
        GLint height_loc = glGetUniformLocation(ShaderID, "heightFactor");
        glUniform1f(height_loc, heightFactor);

        GLint width_loc = glGetUniformLocation(ShaderID, "imageWidth");
        glUniform1f(width_loc, imageWidth);

        GLint imageHeight_loc = glGetUniformLocation(ShaderID, "imageHeight");
        glUniform1f(imageHeight_loc, imageHeight);

        GLint mvp_loc = glGetUniformLocation(ShaderID, "MVP");
        glUniformMatrix4fv(mvp_loc, 1, GL_FALSE, &mvpWorld[0][0]);
    }
}



void EclipseMap::Render(const char *coloredTexturePath, const char *greyTexturePath, const char *moonTexturePath) {
    // Open window
    GLFWwindow *window = openWindow(windowName, screenWidth, screenHeight);

    // Moon commands
    // Load shaders
    GLuint moonShaderID = initShaders("moonShader.vert", "moonShader.frag");

    initMoonColoredTexture(moonTexturePath, moonShaderID);

    // TODO: Set moonVertices
    
    setupVertices(&moonVertices, &moonIndices, moonRadius, horizontalSplitCount, verticalSplitCount);
    
    
    // TODO: Configure Buffers

    glGenVertexArrays(1, &moonVAO);
	glBindVertexArray(moonVAO);

    moonVBO = 0;
    moonNormal = 0;
    moonText = 0;
    glGenBuffers(1, &moonVBO);
    glGenBuffers(1, &moonNormal);
    glGenBuffers(1, &moonText);
    glBindBuffer(GL_ARRAY_BUFFER, moonVBO);
    glBufferData(GL_ARRAY_BUFFER, moonVertices.size() * sizeof(GL_FLOAT), &moonVertices[0], GL_STATIC_DRAW);
    
    glBindBuffer(GL_ARRAY_BUFFER, moonNormal);
    glBufferData(GL_ARRAY_BUFFER, (moonVertices.size() - 3) * sizeof(GL_FLOAT), &moonVertices[3], GL_STATIC_DRAW);
    
    glBindBuffer(GL_ARRAY_BUFFER, moonText);
    glBufferData(GL_ARRAY_BUFFER, (moonVertices.size() - 6) * sizeof(GL_FLOAT), &moonVertices[6], GL_STATIC_DRAW);

    moonEBO = 0;
    glGenVertexArrays(1, &moonEBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, moonEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, moonIndices.size() * sizeof(unsigned int), &moonIndices[0], GL_STATIC_DRAW);

    


    // World commands
    // Load shaders
    GLuint worldShaderID = initShaders("worldShader.vert", "worldShader.frag");

    initColoredTexture(coloredTexturePath, worldShaderID);
    initGreyTexture(greyTexturePath, worldShaderID);

    // TODO: Set worldVertices
    
    setupVertices(&worldVertices, &worldIndices, radius, horizontalSplitCount, verticalSplitCount);


    // TODO: Configure Buffers

    worldVBO = 0;
    worldNormal = 0;
    worldText = 0;
    glGenBuffers(1, &worldVBO);
    glGenBuffers(1, &worldNormal);
    glGenBuffers(1, &worldText);
    glBindBuffer(GL_ARRAY_BUFFER, worldVBO);
    glBufferData(GL_ARRAY_BUFFER, worldVertices.size() * sizeof(GL_FLOAT), &worldVertices[0], GL_STATIC_DRAW);
    
    glBindBuffer(GL_ARRAY_BUFFER, worldNormal);
    glBufferData(GL_ARRAY_BUFFER, (worldVertices.size() - 3) * sizeof(GL_FLOAT), &worldVertices[3], GL_STATIC_DRAW);
    
    glBindBuffer(GL_ARRAY_BUFFER, worldText);
    glBufferData(GL_ARRAY_BUFFER, (worldVertices.size() - 6) * sizeof(GL_FLOAT), &worldVertices[6], GL_STATIC_DRAW);

    worldEBO = 0;
    glGenVertexArrays(1, &worldEBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, worldEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, worldIndices.size() * sizeof(unsigned int), &worldIndices[0], GL_STATIC_DRAW);
    
    // Enable depth test
    glEnable(GL_DEPTH_TEST);
    // Main rendering loop

    int current_width, current_height;
    do {
        glfwGetWindowSize(window, &current_width, &current_height);
        screenWidth = current_width;
        screenHeight = current_height;
        aspectRatio = current_width / current_height;
        glViewport(0, 0, screenWidth, screenHeight);

        glClearStencil(0);
        glClearDepth(1.0f);
        glClearColor(0, 0, 0, 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);


        // TODO: Handle key presses
        handleKeyPress(window);

        // TODO: Manipulate rotation variables
        modelMoon = glm::rotate(0.02f, glm::vec3(0.f, 0.0f, 1.0f)) * modelMoon * glm::rotate(0.5f/horizontalSplitCount, glm::vec3(0.f, 0.0f, 1.0f));
        modelWorld = modelWorld *glm::rotate(0.5f/horizontalSplitCount, glm::vec3(0.f, 0.0f, 1.0f));

        //Arranges camera, set's shader uniform variables except texture variables
        
        
        // TODO: Bind textures

        // TODO: Use moonShaderID program
        glUseProgram(moonShaderID);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D,moonTextureColor);
        
        cameraPositionSet();

        // TODO: Update camera at every frame
        
        // TODO: Update uniform variables at every frame
        cameraSetGL(moonShaderID);
        setOthers(moonShaderID, true);
        

        // TODO: Bind moon vertex array
        bindBuffer(moonVBO, moonEBO, moonNormal, moonText);

        // TODO: Draw moon object
        glDrawElements(GL_TRIANGLES, moonIndices.size(), GL_UNSIGNED_INT, 0);
        
        /*************************/
        
        // TODO: Use worldShaderID program
        glUseProgram(worldShaderID);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D,textureColor);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D,textureGrey);


        // TODO: Update camera at every frame
        
        

        cameraSetGL(worldShaderID);
        setOthers(worldShaderID, false);

        // TODO: Update uniform variables at every frame

        // TODO: Bind world vertex array
        bindBuffer(worldVBO, worldEBO, worldNormal, worldText);

        glDrawElements(GL_TRIANGLES, worldIndices.size(), GL_UNSIGNED_INT, 0);

        // Swap buffers and poll events
        
        glfwSwapBuffers(window);
        glfwPollEvents();
    } while (!glfwWindowShouldClose(window));

    // Delete buffers
    glDeleteBuffers(1, &moonVAO);
    glDeleteBuffers(1, &moonVBO);
    glDeleteBuffers(1, &moonEBO);
    glDeleteBuffers(1, &moonNormal);
    glDeleteBuffers(1, &moonText);


    // Delete buffers
    glDeleteBuffers(1, &worldVAO);
    glDeleteBuffers(1, &worldVBO);
    glDeleteBuffers(1, &worldEBO);
    glDeleteBuffers(1, &worldNormal);
    glDeleteBuffers(1, &worldText);

    glDeleteProgram(moonShaderID);
    glDeleteProgram(worldShaderID);

    // Close window
    glfwTerminate();
}


bool is_pressed = false;

void EclipseMap::handleKeyPress(GLFWwindow *window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
    else if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {
        //increase height factor by 10
        heightFactor += 10;
   
    }
    else if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS) {
        //decrease height factor by 10
        heightFactor -= 10;

    }
    else if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        //rotate gaze around up vector by 0.05 unit
        glm::mat4 rotator = glm::rotate(0.05f, cameraUp);
        cameraDirection = glm::vec3(rotator * glm::vec4(cameraDirection, 1.0f)) ;
    }
    else if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        //rotate gaze around up vector by 0.05 unit
        glm::mat4 rotator = glm::rotate(-0.05f, cameraUp);
        cameraDirection = glm::vec3(rotator * glm::vec4(cameraDirection, 1.0f)) ;

    }
    else if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        //rotate gaze around left vector by 0.05 unit
        glm::mat4 rotator = glm::rotate(-0.05f, cameraLeft);
        cameraDirection = glm::vec3(rotator * glm::vec4(cameraDirection, 1.0f)) ;

  
    }
    else if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        //rotate gaze around left vector by 0.05 unit
        glm::mat4 rotator = glm::rotate(0.05f, cameraLeft);
        cameraDirection = glm::vec3(rotator * glm::vec4(cameraDirection, 1.0f)) ;
    }
    else if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS) {
        //increase speed by 0.01   
        cameraSpeed += 0.01;
    }
    else if (glfwGetKey(window, GLFW_KEY_H) == GLFW_PRESS) {
        //decrease speed by 0.01
        cameraSpeed -= 0.01;
    }
    else if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS) {
        //plane will stop, speed = 0
        cameraSpeed = 0.0;   
    }
    else if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS) {
        //plane wil get to initial state, I   (0,4k, 4k)
        cameraPosition = cameraStartPosition;
        cameraDirection = cameraStartDirection;
        cameraSpeed = cameraStartSpeed;
    }
    else if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {
        //fullscreen
        if (!is_pressed) {
            pKeyPressed = ! pKeyPressed;
            GLFWmonitor* monitor = glfwGetPrimaryMonitor();
            const GLFWvidmode *mode = glfwGetVideoMode(monitor);

            if ( pKeyPressed ) {
                int current_height, current_width;

                glfwSetWindowMonitor(window, monitor,  0, 0, mode->width, mode->height, mode->refreshRate);

            }
            else {
                glfwSetWindowMonitor(window, nullptr, 1, 31, defaultScreenWidth, defaultScreenHeight, mode->refreshRate);
            }
        is_pressed = true;
        }

    }
    else if (glfwGetKey(window, GLFW_KEY_P) == GLFW_RELEASE) {
        is_pressed = false;
    }
    

}

GLFWwindow *EclipseMap::openWindow(const char *windowName, int width, int height) {
    if (!glfwInit()) {
        getchar();
        return 0;
    }

    const GLFWvidmode *mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GLFWwindow *window = glfwCreateWindow(width, height, windowName, NULL, NULL);
    glfwSetWindowMonitor(window, NULL, 1, 31, screenWidth, screenHeight, mode->refreshRate);

    if (window == NULL) {
        getchar();
        glfwTerminate();
        return 0;
    }

    glfwMakeContextCurrent(window);

    glewExperimental = true;
    if (glewInit() != GLEW_OK) {
        getchar();
        glfwTerminate();
        return 0;
    }

    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    glClearColor(0, 0, 0, 0);

    return window;
}

void EclipseMap::initColoredTexture(const char *filename, GLuint shader) {
    int width, height;
    glGenTextures(1, &textureColor);
    cout << shader << endl;
    glBindTexture(GL_TEXTURE_2D, textureColor);
    // set the texture wrapping parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                    GL_CLAMP_TO_EDGE);    // set texture wrapping to GL_REPEAT (default wrapping method)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    // set texture filtering parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    unsigned char *raw_image = NULL;
    int bytes_per_pixel = 3;   /* or 1 for GRACYSCALE images */
    int color_space = JCS_RGB; /* or JCS_GRAYSCALE for grayscale images */

    /* these are standard libjpeg structures for reading(decompression) */
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    /* libjpeg data structure for storing one row, that is, scanline of an image */
    JSAMPROW row_pointer[1];

    FILE *infile = fopen(filename, "rb");
    unsigned long location = 0;
    int i = 0, j = 0;

    if (!infile) {
        printf("Error opening jpeg file %s\n!", filename);
        return;
    }
    printf("Texture filename = %s\n", filename);

    /* here we set up the standard libjpeg error handler */
    cinfo.err = jpeg_std_error(&jerr);
    /* setup decompression process and source, then read JPEG header */
    jpeg_create_decompress(&cinfo);
    /* this makes the library read from infile */
    jpeg_stdio_src(&cinfo, infile);
    /* reading the image header which contains image information */
    jpeg_read_header(&cinfo, TRUE);
    /* Start decompression jpeg here */
    jpeg_start_decompress(&cinfo);

    /* allocate memory to hold the uncompressed image */
    raw_image = (unsigned char *) malloc(cinfo.output_width * cinfo.output_height * cinfo.num_components);
    /* now actually read the jpeg into the raw buffer */
    row_pointer[0] = (unsigned char *) malloc(cinfo.output_width * cinfo.num_components);
    /* read one scan line at a time */
    while (cinfo.output_scanline < cinfo.image_height) {
        jpeg_read_scanlines(&cinfo, row_pointer, 1);
        for (i = 0; i < cinfo.image_width * cinfo.num_components; i++)
            raw_image[location++] = row_pointer[0][i];
    }

    height = cinfo.image_height;
    width = cinfo.image_width;


    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, raw_image);
   

    imageWidth = width;
    imageHeight = height;

    glGenerateMipmap(GL_TEXTURE_2D);

    glUseProgram(shader); // don't forget to activate/use the shader before setting uniforms!
    // either set it manually like so:

    glUniform1i(glGetUniformLocation(shader, "TexColor"), 0);
    /* wrap up decompression, destroy objects, free pointers and close open files */
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    free(row_pointer[0]);
    free(raw_image);
    fclose(infile);

}

void EclipseMap::initGreyTexture(const char *filename, GLuint shader) {

    glGenTextures(1, &textureGrey);
    glBindTexture(GL_TEXTURE_2D, textureGrey);
    // set the texture wrapping parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                    GL_CLAMP_TO_EDGE);    // set texture wrapping to GL_REPEAT (default wrapping method)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    // set texture filtering parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    int width, height;

    unsigned char *raw_image = NULL;
    int bytes_per_pixel = 3;   /* or 1 for GRACYSCALE images */
    int color_space = JCS_RGB; /* or JCS_GRAYSCALE for grayscale images */

    /* these are standard libjpeg structures for reading(decompression) */
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    /* libjpeg data structure for storing one row, that is, scanline of an image */
    JSAMPROW row_pointer[1];

    FILE *infile = fopen(filename, "rb");
    unsigned long location = 0;
    int i = 0, j = 0;

    if (!infile) {
        printf("Error opening jpeg file %s\n!", filename);
        return;
    }
    printf("Texture filename = %s\n", filename);

    /* here we set up the standard libjpeg error handler */
    cinfo.err = jpeg_std_error(&jerr);
    /* setup decompression process and source, then read JPEG header */
    jpeg_create_decompress(&cinfo);
    /* this makes the library read from infile */
    jpeg_stdio_src(&cinfo, infile);
    /* reading the image header which contains image information */
    jpeg_read_header(&cinfo, TRUE);
    /* Start decompression jpeg here */
    jpeg_start_decompress(&cinfo);

    /* allocate memory to hold the uncompressed image */
    raw_image = (unsigned char *) malloc(cinfo.output_width * cinfo.output_height * cinfo.num_components);
    /* now actually read the jpeg into the raw buffer */
    row_pointer[0] = (unsigned char *) malloc(cinfo.output_width * cinfo.num_components);
    /* read one scan line at a time */
    while (cinfo.output_scanline < cinfo.image_height) {
        jpeg_read_scanlines(&cinfo, row_pointer, 1);
        for (i = 0; i < cinfo.image_width * cinfo.num_components; i++)
            raw_image[location++] = row_pointer[0][i];
    }

    height = cinfo.image_height;
    width = cinfo.image_width;

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, raw_image);
  



    glGenerateMipmap(GL_TEXTURE_2D);

    glUseProgram(shader); // don't forget to activate/use the shader before setting uniforms!
    // either set it manually like so:

    glUniform1i(glGetUniformLocation(shader, "TexGrey"), 1);
    /* wrap up decompression, destroy objects, free pointers and close open files */
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    free(row_pointer[0]);
    free(raw_image);
    fclose(infile);

}

void EclipseMap::initMoonColoredTexture(const char *filename, GLuint shader) {
    int width, height;
    glGenTextures(1, &moonTextureColor);
    cout << shader << endl;
    glBindTexture(GL_TEXTURE_2D, moonTextureColor);
    // set the texture wrapping parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                    GL_CLAMP_TO_EDGE);    // set texture wrapping to GL_REPEAT (default wrapping method)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    // set texture filtering parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    unsigned char *raw_image = NULL;
    int bytes_per_pixel = 3;   /* or 1 for GRACYSCALE images */
    int color_space = JCS_RGB; /* or JCS_GRAYSCALE for grayscale images */

    /* these are standard libjpeg structures for reading(decompression) */
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    /* libjpeg data structure for storing one row, that is, scanline of an image */
    JSAMPROW row_pointer[1];

    FILE *infile = fopen(filename, "rb");
    unsigned long location = 0;
    int i = 0, j = 0;

    if (!infile) {
        printf("Error opening jpeg file %s\n!", filename);
        return;
    }
    printf("Texture filename = %s\n", filename);

    /* here we set up the standard libjpeg error handler */
    cinfo.err = jpeg_std_error(&jerr);
    /* setup decompression process and source, then read JPEG header */
    jpeg_create_decompress(&cinfo);
    /* this makes the library read from infile */
    jpeg_stdio_src(&cinfo, infile);
    /* reading the image header which contains image information */
    jpeg_read_header(&cinfo, TRUE);
    /* Start decompression jpeg here */
    jpeg_start_decompress(&cinfo);

    /* allocate memory to hold the uncompressed image */
    raw_image = (unsigned char *) malloc(cinfo.output_width * cinfo.output_height * cinfo.num_components);
    /* now actually read the jpeg into the raw buffer */
    row_pointer[0] = (unsigned char *) malloc(cinfo.output_width * cinfo.num_components);
    /* read one scan line at a time */
    while (cinfo.output_scanline < cinfo.image_height) {
        jpeg_read_scanlines(&cinfo, row_pointer, 1);
        for (i = 0; i < cinfo.image_width * cinfo.num_components; i++)
            raw_image[location++] = row_pointer[0][i];
    }

    height = cinfo.image_height;
    width = cinfo.image_width;


    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, raw_image);
   

    imageWidth = width;
    imageHeight = height;

    glGenerateMipmap(GL_TEXTURE_2D);

    glUseProgram(shader); // don't forget to activate/use the shader before setting uniforms!
    // either set it manually like so:

    glUniform1i(glGetUniformLocation(shader, "MoonTexColor"), 2);
    /* wrap up decompression, destroy objects, free pointers and close open files */
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    free(row_pointer[0]);
    free(raw_image);
    fclose(infile);

}
