
#version 330 core
out vec4 FragColor;
  
in vec2 TexCoords;

uniform sampler2D screenTexture;

void main()
{ 
    float col = float(texture(screenTexture, TexCoords))/80000;
    FragColor = vec4(0, 0, col, 1);
}

