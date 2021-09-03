#version 410 core

// Output for on-screen color
layout(location = 0) out float outColor;

// Interpolated output data from vertex shader
in float fragId;
in vec3 fragBary;

void main() 
{
	outColor = fragId;
}

	//if (secondpass==1) col = float(texture2D(rsttex, gl_FragCoord.xy))/2.f;
	//{
	//	float idValue = float(texture2D(rsttex, gl_FragCoord.xy));
	//	if (fragId == idValue) discard;
	//}
	//if (fragBary.x < 1E-4 || fragBary.y < 1E-4 ||fragBary.z < 1E-4) col = -1.f;
	 //depth = gl_FragCoord.z
	 //uniform int secondpass;
//uniform sampler2D rsttex;