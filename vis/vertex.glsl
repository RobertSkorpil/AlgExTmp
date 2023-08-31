#version 330 core
layout (location = 0) in vec4 aPos;
layout (location = 1) in vec4 aVel;

out vec4 my_color;

float resize(float f, float m)
{
	return (f / m);
}

void main()
{
	float alpha = 0;//3.14159 / 2;
	float beta = 0.;
	mat4x4 iso1, iso2;
	iso1[0] = vec4(1, 0, 0, 0);
	iso1[1] = vec4(0, cos(alpha), sin(alpha), 0);
	iso1[2] = vec4(0, -sin(alpha), cos(alpha), 0);
	iso1[3] = vec4(0, 0, 0, 1);

	iso2[0] = vec4(cos(beta), 0, -sin(beta), 0);
	iso2[1] = vec4(0, 1, 0, 0);
	iso2[2] = vec4(sin(beta), 0, cos(beta), 0);
	iso2[3] = vec4(0, 0, 0, 1);

	mat4x4 flip;
	flip[0] = vec4(0, 1, 0, 0);
	flip[1] = vec4(0, 0, 1, 0);
	flip[2] = vec4(1, 0, 0, 0);
	flip[3] = vec4(0, 0, 0, 1);

	mat4x4 iso = iso2 * iso1;
	vec4 scaled = vec4(aPos[1], aPos[2], aPos[3], 30);
	gl_Position = scaled * iso;
	my_color = aVel;
}