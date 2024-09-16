// LEARNOPENGL - GOOD SITE TO learn how to draw things properly
// clip space represents the space in which the GPU operates
// everything within a clip space gets drawn into the screen, everything else doesn't
// this clip space is a 2x2 cube centered at the origin, middle of cube
// towards the top is pos y axis, to the right is pos x axis, and out of the screen is pos z axis
// basically, every vertice you provide for this assignment must be between (-1,-1,-1) and (1,1,1)
// opengl shaders are in glsl (dont need to know about for assignment #1)
// FOLLOWED PARAMETERS CODE in tutorial to deal with the fact that we have multiple shapes that we gotta draw, not just 1

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include "Geometry.h"
#include "GLDebug.h"
#include "Log.h"
#include "ShaderProgram.h"
#include "Shader.h"
#include "Window.h"
#include <cstdlib>  // For rand() and RAND_MAX


struct Parameters { // struct for parameters for user input example
	float tStep = 0.1f;
	float uStep = 0.1f;

	float t = 1.0f;
	float u = 0.0f;

	bool isDifferent(Parameters p) { // for efficiency reasons, will check to see if any change was made (with the arrows from the user controlling them)
		return p.t != t || p.u != u;
	}

	int iterations = 0; // New parameter for Sierpinski iterations

};

// EXAMPLE CALLBACKS
class MyCallbacks : public CallbackInterface { // can be used to deal with user input

public:
	MyCallbacks(ShaderProgram& shader) : shader(shader) {}

	virtual void keyCallback(int key, int scancode, int action, int mods) {
		if (action == GLFW_PRESS || action == GLFW_REPEAT) { // if we press the right arrow on keyboard, will reload shader
			if (key == GLFW_KEY_R) { // if we press R, the image will reload
				shader.recompile();
			}
			if (key == GLFW_KEY_UP) { // if we press the up arrow the 't' parameter will increase a little
				parameters.t += parameters.tStep;
				parameters.iterations = std::min(parameters.iterations + 1, 10); // Cap at 10 iterations
			}

			if (key == GLFW_KEY_DOWN) { // if we press the down arrow the 't' parameter will decrease a little
				parameters.t -= parameters.tStep;
				parameters.iterations = std::max(parameters.iterations - 1, 0); // Minimum 0 iterations
			}

			if (key == GLFW_KEY_LEFT) { // if we press the LEFT arrow the 'u' parameter will decrease a little
				parameters.u -= parameters.uStep;
			}

			if (key == GLFW_KEY_RIGHT) { // if we press the RIGHT arrow the 'u' parameter will increase a little
				parameters.u += parameters.uStep;
			}

		}
	}

	Parameters getParameters() { // getter function to get parameters
		return parameters;
	}

private:
	ShaderProgram& shader;
	Parameters parameters; // parameters is a private variable
};

class MyCallbacks2 : public CallbackInterface {

public:
	MyCallbacks2() {}

	virtual void keyCallback(int key, int scancode, int action, int mods) {
		if (key == GLFW_KEY_R && action == GLFW_PRESS) {
			std::cout << "called back" << std::endl;
		}
	}
};
// END EXAMPLES (the above things are just example pieces of code)

/*
std::vector<glm::vec3> sqSpiral(int i){
	std::vector<glm::vec3> ret;
	float step = 0.2;
	bool axis = false;
	bool direction = false;
	glm::vec3 currentPos(0., 0., 0.);

	for (int it = 0; it < i; it++){
		glm::vec3 newpos;
		if(direction){

			if(axis){
				newpos = currentPos + glm::vec3(0., step, 0.);
			}
			else{
				newpos = currentPos + glm::vec3(step, 0., 0.);
			}
		}
		else{
			if (axis){
				newpos = currentPos + glm::vec3(0., -step, 0.)
			}
			else{
				newpos = currentPos + glm::vec3(-step, 0., 0.)
			}
		}
	}
	ret.push_back(newpos);
	currentPos = newpos;

} // then put colours in, etc, etc
*/


CPU_Geometry generateSin(Parameters p) { // originally from main, but now put into own function for cpuGeom stuff. Done for efficiency, will only reload if some difference is detected!

	CPU_Geometry cpuGeom;
	for (float x = -1.0f; x <= 1.0f; x += 0.01f) {
		cpuGeom.verts.push_back(glm::vec3(x, sin(x*p.t)*p.u, 0.0));
		cpuGeom.cols.push_back(glm::vec3(0.6, 0.2, 0.0));
	}
	return cpuGeom;

}

CPU_Geometry generateSierpinski(int iterations) {
	CPU_Geometry cpuGeom;

	// Initial triangle vertices (an equilateral triangle)
	glm::vec3 p1(-0.5f, -0.5f, 0.f);
	glm::vec3 p2(0.5f, -0.5f, 0.f);
	glm::vec3 p3(0.f, 0.5f, 0.f);


	/*
	// function to subdivide the sierpinksi triangle fractal (from notes)
	void divide_triangle(glm::vec3 pointA, glm::vec3 pointB, glm::vec3 pointC, int m) {
		glm::vec3 v0, v1, v2;
		int j;
		if (m > 0) {
			for (j = 0; j < 2; j++)
				v0[j] = (pointA[j] + pointB[j]) * 0.5;
			for (j = 0; j < 2; j++)
				v1[j] = (pointA[j] + pointC[j]) * 0.5;
			for (j = 0; j < 2; j++)
				v2[j] = (pointB[j] + pointC[j]) * 0.5;
			divide_triangle(pointA, v0, v1, m - 1);
			divide_triangle(pointC, v1, v2, m - 1);
			divide_triangle(pointB, v2, v0, m - 1);
		}
		else
			draw_triangle(pointA, pointB, pointC);
	};
	*/

	// Function to draw a single triangle by pushing its vertices
	auto draw_triangle = [&](const glm::vec3& a, const glm::vec3& b, const glm::vec3& c) {
		cpuGeom.verts.push_back(a);
		cpuGeom.verts.push_back(b);
		cpuGeom.verts.push_back(c);

		// Generate a random color
		glm::vec3 color(static_cast<float>(rand()) / RAND_MAX,
			static_cast<float>(rand()) / RAND_MAX,
			static_cast<float>(rand()) / RAND_MAX);

		// Add the same color for all three vertices of the triangle
		cpuGeom.cols.push_back(color);
		cpuGeom.cols.push_back(color);
		cpuGeom.cols.push_back(color);

		};
		/*	// Add colors for each vertex
		cpuGeom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f)); // Color for vertex a
		cpuGeom.cols.push_back(glm::vec3(0.0f, 1.0f, 0.0f)); // Color for vertex b
		cpuGeom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f)); // Color for vertex c
		};
		*/

	// Recursive function to subdivide the triangle
	std::function<void(glm::vec3, glm::vec3, glm::vec3, int)> divide_triangle = [&](glm::vec3 a, glm::vec3 b, glm::vec3 c, int m) {
		if (m > 0) {
			// Calculate midpoints of each side
			glm::vec3 v0 = (a + b) * 0.5f;
			glm::vec3 v1 = (a + c) * 0.5f;
			glm::vec3 v2 = (b + c) * 0.5f;

			// Recursively divide the triangles
			divide_triangle(a, v0, v1, m - 1);
			divide_triangle(c, v1, v2, m - 1);
			divide_triangle(b, v2, v0, m - 1);
		}
		else {
			// Draw the triangle when the recursion depth is zero
			draw_triangle(a, b, c);
		}

		};

	// Start the recursive subdivision
	divide_triangle(p1, p2, p3, iterations);

	return cpuGeom;
}


int main() {
	Log::debug("Starting main");

	// WINDOW
	glfwInit(); // we are using glfw, a openGL library
	Window window(800, 800, "CPSC 453 A1"); // can set callbacks at construction if desired

	GLDebug::enable();

	// SHADERS
	ShaderProgram shader("shaders/test.vert", "shaders/test.frag");

	// CALLBACKS
	auto callbacks = std::make_shared<MyCallbacks>(shader); // can also update callbacks to new ones
	window.setCallbacks(callbacks);
	//window.setCallbacks(std::make_shared<MyCallbacks>(shader)); // can also update callbacks to new ones

	// GEOMETRY
	CPU_Geometry cpuGeom;
	GPU_Geometry gpuGeom;

	/*
	for (float x = -1.0f; x <= 1.0f; x += 0.01f) {
		cpuGeom.verts.push_back(glm::vec3(x, sin(x*10)*0.5, 0.0));
		cpuGeom.cols.push_back(glm::vec3(cos(x), sin(x), 0.0));
	}
	*/


	/*
	// vertices
	cpuGeom.verts.push_back(glm::vec3(-0.5f, -0.5f, 0.f)); // bottom left
	cpuGeom.verts.push_back(glm::vec3(0.5f, -0.5f, 0.f)); // bottom right
	cpuGeom.verts.push_back(glm::vec3(0.f, 0.5f, 0.f)); // top of triangle

	// colours (these should be in linear space)
	cpuGeom.cols.push_back(glm::vec3(1.f, 0.f, 0.f));
	cpuGeom.cols.push_back(glm::vec3(0.f, 1.f, 0.f));
	cpuGeom.cols.push_back(glm::vec3(0.f, 0.f, 1.f));

	*/

	Parameters p;
	// cpuGeom = generateSin(p);
	cpuGeom = generateSierpinski(p.iterations);
	// uploading data to gpu from cpu (note here we're uploading BEFORE THE LOOP, because in this example we're just redrawing the same traingle over and over again)
	gpuGeom.setVerts(cpuGeom.verts);
	gpuGeom.setCols(cpuGeom.cols);

	// RENDER LOOP
	// Keeps running until user decides to press esc or something

	// this loop runs EVERY frame
	while (!window.shouldClose()) {
		glfwPollEvents();

		Parameters newP = callbacks->getParameters();

		/*
		if (newP.isDifferent(p)) { // if difference detected, then we will reupload things
			p = newP;
			cpuGeom = generateSin(p);
			gpuGeom.setVerts(cpuGeom.verts);
			gpuGeom.setCols(cpuGeom.cols);
		}
		*/

		if (newP.iterations != p.iterations) {
			p = newP;
			cpuGeom = generateSierpinski(p.iterations);
			gpuGeom.setVerts(cpuGeom.verts);
			gpuGeom.setCols(cpuGeom.cols);
		}

		shader.use(); // tells us to use the shader we loaded in earlier on
		gpuGeom.bind(); // tells us to use the gpu geometry we loaded in earlier on

		glEnable(GL_FRAMEBUFFER_SRGB); // tells us to use the srgb color space, flag used (OpenGL is a global state machine)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clearing the current buffers (2, have to clear both). Two buffers, each split into 3 parts, R, G, B, R, G, B
		glDrawArrays(GL_TRIANGLES, 0, GLsizei(cpuGeom.verts.size())); // drawing on the buffer
		glDisable(GL_FRAMEBUFFER_SRGB); // disable sRGB for things like imgui

		window.swapBuffers(); // we're using a dual buffer system, so once we're done drawing on one buffer we swap it
	}

	// NOTE: openGL likes a counter-clockwise order format, so try to draw vertices via counter-clockwise order
	// theres several things you can draw. Some examples:
	// GL_POINTS (just draws single points) 1:1 vertices
	// GL_LINES (draws lines between two points) 2:1 vertices (NOTE EXAMPLE OF BEHAVIOUR: connects 1st and 2nd point, then connects 3rd and 4th point, but wouldn't connect 2nd and 3rd point together)
	// GL_LINE_STRIP (n-1 segments, if 5 points, 4 lines, just connecting all points together (no loop))
	// GL_LINE_LOOP (connects all dots together like a triangle, square, etc) n segments, similar to strip but forms a loop shape
	// GL_TRIANGLES (takes 3 vertices, produces 1 triangle. Fills in everything in betwen the triangle)

	glfwTerminate();
	return 0;
}
