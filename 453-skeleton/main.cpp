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

	bool isDifferent(Parameters p) { // parameter for sin wave, for efficiency reasons, will check to see if any change was made (with the arrows from the user controlling them)
		return p.t != t || p.u != u;
	}

	int iterations = 0; // New parameter for Sierpinski iterations

};

// EXAMPLE CALLBACKS
class MyCallbacks : public CallbackInterface { // can be used to deal with user input

public:
	MyCallbacks(ShaderProgram& shader) : shader(shader), fractalType(0){}

	virtual void keyCallback(int key, int scancode, int action, int mods) {
		if (action == GLFW_PRESS || action == GLFW_REPEAT) { // if we press the right arrow on keyboard, will reload shader

			if (key == GLFW_KEY_0) {  // Load Sin Shape from tutorial
				fractalType = 0;
			}

			if (key == GLFW_KEY_1) {  // Load Sierpinski Triangle
				fractalType = 1;
			}

			if (key == GLFW_KEY_2) {  // Load Pythagoras Tree
				fractalType = 2;
			}
			

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

	int getFractalType() const {
		return fractalType;
	}

private:
	ShaderProgram& shader;
	Parameters parameters; // parameters is a private variable
	int fractalType; // to be able to switch fractals
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


CPU_Geometry generateSin(Parameters p) { // originally from main, but now put into own function for cpuGeom stuff. Done for efficiency, will only reload if some difference is detected!

	CPU_Geometry cpuGeom;
	for (float x = -1.0f; x <= 1.0f; x += 0.01f) {
		cpuGeom.verts.push_back(glm::vec3(x, sin(x*p.t)*p.u, 0.0));
		cpuGeom.cols.push_back(glm::vec3(sin(x), cos(x), 0.0));
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

		// Generate a random color for each iteration, this bascially separates each triangle by a separate colour
		glm::vec3 color(static_cast<float>(rand()) / RAND_MAX,
			static_cast<float>(rand()) / RAND_MAX,
			static_cast<float>(rand()) / RAND_MAX);

		// Add the same color for all three vertices of the triangle
		cpuGeom.cols.push_back(color);
		cpuGeom.cols.push_back(color);
		cpuGeom.cols.push_back(color);

		};

		

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


CPU_Geometry generatePythagorasTree(int iterations) {
	CPU_Geometry cpuGeom;

	// Function to draw a square
	auto draw_square = [&](const glm::vec3& bl, const glm::vec3& br, const glm::vec3& tr, const glm::vec3& tl) {
		// Add vertices in counterclockwise order
		cpuGeom.verts.push_back(bl);
		cpuGeom.verts.push_back(br);
		cpuGeom.verts.push_back(tr);
		cpuGeom.verts.push_back(tr);
		cpuGeom.verts.push_back(tl);
		cpuGeom.verts.push_back(bl);

		// Generate a random color for the square
		glm::vec3 color(static_cast<float>(rand()) / RAND_MAX,
			static_cast<float>(rand()) / RAND_MAX,
			static_cast<float>(rand()) / RAND_MAX);
		for (int i = 0; i < 6; ++i) {
			cpuGeom.cols.push_back(color);
		}
		};

	// Recursive function to build the Pythagoras Tree
	std::function<void(glm::vec3, glm::vec3, glm::vec3, glm::vec3, int)> add_tree = [&](glm::vec3 bl, glm::vec3 br, glm::vec3 tr, glm::vec3 tl, int depth) {
		if (depth == 0) {
			return;
		}

		// Draw the current square
		draw_square(bl, br, tr, tl);

		// Calculate the size of the new squares (should be sqrt(0.5) of the current square's side)
		glm::vec3 top_center = (tr + tl) * 0.5f;
		float new_size = glm::length(tr - tl) * sqrt(0.5f);

		// Vectors for new squares' orientation
		glm::vec3 left_dir = glm::normalize(tl - bl);
		glm::vec3 right_dir = glm::normalize(tr - br);

		// Calculate corners for the left branch (45-degree rotation)
		glm::vec3 new_tl_left = tl;
		glm::vec3 new_tr_left = new_tl_left + left_dir * new_size;
		glm::vec3 new_br_left = new_tl_left + left_dir * new_size + right_dir * new_size;
		glm::vec3 new_bl_left = new_tl_left + right_dir * new_size;

		// Calculate corners for the right branch (-45-degree rotation)
		glm::vec3 new_tr_right = tr;
		glm::vec3 new_tl_right = new_tr_right - right_dir * new_size;
		glm::vec3 new_bl_right = new_tr_right - right_dir * new_size - left_dir * new_size;
		glm::vec3 new_br_right = new_tr_right - left_dir * new_size;

		// Recursively add the branches
		add_tree(new_bl_left, new_br_left, new_tr_left, new_tl_left, depth - 1);
		add_tree(new_bl_right, new_br_right, new_tr_right, new_tl_right, depth - 1);
		};

	// Start the tree with the base square
	glm::vec3 p1(-0.1f, -0.5f, 0.f); // bottom left
	glm::vec3 p2(0.1f, -0.5f, 0.f);  // bottom right
	glm::vec3 p3(0.1f, -0.3f, 0.f);  // top right
	glm::vec3 p4(-0.1f, -0.3f, 0.f); // top left

	// Begin the recursive drawing
	add_tree(p1, p2, p3, p4, iterations);

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
	

	// GEOMETRY
	CPU_Geometry cpuGeom;
	GPU_Geometry gpuGeom;

	int lastFractalType = 0;

	
	for (float x = -1.0f; x <= 1.0f; x += 0.01f) {
		cpuGeom.verts.push_back(glm::vec3(x, sin(x*10)*0.5, 0.0));
		cpuGeom.cols.push_back(glm::vec3(cos(x), sin(x), 0.0));
	}


	Parameters p;
	cpuGeom = generateSin(p); // to initially load something, we'll start with the sin 
	//cpuGeom = generateSierpinski(p.iterations);


	// uploading data to gpu from cpu (note here we're uploading BEFORE THE LOOP, because in this example we're just redrawing the same traingle over and over again)
	gpuGeom.setVerts(cpuGeom.verts);
	gpuGeom.setCols(cpuGeom.cols);

	// RENDER LOOP
	// Keeps running until user decides to press esc or something
	// this loop runs EVERY frame
	while (!window.shouldClose()) {
		glfwPollEvents();

		Parameters newP = callbacks->getParameters();
		int currentFractalType = callbacks->getFractalType();

		if (currentFractalType != lastFractalType) {
			if (currentFractalType == 0) {
				for (float x = -1.0f; x <= 1.0f; x += 0.01f) {
					cpuGeom.verts.push_back(glm::vec3(x, sin(x * 10) * 0.5, 0.0));
					cpuGeom.cols.push_back(glm::vec3(cos(x), sin(x), 0.0));
				}
				cpuGeom = generateSin(p);
			}
			else if (currentFractalType == 1) {
				cpuGeom = generateSierpinski(p.iterations);
			}
			else if (currentFractalType == 2) {
				cpuGeom = generatePythagorasTree(p.iterations);
			}
		}

		if (currentFractalType == 0) {
			if (newP.isDifferent(p)) { // if difference detected, then we will reupload things
				p = newP;
				cpuGeom = generateSin(p);
				gpuGeom.setVerts(cpuGeom.verts);
				gpuGeom.setCols(cpuGeom.cols);
			}
		}

		if (currentFractalType == 1) {
			if (newP.iterations != p.iterations) {
				p = newP;
				cpuGeom = generateSierpinski(p.iterations);
				gpuGeom.setVerts(cpuGeom.verts);
				gpuGeom.setCols(cpuGeom.cols);
			}
		}

		if (currentFractalType == 2) {
			if (newP.iterations != p.iterations) {
				p = newP;
				cpuGeom = generatePythagorasTree(p.iterations);
				gpuGeom.setVerts(cpuGeom.verts);
				gpuGeom.setCols(cpuGeom.cols);
			}
		}

			shader.use(); // tells us to use the shader we loaded in earlier on
			gpuGeom.bind(); // tells us to use the gpu geometry we loaded in earlier on

			glEnable(GL_FRAMEBUFFER_SRGB); // tells us to use the srgb color space, flag used (OpenGL is a global state machine)
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clearing the current buffers (2, have to clear both). Two buffers, each split into 3 parts, R, G, B, R, G, B

			if (currentFractalType == 0) {
				glDrawArrays(GL_LINE_STRIP, 0, GLsizei(cpuGeom.verts.size())); // drawing on the buffer
			}

			if (currentFractalType == 1) {
				glDrawArrays(GL_TRIANGLES, 0, GLsizei(cpuGeom.verts.size())); // drawing on the buffer
			}

			if (currentFractalType == 2) {
				glDrawArrays(GL_TRIANGLES, 0, GLsizei(cpuGeom.verts.size())); // drawing on the buffer
			}

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
