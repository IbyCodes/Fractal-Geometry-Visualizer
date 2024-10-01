// LEARNOPENGL - GOOD SITE TO learn how to draw things properly
// clip space represents the space in which the GPU operates
// everything within a clip space gets drawn into the screen, everything else doesn't
// this clip space is a 2x2 cube centered at the origin, middle of cube
// towards the top is pos y axis, to the right is pos x axis, and out of the screen is pos z axis
// basically, every vertice you provide for this assignment must be between (-1,-1,-1) and (1,1,1)
// opengl shaders are in glsl (dont need to know about for assignment #1)
// FOLLOWED PARAMETERS CODE in tutorial to deal with the fact that we have multiple shapes that we gotta draw, not just 1

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"


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
#include <glm/glm.hpp>
#include <cmath>
#define M_PI 3.14159265358979323846 // PI
#include <vector>
#include <functional>
#include <glm/gtc/matrix_transform.hpp> // For transformations
#include <stack>


// Variables to manage the GUI state
int currentSceneGUI = 0; // Scene selector
int fractalDepthGUI = 0; // Depth/iteration of fractal shapes



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

			if (key == GLFW_KEY_3) {  // Load Koch Snowflake
				fractalType = 3;
			}

			if (key == GLFW_KEY_4) {
				fractalType = 4;
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

	void setParameters(const Parameters& newParameters) {
		parameters = newParameters;
	}


	int getFractalType() const {
		return fractalType;
	}

	void setFractalType(int const fractalInput) {
		fractalType = fractalInput;
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

// Function to blend between two colors based on a ratio
glm::vec3 blend(const glm::vec3& c1, const glm::vec3& c2, float ratio) {
	return c1 * (1 - ratio) + c2 * ratio;
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

	// Defining three base colors
	glm::vec3 colorTop = glm::vec3(0.0f, 0.5f, 1.0f);    // Blue for top section
	glm::vec3 colorLeft = glm::vec3(0.f, 0.2f, 0.f);   // Green for bottom-left
	glm::vec3 colorRight = glm::vec3(1.0f, 0.5f, 0.5f);  // Pink for bottom-right


	// Function to draw a single triangle by pushing its vertices
	auto draw_triangle = [&](const glm::vec3& a, const glm::vec3& b, const glm::vec3& c) {
		cpuGeom.verts.push_back(a);
		cpuGeom.verts.push_back(b);
		cpuGeom.verts.push_back(c);

		// Calculate the centroid to determine which region the triangle belongs to
		glm::vec3 centroid = (a + b + c) / 3.0f;

		// Assign color based on the region of the triangle

		glm::vec3 color;
		if (centroid.y > 0.0f) {
			color = colorTop;  // Top region
		}
		else if (centroid.x < 0.0f) {
			color = colorLeft; // Bottom-left region
		}
		else {
			color = colorRight; // Bottom-right region
		}
		

		// Add the same color for all three vertices of the triangle
		cpuGeom.cols.push_back(color);
		cpuGeom.cols.push_back(color);
		cpuGeom.cols.push_back(color);

		};

		

	// Recursive function to subdivide the triangle (from the notes)
	std::function<void(glm::vec3, glm::vec3, glm::vec3, int)> divide_triangle = [&](glm::vec3 a, glm::vec3 b, glm::vec3 c, int m) {
		if (m > 0) { // where m is the iteration level
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



// OMG THIS WORKS
CPU_Geometry generatePythagorasTree(int iterations) {
	CPU_Geometry cpuGeom;

	// Function to draw a square with the given color based on depth
	auto draw_square = [&](const glm::vec3& bl, const glm::vec3& br, const glm::vec3& tr, const glm::vec3& tl, int current_iteration) {
		// Adding vertices in counterclockwise order
		cpuGeom.verts.push_back(bl); // bottom left
		cpuGeom.verts.push_back(br); // bottom right
		cpuGeom.verts.push_back(tr); // top right
		cpuGeom.verts.push_back(tr); // top right
		cpuGeom.verts.push_back(tl); // top left
		cpuGeom.verts.push_back(bl); // bottom left

		// Set colors based on the iteration level
		glm::vec3 color;
		if (current_iteration == 0) {
			color = glm::vec3(0.30f, 0.16f, 0.14f);  // Brown for the first stem (initial trunk)
		}
		else if (current_iteration == 1) { // dark green
			color = glm::vec3(0.f, 0.1f, 0.f); 
		}
		else if (current_iteration == 2){ // lighter then last green
			color = glm::vec3(0.f, 0.2f, 0.f);
		}
		else if (current_iteration == 3) { // lighter then last green
			color = glm::vec3(0.f, 0.4f, 0.f);
		}
		else if (current_iteration == 4) { // lighter then last green
			color = glm::vec3(0.f, 0.7f, 0.f);
		}
		else { // lightest green in tree
			color = glm::vec3(0.f, 1.f, 0.f);
		}

		// Assign the chosen color to the square
		for (int i = 0; i < 6; ++i) {
			cpuGeom.cols.push_back(color);
		}
		};

	// Recursive function to build the Pythagoras Tree
	std::function<void(glm::vec3, glm::vec3, glm::vec3, glm::vec3, int, float)> add_tree =
		[&](glm::vec3 bl, glm::vec3 br, glm::vec3 tr, glm::vec3 tl, int depth, float oldSideLength) {
		int current_iteration = iterations - depth;  // Calculate the "current iteration" level

		if (depth < 0) {
			return;
		}

		draw_square(bl, br, tr, tl, current_iteration);  // Draw the current square with appropriate color

		if (depth == 0) return;  // Stop when depth is zero

		// Convert angles to radians
		float angle_radLeft = -45.0f * (M_PI / 180.0f);
		float angle_radRight = 45.0f * (M_PI / 180.0f);

		// Calculate the side length for the next squares
		float newSideLength = (sqrt(2) / 2) * oldSideLength;

		// Generate direction vectors for the current square
		glm::vec3 v = br - bl;                // Base vector (bottom side of square)
		glm::vec3 v_perpendicular = tl - bl;  // Perpendicular vector (left side of square)

		// Rotation coefficients for left square
		glm::vec3 left_v = cos(-angle_radLeft) * v + sin(-angle_radLeft) * v_perpendicular;
		glm::vec3 left_v_perpendicular = -sin(-angle_radLeft) * v + cos(-angle_radLeft) * v_perpendicular;

		// Calculating new square vertices for the left square
		glm::vec3 new_bl_left = tl;
		glm::vec3 new_br_left = new_bl_left + left_v * (newSideLength / glm::length(left_v));
		glm::vec3 new_tr_left = new_br_left + left_v_perpendicular * (newSideLength / glm::length(left_v_perpendicular));
		glm::vec3 new_tl_left = new_bl_left + left_v_perpendicular * (newSideLength / glm::length(left_v_perpendicular));

		// Rotation coefficients for right square
		glm::vec3 right_v = cos(-angle_radRight) * v + sin(-angle_radRight) * v_perpendicular;
		glm::vec3 right_v_perpendicular = -sin(-angle_radRight) * v + cos(-angle_radRight) * v_perpendicular;

		// Calculate new square vertices for the right square
		glm::vec3 new_bl_right = tr - right_v * (newSideLength / glm::length(right_v));
		glm::vec3 new_br_right = tr;
		glm::vec3 new_tr_right = new_br_right + right_v_perpendicular * (newSideLength / glm::length(right_v_perpendicular));
		glm::vec3 new_tl_right = new_bl_right + right_v_perpendicular * (newSideLength / glm::length(right_v_perpendicular));

		// Recursively adding branches to the left and right squares
		add_tree(new_bl_left, new_br_left, new_tr_left, new_tl_left, depth - 1, newSideLength);
		add_tree(new_bl_right, new_br_right, new_tr_right, new_tl_right, depth - 1, newSideLength);
		};

	// Starting the tree with the base square
	glm::vec3 p1(-0.1f, -0.5f, 0.f); // bottom left
	glm::vec3 p2(0.1f, -0.5f, 0.f);  // bottom right
	glm::vec3 p3(0.1f, -0.3f, 0.f);  // top right
	glm::vec3 p4(-0.1f, -0.3f, 0.f); // top left

	// Begin the recursive drawing
	add_tree(p1, p2, p3, p4, iterations, p3.x - p4.x);

	return cpuGeom;
}



CPU_Geometry generateKochSnowflake(int depth) {
	CPU_Geometry cpuGeom;

	// Initial triangle vertices (equilateral triangle)
	glm::vec3 p1(-0.5f, -0.5f, 0.f);  // bottom-left
	glm::vec3 p2(0.5f, -0.5f, 0.f);   // bottom-right
	glm::vec3 p3(0.f, sqrt(3.0f) / 2.0f - 0.5f, 0.f);  // top vertex of the equilateral triangle


	// Recursive function to divide and create the Koch snowflake with changing colors (sort of similar to Serpinsi triangle)
	std::function<void(glm::vec3, glm::vec3, int)> divideKoch = [&](glm::vec3 a, glm::vec3 b, int m) {
		if (m > 0) {
			// Calculate 1/3 and 2/3 points along the line
			glm::vec3 P1 = a + (b - a) / 3.0f;           // 1/3 of side point
			glm::vec3 P2 = a + (b - a) * 2.0f / 3.0f;    // 2/3 of side point

			// Direction of the line from P1 to P2
			glm::vec3 dir = P2 - P1;

			// Perpendicular vector to the segment, rotated 60 degrees outward
			glm::vec3 Ppeak(
				P1.x + dir.x * 0.5f + dir.y * sqrt(3.0f) / 2.0f,   // Rotate counter-clockwise
				P1.y - dir.x * sqrt(3.0f) / 2.0f + dir.y * 0.5f,
				0.0f
			);

			// Recursively subdividing each segment
			divideKoch(a, P1, m - 1);   // Segment 1
			divideKoch(P1, Ppeak, m - 1);  // Segment 2 with peak
			divideKoch(Ppeak, P2, m - 1);  // Segment 3
			divideKoch(P2, b, m - 1);   // Segment 4
		}
		else {
			// If recursion depth is reached, add final vertices to geometry
			cpuGeom.verts.push_back(a);
			cpuGeom.verts.push_back(b);
			glm::vec3 yellowColor = glm::vec3(1.0f, 1.0f, 0.0f);
			cpuGeom.cols.push_back(yellowColor); // first half of star
			cpuGeom.cols.push_back(yellowColor); // 2nd half of star
		}
		};

	// Start the recursive subdivision for each side of the initial triangle
	divideKoch(p1, p2, depth);  // First side
	divideKoch(p2, p3, depth);  // Second side
	divideKoch(p3, p1, depth);  // Third side

	return cpuGeom;
}


// sources used to help for this: https://rosettacode.org/wiki/Dragon_curve, https://www.csharphelper.com/howtos/howto_heighway_dragon.html, https://tfetimes.com/c-dragon-curve/
// 
// Defining directions for easy understanding
enum Directions { NORTH = 0, EAST, SOUTH, WEST };

// Generates a Dragon Curve fractal
CPU_Geometry generateDragonCurve(int iterations) {
	CPU_Geometry cpuGeometry;
	std::vector<glm::vec2> pointsList; // will store dragon curve points

	// Starting points for the curve
	pointsList.push_back(glm::vec2(-0.7, 0.2)); // Initial point
	pointsList.push_back(glm::vec2(0.5, 0.2));  // Final point

	// Defining a color palette for the segments (I need to figure out now to color it properly) 
	std::vector<glm::vec3> colors = {
		glm::vec3(1.0, 0.0, 0.0), // Red
		glm::vec3(0.0, 1.0, 0.0), // Green
		glm::vec3(0.0, 0.0, 1.0), // Blue
		glm::vec3(1.0, 1.0, 0.0), // Yellow
		glm::vec3(0.0, 1.0, 1.0), // Cyan
		glm::vec3(1.0, 0.0, 1.0), // Magenta
		glm::vec3(0.5, 0.5, 0.5), // Gray
		glm::vec3(1.0, 0.5, 0.0), // Orange
		glm::vec3(0.5, 0.0, 0.5)  // Purple
	};

	Directions direction = NORTH; // Initialize with North direction
	bool turnRight = true;

	// Defining lambda functions for computing left and right points
	auto computeLeftPoint = [](glm::vec2 startPoint, glm::vec2 endPoint, Directions direction) -> glm::vec2 {
		glm::vec2 midPoint((startPoint.x + endPoint.x) / 2, (startPoint.y + endPoint.y) / 2);
		glm::vec2 offset = midPoint - startPoint;

		switch (direction) {
		case NORTH:
			return glm::vec2(midPoint.x, midPoint.y + offset.x);
		case EAST:
			return glm::vec2(startPoint.x, endPoint.y);
		case SOUTH:
			return glm::vec2(midPoint.x - offset.y, midPoint.y);
		case WEST:
			return glm::vec2(endPoint.x, startPoint.y);
		default:
			return midPoint;
		}
		};

	auto computeRightPoint = [](glm::vec2 startPoint, glm::vec2 endPoint, Directions direction) -> glm::vec2 {
		glm::vec2 midPoint((startPoint.x + endPoint.x) / 2, (startPoint.y + endPoint.y) / 2);
		glm::vec2 offset = midPoint - startPoint;

		switch (direction) {
		case NORTH:
			return glm::vec2(midPoint.x, midPoint.y - offset.x);
		case EAST:
			return glm::vec2(endPoint.x, startPoint.y);
		case SOUTH:
			return glm::vec2(midPoint.x + offset.y, midPoint.y);
		case WEST:
			return glm::vec2(startPoint.x, endPoint.y);
		default:
			return midPoint;
		}
		};

	// Loop through the number of iterations to build the curve
	while (iterations > 1) {
		for (int index = 0; index < pointsList.size() - 1; index += 2) {
			if (!turnRight) {
				// Insert the point on the left side
				pointsList.insert(pointsList.begin() + index + 1, computeLeftPoint(pointsList.at(index), pointsList.at(index + 1), direction));
				turnRight = true;
			}
			else {
				// Insert the point on the right side
				pointsList.insert(pointsList.begin() + index + 1, computeRightPoint(pointsList.at(index), pointsList.at(index + 1), direction));
				turnRight = false;
			}
			// Update direction after two points are inserted
			direction = static_cast<Directions>((direction + 2) % 4);
		}
		// Increment direction for next iteration
		direction = static_cast<Directions>((direction + 1) % 4);
		iterations--;
	}

	// Safeguard against divide-by-zero issue
	int pointsPerColor = (pointsList.size() / colors.size()) > 0 ? (pointsList.size() / colors.size()) : 1;

	// Apply colors and handle X-axis mirroring for visual symmetry (I wasn't sure if I had to show the curve exactly like the assignment or not, so I decided to mirror my points just in case)
	for (int i = 0; i < pointsList.size() - 1; i++) {
		glm::vec2 startPoint = pointsList[i];
		glm::vec2 endPoint = pointsList[i + 1];

		// Invert X-coordinates for visualization
		startPoint.x = -startPoint.x;
		endPoint.x = -endPoint.x;

		// Add the starting point and color for the segment based on the adjusted iteration index
		cpuGeometry.verts.push_back(glm::vec3(startPoint.x, startPoint.y, 0.f));
		cpuGeometry.cols.push_back(colors[(i / pointsPerColor) % colors.size()]); // Assign color per segment, considering the corrected pointsPerColor

		// Add the ending point and the same color (to form a complete line segment)
		cpuGeometry.verts.push_back(glm::vec3(endPoint.x, endPoint.y, 0.f));
		cpuGeometry.cols.push_back(colors[(i / pointsPerColor) % colors.size()]); // Maintain the same color for the complete line segment
	}

	return cpuGeometry;
}




int main() {
	Log::debug("Starting main");

	// WINDOW
	glfwInit(); // we are using glfw, a openGL library

	Window window(800, 800, "CPSC 453 A1"); // can set callbacks at construction if desired

	
	window.setupImGui();  // SOURCE USED: https://www.youtube.com/watch?v=VRwhNKoxUtk&ab_channel=VictorGordan


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


		window.startImGuiFrame();


		// Scene Selector
		ImGui::Begin("Fractal Menu by Mohammad Khan");
		const char* scenes[] = { "Sin Wave (part of tutorial template)", "Sierpinski Triangle", "Pythagoras Tree", "Koch Snowflake", "Dragon Curve" };
		ImGui::Text("Select Scene:");
		bool sceneChangedGUI = ImGui::Combo("##scene_selector", &currentSceneGUI, scenes, 5);
		// Fractal Depth Slider
		ImGui::Text("Fractal Iteration:");
		bool depthChangedGUI = ImGui::SliderInt("##fractal_depth", &fractalDepthGUI, 0, 10); // Setting iteration range from 0 to 10
		ImGui::End();
		
	
		Parameters newP = callbacks->getParameters();
		int currentFractalType = callbacks->getFractalType();


		if (depthChangedGUI) {
			std::cout << "Depth changed by GUI to: " << fractalDepthGUI << std::endl;
			newP.iterations = fractalDepthGUI; // Cap at 10 iterations // FOR GUI, I THINK THIS IS WORKING AS EXPECTED.	

		}
		else {
			fractalDepthGUI = newP.iterations;
		}

		if (sceneChangedGUI) {
			std::cout << "Scene changed by GUI to: " << currentSceneGUI << std::endl;
			currentFractalType = currentSceneGUI;
			newP.iterations = 0; // to avoid crashing
		}
		else {
			currentSceneGUI = currentFractalType;
		}

		callbacks->setParameters(newP);  // FINALLY WORKING GUI INPUTS
		callbacks->setFractalType(currentSceneGUI); // FINALLY WORKING GUI INPUTS


		if (depthChangedGUI || sceneChangedGUI || currentFractalType != lastFractalType) {
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
			else if (currentFractalType == 3) {
		
				cpuGeom = generateKochSnowflake(p.iterations);
			}

			else if (currentFractalType == 4) {
			
				cpuGeom = generateDragonCurve(p.iterations);
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

		if (currentFractalType == 3) {
			if (newP.iterations != p.iterations) {
				p = newP;
				cpuGeom = generateKochSnowflake(p.iterations);
				gpuGeom.setVerts(cpuGeom.verts);
				gpuGeom.setCols(cpuGeom.cols);
				
			}
		}

		if (currentFractalType == 4) {
			if (newP.iterations != p.iterations) {
				p = newP;
				cpuGeom = generateDragonCurve(p.iterations);
				gpuGeom.setVerts(cpuGeom.verts);
				gpuGeom.setCols(cpuGeom.cols);
			}
		}

			shader.use(); // tells us to use the shader we loaded in earlier on
			gpuGeom.bind(); // tells us to use the gpu geometry we loaded in earlier on

			glEnable(GL_FRAMEBUFFER_SRGB); // tells us to use the srgb color space, flag used (OpenGL is a global state machine)
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clearing the current buffers (2, have to clear both). Two buffers, each split into 3 parts, R, G, B, R, G, B



			if (currentFractalType == 0) {
				glDrawArrays(GL_LINE_STRIP, 0, GLsizei(cpuGeom.verts.size())); // drawing on the buffer (sin wave)
			}

			if (currentFractalType == 1) {
				glDrawArrays(GL_TRIANGLES, 0, GLsizei(cpuGeom.verts.size())); // drawing on the buffer (triangle fractal)
			}

			if (currentFractalType == 2) {
				glDrawArrays(GL_TRIANGLES, 0, GLsizei(cpuGeom.verts.size())); // drawing on the buffer (tree)
			}

			if (currentFractalType == 3) {
				glDrawArrays(GL_LINE_LOOP , 0, GLsizei(cpuGeom.verts.size())); // drawing on the buffer (snowflake)
			}

			if (currentFractalType == 4) {
				glDrawArrays(GL_LINE_STRIP, 0, GLsizei(cpuGeom.verts.size())); // drawing on the buffer (dragon curve)
			}

			glDisable(GL_FRAMEBUFFER_SRGB); // disable sRGB for things like imgui

			window.renderImGui();


			window.swapBuffers(); // we're using a dual buffer system, so once we're done drawing on one buffer we swap it
		}


	// NOTE: openGL likes a counter-clockwise order format, so try to draw vertices via counter-clockwise order
	// theres several things you can draw. Some examples:
	// GL_POINTS (just draws single points) 1:1 vertices
	// GL_LINES (draws lines between two points) 2:1 vertices (NOTE EXAMPLE OF BEHAVIOUR: connects 1st and 2nd point, then connects 3rd and 4th point, but wouldn't connect 2nd and 3rd point together)
	// GL_LINE_STRIP (n-1 segments, if 5 points, 4 lines, just connecting all points together (no loop))
	// GL_LINE_LOOP (connects all dots together like a triangle, square, etc) n segments, similar to strip but forms a loop shape
	// GL_TRIANGLES (takes 3 vertices, produces 1 triangle. Fills in everything in betwen the triangle)


	window.shutdownImGui();

	glfwTerminate();
	return 0;
}
