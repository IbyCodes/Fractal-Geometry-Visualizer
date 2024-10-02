// CPSC 453 Assignment 1
// Mohammad Khan 30103764

// NOTES FOR MYSELF:
// LEARNOPENGL - GOOD SITE TO learn how to draw things properly
// clip space represents the space in which the GPU operates
// everything within a clip space gets drawn into the screen, everything else doesn't
// this clip space is a 2x2 cube centered at the origin, middle of cube
// towards the top is pos y axis, to the right is pos x axis, and out of the screen is pos z axis
// basically, every vertice you provide for this assignment must be between (-1,-1,-1) and (1,1,1)
// opengl shaders are in glsl (dont need to know about for assignment #1)
// FOLLOWED PARAMETERS CODE in tutorial video to deal with the fact that we have multiple shapes that we gotta draw, not just 1

// includes for the GUI portion of the code
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

// PLEASE NOTE: some code is to do with the sin wave from the tutorial video on D2L, I decided to build on top of that code and kept that shape in my assignment! You may ignore it if you wish.


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

	int iterations = 0; // New parameter for Sierpinski shape and all other iterations for other shapes

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

			if (key == GLFW_KEY_4) { // Load Dragon Curve
				fractalType = 4;
			}
			

			if (key == GLFW_KEY_R) { // if we press R, the image will reload
				shader.recompile();
			}
			if (key == GLFW_KEY_UP) { // if we press the up arrow the 't' parameter will increase a little 
				parameters.t += parameters.tStep; // (for sin wave)
				parameters.iterations = std::min(parameters.iterations + 1, 10); // Cap at 10 iterations
			}

			if (key == GLFW_KEY_DOWN) { // if we press the down arrow the 't' parameter will decrease a little
				parameters.t -= parameters.tStep;  // (for sin wave)
				parameters.iterations = std::max(parameters.iterations - 1, 0); // Minimum 0 iterations
			}

			if (key == GLFW_KEY_LEFT) { // if we press the LEFT arrow the 'u' parameter will decrease a little (only works for sin wave)
				parameters.u -= parameters.uStep;
			}

			if (key == GLFW_KEY_RIGHT) { // if we press the RIGHT arrow the 'u' parameter will increase a little (only works for sin wave)
				parameters.u += parameters.uStep;
			}

		}
	}

	Parameters getParameters() { // getter function to get parameters
		return parameters;
	}

	void setParameters(const Parameters& newParameters) { // setter function to set parameters
		parameters = newParameters;
	}


	int getFractalType() const { // getter function to get the current fractal in window
		return fractalType;
	}

	void setFractalType(int const fractalInput) { // setter function to set the current fractal in window
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


CPU_Geometry generateSin(Parameters p) { // originally from main, but now put into own function for cpuGeom stuff. Done for efficiency, will only reload if some difference is detected! (you may choose to ignore this code, its for the sin wave from the tutorial)

	CPU_Geometry cpuGeom;
	for (float x = -1.0f; x <= 1.0f; x += 0.01f) {
		cpuGeom.verts.push_back(glm::vec3(x, sin(x*p.t)*p.u, 0.0));
		cpuGeom.cols.push_back(glm::vec3(sin(x), cos(x), 0.0));
	}
	return cpuGeom;

}



// Function to blend between two colors based on a ratio (will be used for sierpinksi)
glm::vec3 blend(const glm::vec3& c1, const glm::vec3& c2, float ratio) {
	return c1 * (1 - ratio) + c2 * ratio;
}

CPU_Geometry generateSierpinski(int iterations) { // function to generate the sierpinksi triangle
	CPU_Geometry cpuGeom;

	
	/*
	// function to subdivide the sierpinksi triangle fractal (followed this pseudocode for shape, from CPSC 453 notes on D2L)
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

		// Calculating the centroid to determine which region the triangle belongs to
		glm::vec3 centroid = (a + b + c) / 3.0f;

		// Assigning color based on the region of the triangle
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
		

		// Adding the same color for all three vertices of the triangle
		cpuGeom.cols.push_back(color);
		cpuGeom.cols.push_back(color);
		cpuGeom.cols.push_back(color);

		};

		

	// Recursive function to subdivide the triangle (from the notes)
	std::function<void(glm::vec3, glm::vec3, glm::vec3, int)> divide_triangle = [&](glm::vec3 a, glm::vec3 b, glm::vec3 c, int m) {
		if (m > 0) { // where m is the iteration level
			// Calculating midpoints of each side
			glm::vec3 v0 = (a + b) * 0.5f;
			glm::vec3 v1 = (a + c) * 0.5f;
			glm::vec3 v2 = (b + c) * 0.5f;

			// Recursively dividing the triangles
			divide_triangle(a, v0, v1, m - 1);
			divide_triangle(c, v1, v2, m - 1);
			divide_triangle(b, v2, v0, m - 1);
		}
		else {
			// Drawing the triangle when the recursion depth is zero
			draw_triangle(a, b, c);
		}

		};

	// Initial triangle vertices (an equilateral triangle)
	glm::vec3 p1(-0.5f, -0.5f, 0.f);
	glm::vec3 p2(0.5f, -0.5f, 0.f);
	glm::vec3 p3(0.f, 0.5f, 0.f);

	// Starting the recursive subdivision
	divide_triangle(p1, p2, p3, iterations);

	return cpuGeom;
}




CPU_Geometry generatePythagorasTree(int iterations) { // Main function to generate the Pythagoras Tree
	CPU_Geometry cpuGeom;

	// Function to draw a square with the given color based on depth (sort of similar to the serpinksi triangle coloring)
	auto draw_square = [&](const glm::vec3& bl, const glm::vec3& br, const glm::vec3& tr, const glm::vec3& tl, int current_iteration) {
		// Adding vertices in counterclockwise order
		cpuGeom.verts.push_back(bl); // bottom left
		cpuGeom.verts.push_back(br); // bottom right
		cpuGeom.verts.push_back(tr); // top right
		cpuGeom.verts.push_back(tr); // top right
		cpuGeom.verts.push_back(tl); // top left
		cpuGeom.verts.push_back(bl); // bottom left

		// Setting colors based on the iteration level
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
		else if (current_iteration == 3) { // even lighter then last green
			color = glm::vec3(0.f, 0.4f, 0.f);
		}
		else if (current_iteration == 4) { // even lighter then last green
			color = glm::vec3(0.f, 0.7f, 0.f);
		}
		else { // lightest green in tree, will keep this color for further iterations as well
			color = glm::vec3(0.f, 1.f, 0.f);
		}

		// Assigning the chosen color to the square
		for (int i = 0; i < 6; ++i) {
			cpuGeom.cols.push_back(color);
		}
		};

	// Recursive function to build the Pythagoras Tree
	std::function<void(glm::vec3, glm::vec3, glm::vec3, glm::vec3, int, float)> add_tree =
		[&](glm::vec3 bl, glm::vec3 br, glm::vec3 tr, glm::vec3 tl, int depth, float oldSideLength) {

		int current_iteration = iterations - depth;  // Calculating the "current iteration" level

		if (depth < 0) {
			return;
		}

		draw_square(bl, br, tr, tl, current_iteration);  // Drawing the current square with appropriate color

		if (depth == 0) return;  // Stopping when depth is zero

		// Converting angles to radians (I had no big reason to use rad instead of angles, I just searched it up and it seemed like C++ preferred rad over angles so I stuck with them)
		float angle_radLeft = -45.0f * (M_PI / 180.0f);
		float angle_radRight = 45.0f * (M_PI / 180.0f);

		// Calculating the side length for the next squares
		float newSideLength = (sqrt(2) / 2) * oldSideLength;

		// Generating direction vectors for the current square (I'm using the perpendicular tracing method. Source for guidance for this portion of my code: Riley James TA)
		glm::vec3 v = br - bl;                // Base vector (bottom side of square)
		glm::vec3 v_perpendicular = tl - bl;  // Perpendicular vector (left side of square)

		// Rotation coefficients for left square (had to do static cast to float for the code to work in QTCreator)
		glm::vec3 left_v = static_cast<float>(cos(-angle_radLeft)) * v + static_cast<float>(sin(-angle_radLeft)) * v_perpendicular;
		glm::vec3 left_v_perpendicular = static_cast<float>(-sin(-angle_radLeft)) * v + static_cast<float>(cos(-angle_radLeft)) * v_perpendicular;

		// Calculating new square vertices for the left square
		glm::vec3 new_bl_left = tl;
		glm::vec3 new_br_left = new_bl_left + left_v * (newSideLength / glm::length(left_v));
		glm::vec3 new_tr_left = new_br_left + left_v_perpendicular * (newSideLength / glm::length(left_v_perpendicular));
		glm::vec3 new_tl_left = new_bl_left + left_v_perpendicular * (newSideLength / glm::length(left_v_perpendicular));

		// Rotation coefficients for right square (had to do static cast to float for the code to work in QTCreator)
		glm::vec3 right_v = static_cast<float>(cos(-angle_radRight)) * v + static_cast<float>(sin(-angle_radRight)) * v_perpendicular;
		glm::vec3 right_v_perpendicular = static_cast<float>(- sin(-angle_radRight)) * v + static_cast<float>(cos(-angle_radRight)) * v_perpendicular;

		// Calculating new square vertices for the right square
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

	// Starting the recursive drawing
	add_tree(p1, p2, p3, p4, iterations, p3.x - p4.x);

	return cpuGeom;
}



CPU_Geometry generateKochSnowflake(int depth) { // main function to generate the Koch Snowflake
	CPU_Geometry cpuGeom;

	// Recursive function to divide and create the Koch snowflake + color it
	std::function<void(glm::vec3, glm::vec3, int)> divideKoch = [&](glm::vec3 a, glm::vec3 b, int m) {

		if (m > 0) {
			// Calculating 1/3 and 2/3 points along the line
			glm::vec3 P1 = a + (b - a) / 3.0f;           // 1/3 of side point
			glm::vec3 P2 = a + (b - a) * 2.0f / 3.0f;    // 2/3 of side point

			// Direction of the line from P1 to P2
			glm::vec3 dir = P2 - P1;

			// Perpendicular vector to the segment, rotated 60 degrees outward
			glm::vec3 Ppeak(
				P1.x + dir.x * 0.5f + dir.y * sqrt(3.0f) / 2.0f,   // Rotating counter-clockwise
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
			// If recursion depth is reached, adding final vertices to shape
			cpuGeom.verts.push_back(a);
			cpuGeom.verts.push_back(b);
			glm::vec3 yellowColor = glm::vec3(1.0f, 1.0f, 0.0f); // I stuck with just 1 color because the assignment didn't specify you needed any alternating colors
			cpuGeom.cols.push_back(yellowColor); // first half of the star being colored
			cpuGeom.cols.push_back(yellowColor); // 2nd half of the star being colored
		}
		};

	// Initial triangle vertices (equilateral triangle)
	glm::vec3 p1(-0.5f, -0.5f, 0.f);  // bottom-left
	glm::vec3 p2(0.5f, -0.5f, 0.f);   // bottom-right
	glm::vec3 p3(0.f, sqrt(3.0f) / 2.0f - 0.5f, 0.f);  // top vertex of the EQUILATERAL triangle (this was important as stated in the assignment description)

	// Starting the recursive subdivision for each side of the initial triangle
	divideKoch(p1, p2, depth);  // First side
	divideKoch(p2, p3, depth);  // Second side
	divideKoch(p3, p1, depth);  // Third side

	return cpuGeom;
}



// sources used to help for this dragon curve: https://rosettacode.org/wiki/Dragon_curve, https://www.csharphelper.com/howtos/howto_heighway_dragon.html, https://tfetimes.com/c-dragon-curve/

// Defining directions for easy understanding
enum Directions { NORTH = 0, EAST, SOUTH, WEST };


CPU_Geometry generateDragonCurve(int iterations) { // main function to draw the dragon curve
	CPU_Geometry cpuGeometry;
	std::vector<glm::vec2> pointsList; // will store dragon curve points

	// Starting points for the curve
	pointsList.push_back(glm::vec2(-0.7, 0.2)); // Initial point (starting)
	pointsList.push_back(glm::vec2(0.5, 0.2));  // Final point (starting)

	// Defining a color palette for the segments (I'm not too sure if this is coloring method that's expected or not, but stuck with this) 
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

	Directions direction = NORTH; // Initializing direction with North direction
	bool turnRight = true; // a boolean variable that will determine if we should flip right or not


	// Defining functions for computing left and right points of the curve

	auto computeLeftPoint = [](glm::vec2 startPoint, glm::vec2 endPoint, Directions direction) -> glm::vec2 { // will compute all left dragon curve points
		glm::vec2 midPoint((startPoint.x + endPoint.x) / 2, (startPoint.y + endPoint.y) / 2);
		glm::vec2 offset = midPoint - startPoint;

		switch (direction) { // easy switch to decide which vector should be drawn depending on direction
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
	 
	auto computeRightPoint = [](glm::vec2 startPoint, glm::vec2 endPoint, Directions direction) -> glm::vec2 { // will compute all right dragon curve points
		glm::vec2 midPoint((startPoint.x + endPoint.x) / 2, (startPoint.y + endPoint.y) / 2);
		glm::vec2 offset = midPoint - startPoint;

		switch (direction) { // easy switch to decide which vector should be drawn depending on direction
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

	// Looping through the number of iterations to build the curve
	while (iterations > 0) {
		for (int index = 0; index < pointsList.size() - 1; index += 2) {
			if (!turnRight) {
				// Inserting the point on the left side if we're not to turn right
				pointsList.insert(pointsList.begin() + index + 1, computeLeftPoint(pointsList.at(index), pointsList.at(index + 1), direction));
				turnRight = true;
			}
			else {
				// Inserting the point on the right side
				pointsList.insert(pointsList.begin() + index + 1, computeRightPoint(pointsList.at(index), pointsList.at(index + 1), direction));
				turnRight = false; // next turn should be left 
			}
			// Updating direction after two points are inserted
			direction = static_cast<Directions>((direction + 2) % 4);
		}
		// Incrementing direction for the next iteration
		direction = static_cast<Directions>((direction + 1) % 4);
		iterations--;
	}

	// Safeguard against divide-by-zero issue for the color (was running into this error sometimes)
	int pointsPerColor = (pointsList.size() / colors.size()) > 0 ? (pointsList.size() / colors.size()) : 1;

	// Applying colors and handle X-axis mirroring for visual symmetry (I wasn't sure if I had to show the curve exactly like the assignment or not, so I decided to mirror my points just in case)
	for (int i = 0; i < pointsList.size() - 1; i++) {
		glm::vec2 startPoint = pointsList[i];
		glm::vec2 endPoint = pointsList[i + 1];

		// Inverting X-coordinates for visualization (so it looks exactly like the assignment description)
		startPoint.x = -startPoint.x;
		endPoint.x = -endPoint.x;

		// Adding the starting point and color for the segment based on the adjusted iteration index
		cpuGeometry.verts.push_back(glm::vec3(startPoint.x, startPoint.y, 0.f));
		cpuGeometry.cols.push_back(colors[(i / pointsPerColor) % colors.size()]); // Assigning color per segment, considering the corrected pointsPerColor

		// Adding the ending point and the same color (to form a complete line segment)
		cpuGeometry.verts.push_back(glm::vec3(endPoint.x, endPoint.y, 0.f));
		cpuGeometry.cols.push_back(colors[(i / pointsPerColor) % colors.size()]); // Maintaining the same color for the complete line segment
	}

	return cpuGeometry;
}




int main() {
	Log::debug("Starting main");

	// WINDOW
	glfwInit(); // we are using glfw, a openGL library

	Window window(800, 800, "CPSC 453 A1"); // can set callbacks at construction if desired

	
	window.setupImGui();  // SOURCE USED for all GUI code: https://www.youtube.com/watch?v=VRwhNKoxUtk&ab_channel=VictorGordan


	GLDebug::enable();

	// SHADERS
	ShaderProgram shader("shaders/test.vert", "shaders/test.frag");

	// CALLBACKS
	auto callbacks = std::make_shared<MyCallbacks>(shader); // can also update callbacks to new ones
	window.setCallbacks(callbacks);


	// GEOMETRY
	CPU_Geometry cpuGeom;
	GPU_Geometry gpuGeom;

	int lastFractalType = 0; // to keep track of if fractal has changed

	
	for (float x = -1.0f; x <= 1.0f; x += 0.01f) { // part of initial sin wave code (from tutorial, can be ignored)
		cpuGeom.verts.push_back(glm::vec3(x, sin(x*10)*0.5, 0.0));
		cpuGeom.cols.push_back(glm::vec3(cos(x), sin(x), 0.0));
	}


	Parameters p;
	cpuGeom = generateSin(p); // to initially load something, we'll start with the sin 


	
	// uploading data to gpu from cpu (note here we're uploading BEFORE THE LOOP, because in this example we're just redrawing the same triangle over and over again)
	// (this is a note to myself, can be ignored)
	gpuGeom.setVerts(cpuGeom.verts);
	gpuGeom.setCols(cpuGeom.cols);
	


	// RENDER LOOP
	// Keeps running until user decides to press esc or something
	// this loop runs EVERY frame
	while (!window.shouldClose()) {
		glfwPollEvents();


		window.startImGuiFrame(); // from window.cpp 


		// Scene Selector 
		ImGui::Begin("Fractal Menu by Mohammad Khan"); // title of my GUI
		const char* scenes[] = { "Sin Wave (part of tutorial template)", "Sierpinski Triangle", "Pythagoras Tree", "Koch Snowflake", "Dragon Curve" }; // all scenes in menu
		ImGui::Text("Select Scene:");
		bool sceneChangedGUI = ImGui::Combo("##scene_selector", &currentSceneGUI, scenes, 5);
		// Fractal Depth Slider
		ImGui::Text("Fractal Iteration:");
		bool depthChangedGUI = ImGui::SliderInt("##fractal_depth", &fractalDepthGUI, 0, 10); // Setting iteration range from 0 to 10 for ALL shapes (meets requirements for assignment)
		ImGui::End();
		
	
		Parameters newP = callbacks->getParameters();
		int currentFractalType = callbacks->getFractalType();


		if (depthChangedGUI) { // if the depth is changed via the GUI
			std::cout << "Depth changed by GUI to: " << fractalDepthGUI << std::endl;
			newP.iterations = fractalDepthGUI; // Cap at 10 iterations 

		}
		else { // the GUI should be updated if the depth is changed via keyboard controls
			fractalDepthGUI = newP.iterations;
		}

		if (sceneChangedGUI) { // if the scene is changed via the GUI
			std::cout << "Scene changed by GUI to: " << currentSceneGUI << std::endl;
			currentFractalType = currentSceneGUI;
			newP.iterations = 0; // to avoid crashing/bogging down of system, will start new shape at 0th iteration
		}
		else { // the GUI should be updated if the scene is changed via keyboard controls
			currentSceneGUI = currentFractalType;
		}

		callbacks->setParameters(newP);  // to update to the new parameters
		callbacks->setFractalType(currentSceneGUI); // to update to the new scene


		if (depthChangedGUI || sceneChangedGUI || currentFractalType != lastFractalType) {

			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

			cpuGeom.verts.clear(); // to avoid any weird shapes in between switching shapes
			cpuGeom.cols.clear(); // to avoid any weird shapes in between switching shapes



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


	// NOTES for myself: openGL likes a counter-clockwise order format, so try to draw vertices via counter-clockwise order
	// theres several things you can draw. Some examples:
	// GL_POINTS (just draws single points) 1:1 vertices
	// GL_LINES (draws lines between two points) 2:1 vertices (NOTE EXAMPLE OF BEHAVIOUR: connects 1st and 2nd point, then connects 3rd and 4th point, but wouldn't connect 2nd and 3rd point together)
	// GL_LINE_STRIP (n-1 segments, if 5 points, 4 lines, just connecting all points together (no loop))
	// GL_LINE_LOOP (connects all dots together like a triangle, square, etc) n segments, similar to strip but forms a loop shape
	// GL_TRIANGLES (takes 3 vertices, produces 1 triangle. Fills in everything in betwen the triangle)


	window.shutdownImGui(); // shuts down GUI

	glfwTerminate();

	return 0;
}
