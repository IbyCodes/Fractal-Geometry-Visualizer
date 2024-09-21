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
#include <glm/glm.hpp>
#include <cmath>



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




CPU_Geometry generatePythagorasTree(int iterations) {
	CPU_Geometry cpuGeom;
	
	// Function to draw a square
	auto draw_square = [&](const glm::vec3& bl, const glm::vec3& br, const glm::vec3& tr, const glm::vec3& tl) {
		// Add vertices in counterclockwise order
		 // First triangle (bottom left, bottom right, top left)

		
		cpuGeom.verts.push_back(bl); // first triangle
		cpuGeom.verts.push_back(br);
		cpuGeom.verts.push_back(tr);

		cpuGeom.verts.push_back(tr); // 2nd triangle
		cpuGeom.verts.push_back(tl);
		cpuGeom.verts.push_back(bl);
		

		// Generates a random color for the square (I'll fix this later)
		glm::vec3 color(static_cast<float>(rand()) / RAND_MAX,
			static_cast<float>(rand()) / RAND_MAX,
			static_cast<float>(rand()) / RAND_MAX);
		for (int i = 0; i < 6; ++i) {
			cpuGeom.cols.push_back(color);
		}

		};



	// ok so basically, the top of the base square is the length of the hypotenuse for the triangle that is formed in between the two squares added to the top of the base
	// to solve for the length of the sides of the new squares being added, just do root2/2 * hypotenuse


	// Recursive function to build the Pythagoras Tree
	std::function<void(glm::vec3, glm::vec3, glm::vec3, glm::vec3, int)> add_tree = [&](glm::vec3 bl, glm::vec3 br, glm::vec3 tr, glm::vec3 tl, int depth) {
		if (depth <= 0) {
			draw_square(bl, br, tr, tl); // don't want program to crash so we'll just draw a square if this is the case
			return;
		}

		// Draw the current square
		//draw_square(bl, br, tr, tl);



		// centre of rotation: for the left square: its the bottom left
		// for the right square, its the bottom right
		//
		// NOTE: TL TR, BL, BR ARE ALL (X,Y,Z) COORDINATES
		// Calculate the size of the new squares

		float hypotenuseLength = (tr.x - tl.x);
		float new_square_side_length = hypotenuseLength * (sqrt(2.0f) / 2.0f); // s = (root2)/2 * h, this will be one side length of each of the new squares

		glm::vec3 new_bl_left, new_br_left, new_tr_left, new_tl_left;
		glm::vec3 new_bl_right, new_br_right, new_tr_right, new_tl_right;
		
		if (depth % 2 == 1) { // if depth is odd .
			// Calculate vertice sides for the left branch square (-45-degree rotation)
			new_bl_left = tl;

			new_br_left = glm::vec3(tl.x + (hypotenuseLength / 2), (tl.y + (hypotenuseLength / 2)), 0.f);

			 new_tl_left = glm::vec3(tl.x - (hypotenuseLength / 2), (tl.y + (hypotenuseLength / 2)), 0.f);

			new_tr_left= glm::vec3(tl.x, tl.y + hypotenuseLength, 0.f);

		}

		else { // if depth is even

			// Calculate vertice sides for the left branch square (note that now its already at a 45% angle, so i need to change offsets accordingly
			new_bl_left = tl; // good

			new_br_left = glm::vec3(tr.x-new_square_side_length, tl.y + new_square_side_length, 0.f); // good

			new_tl_left = glm::vec3(tl.x - new_square_side_length, tl.y, 0.f); // good

			new_tr_left = glm::vec3(tl.x - new_square_side_length, tl.y + new_square_side_length, 0.f); // good

		}

		// Draw the left branch square
		draw_square(new_bl_left, new_br_left, new_tr_left, new_tl_left);

		if (depth % 2 == 1) { // if depth is odd (right branch)

			// Calculate sides for the right branch square (+45-degree rotation)

			new_bl_right = glm::vec3(tl.x + (hypotenuseLength / 2), (tr.y + (hypotenuseLength / 2)), 0.f);

			new_br_right = tr;

			new_tl_right = glm::vec3(tr.x, tr.y + hypotenuseLength, 0.f);

			new_tr_right = glm::vec3(tr.x + (hypotenuseLength / 2), (tr.y + (hypotenuseLength / 2)), 0.f);



		}


		else {
			new_bl_right = glm::vec3(tl.x + new_square_side_length, tl.y, 0.f); // good?

			new_br_right = (tr); // good

			new_tl_right = glm::vec3(tr.x + new_square_side_length, tr.y + new_square_side_length, 0.f); // good

			new_tr_right = glm::vec3(tr.x + new_square_side_length, tr.y, 0.f); // good


		}

		// Draw the right branch square
		draw_square(new_bl_right, new_br_right, new_tr_right, new_tl_right);


		
		// Recursively add the branches
		add_tree(new_bl_left, new_br_left, new_tr_left, new_tl_left, depth - 1);
		add_tree(new_bl_right, new_br_right, new_tr_right, new_tl_right, depth - 1);

		};



	// Starting the tree with the base square
	glm::vec3 p1(-0.1f, -0.5f, 0.f); // bottom left
	glm::vec3 p2(0.1f, -0.5f, 0.f);  // bottom right
	glm::vec3 p3(0.1f, -0.3f, 0.f);  // top right
	glm::vec3 p4(-0.1f, -0.3f, 0.f); // top left

	// Draw the initial square
	draw_square(p1, p2, p3, p4);

	// Begin the recursive drawing
	add_tree(p1, p2, p3, p4, iterations);

	return cpuGeom;

}








/*
glm::vec3 calculateKochTip(const glm::vec3& A, const glm::vec3& B) {
	// Calculate the length of segment AB
	float L = glm::length(B - A);

	// Calculate point C (1/3 of the way from A to B)
	glm::vec3 C = A + (1.0f / 3.0f) * (B - A);

	// Calculate the angle of the segment AB
	float theta = atan2(B.y - A.y, B.x - A.x);

	// Calculate the tip point D
	glm::vec3 D;
	float height = (glm::sqrt(3.0f) / 6.0f) * L; // height of the equilateral triangle

	D.x = C.x + height * cos(theta + glm::pi<float>() / 3.0f); // 60 degrees in radians
	D.y = C.y + height * sin(theta + glm::pi<float>() / 3.0f);
	D.z = 0.f; // Keep the z-coordinate the same

	return D;
}
*/





//CPU_Geometry generateKochSnowflake(int depth) {
//	CPU_Geometry cpuGeom;
//
//	// Initial triangle vertices (an equilateral triangle)
//	glm::vec3 p1(-0.5f, -0.5f, 0.f);
//	glm::vec3 p2(0.5f, -0.5f, 0.f);
//	glm::vec3 p3(0.f, 0.5f, 0.f);
//
//
//	
//	// Function to draw a single triangle by pushing its vertices
//	auto draw_triangle = [&](const glm::vec3& a, const glm::vec3& b, const glm::vec3& c) {
//		cpuGeom.verts.push_back(a); // bottom left
//		cpuGeom.verts.push_back(b); // bottom right
//		cpuGeom.verts.push_back(c); // top
//
//		/*
//		// Generate a random color for each iteration, this bascially separates each triangle by a separate colour
//		glm::vec3 color(static_cast<float>(rand()) / RAND_MAX,
//			static_cast<float>(rand()) / RAND_MAX,
//			static_cast<float>(rand()) / RAND_MAX);
//			*/
//		glm::vec3 colorOne = glm::vec3(0.0f, 0.0f, 1.0f); // blue
//		glm::vec3 colorTwo = glm::vec3(0.0f, 1.0f, 0.0f); //green
//		glm::vec3 colorThree = glm::vec3(1.0f, 0.0f, 0.0f); // red
//		// Add the same color for all three vertices of the triangle
//		cpuGeom.cols.push_back(colorOne);
//		cpuGeom.cols.push_back(colorTwo);
//		cpuGeom.cols.push_back(colorThree);
//		};
//
//	
//
//	// Recursive function to subdivide the triangle 
//	std::function<void(glm::vec3, glm::vec3, glm::vec3, int)> divideSnowflake = [&](glm::vec3 a, glm::vec3 b, glm::vec3 c, int m) {
//		if (m > 0) { // where m is the iteration level
//
//			// A IS LEFT BOTTOM, B IS RIGHT BOTTOM, C IS TOP
//
//			// Calculate 1/3's of each side
//
//			// NEW TIP calculation
//			//F.x = (P1.x + P2.x) / 2 + (sqrt(3) / 6) * (P1.y - P2.y);
//			// F.y = (P1.y + P2.y) / 2 + (sqrt(3) / 6) * (P2.x - P1.x);
//
//			// LETS TRY THIS AGAIN
//
//
//			// length = sqrt((P2.x - P1.x)**2 + (P2.Y - P1.Y)**2)
//			// 
//			// NEWPEAK.X = ((P1.X + P2.X)/2) + (sqrt((P2.x - P1.x)**2 + (P2.Y - P1.Y)**2) * (sqrt(3)/2)) * cos(60)
//			// NEWPEAK.Y = ((P1.Y + P2.Y)/2) + (sqrt((P2.x - P1.x)**2 + (P2.Y - P1.Y)**2) *(sqrt(3)/2)) * sin(60)
//
//			
//
//
//			// Find points 1/3 and 2/3 the way along each side
//			glm::vec3 P1 = a + (b - a) * (1.0f / 3.0f); // 1/3 of the way from A to B
//			glm::vec3 P2 = a + (b - a) * (2.0f / 3.0f); // 2/3 of the way from A to B
//
//			// Calculate the peak of the triangle
//			//glm::vec3 P3 = P1 + glm::vec3(0.0f, (sqrt(3) / 3) * (P2.x - P1.x), 0.f); // Peak above the base // new tip to draw between A and B
//			// glm::vec3 P3(((P1.x + P2.x) / 2 + (sqrt(3) / 6) * (P1.y - P2.y)), ((P1.y + P2.y) / 2 + (sqrt(3) / 6) * (P2.x - P1.x)), 0.f);  // new tip to draw between A and B
//
//			double x_newOne = ((P1.x + P2.x) / 2) +
//				(sqrt(pow(P2.x - P1.x, 2) + pow(P2.y - P1.y, 2)) * (sqrt(3) / 2)) * cos(60);
//
//			double y_newOne = ((P1.y + P2.y) / 2) +
//				(sqrt(pow(P2.x - P1.x, 2) + pow(P2.y - P1.y, 2)) * (sqrt(3) / 2)) * sin(60);
//
//			double z_newOne = 0.0f;
//
//			glm::vec3 P3(
//				x_newOne,
//				y_newOne,
//				z_newOne
//			);
//
//
//
//			glm::vec3 P4 = b + (c - b) * (1.0f / 3.0f); // 1/3 of the way from B to C
//			glm::vec3 P5 = b + (c - b) * (2.0f / 3.0f); // 2/3 of the way from B to C
//			//glm::vec3 P6 = P4 + glm::vec3(0.0f, (sqrt(3) / 3) * (P5.x - P4.x), 0.f); // Peak above the base// new tip to draw between B and C
//			//glm::vec3 P6(((P4.x + P5.x) / 2 + (sqrt(3) / 6) * (P4.y - P5.y)), ((P4.y + P5.y) / 2 + (sqrt(3) / 6) * (P5.x - P4.x)), 0.f);// new tip to draw between B and C
//			double x_newTwo = ((P4.x + P5.x) / 2) +
//				(sqrt(pow(P5.x - P4.x, 2) + pow(P5.y - P4.y, 2)) * (sqrt(3) / 2)) * cos(60);
//
//			double y_newTwo = ((P4.y + P5.y) / 2) +
//				(sqrt(pow(P5.x - P4.x, 2) + pow(P5.y - P4.y, 2)) * (sqrt(3) / 2)) * sin(60);
//
//			double z_newTwo = 0.0f;
//
//			glm::vec3 P6(
//				x_newTwo,
//				y_newTwo,
//				z_newTwo
//			);
//
//
//			glm::vec3 P7 = c + (a - c) * (1.0f / 3.0f); // 1/3 of the way from C to A
//			glm::vec3 P8 = c + (a - c) * (2.0f / 3.0f); // 2/3 of the way from C to A
//			//glm::vec3 P9 = P7 + glm::vec3(0.0f, (sqrt(3) / 3) * (P8.x - P7.x), 0.f); // Peak above the base // new tip to draw between C and A
//			//glm::vec3 P9(((P7.x + P8.x) / 2 + (sqrt(3) / 6) * (P7.y - P8.y)), ((P7.y + P8.y) / 2 + (sqrt(3) / 6) * (P8.x - P7.x)), 0.f); // new tip to draw between C and A
//			double x_newThree = ((P7.x + P8.x) / 2) +
//				(sqrt(pow(P6.x - P7.x, 2) + pow(P8.y - P7.y, 2)) * (sqrt(3) / 2)) * cos(60);
//
//			double y_newThree = ((P7.y + P8.y) / 2) +
//				(sqrt(pow(P8.x - P7.x, 2) + pow(P8.y - P7.y, 2)) * (sqrt(3) / 2)) * sin(60);
//
//			double z_newThree = 0.0f;
//
//			glm::vec3 P9(
//				x_newThree,
//				y_newThree,
//				z_newThree
//			);
//			
//
//
//		
//			// Recursively add the triangles
//			// Recursively add the triangles
//			divideSnowflake(P1, P2, P3, m - 1);
//			divideSnowflake(P4, P5, P6, m - 1);
//			divideSnowflake(P7, P8, P9, m - 1);
//
//		}
//		else {
//			// Draw the triangle when the recursion depth is zero
//			draw_triangle(a, b, c);
//		}
//
//		};
//
//	// Start the recursive subdivision
//	divideSnowflake(p1, p2, p3, depth);
//
//
//	return cpuGeom;
//}
//
//
//


CPU_Geometry generateKochSnowflake(int depth) {
	CPU_Geometry cpuGeom;

	// Initial triangle vertices (an equilateral triangle)
	glm::vec3 p1(-0.5f, -0.5f, 0.f);
	glm::vec3 p2(0.5f, -0.5f, 0.f);
	glm::vec3 p3(0.f, 0.5f, 0.f);

	// Function to draw a single triangle by pushing its vertices
	auto draw_triangle = [&](const glm::vec3& a, const glm::vec3& b, const glm::vec3& c) {
		cpuGeom.verts.push_back(a); // bottom left
		cpuGeom.verts.push_back(b); // bottom right
		cpuGeom.verts.push_back(c); // top

		glm::vec3 colorOne = glm::vec3(0.0f, 0.0f, 1.0f); // blue
		glm::vec3 colorTwo = glm::vec3(0.0f, 1.0f, 0.0f); // green
		glm::vec3 colorThree = glm::vec3(1.0f, 0.0f, 0.0f); // red

		// Add the same color for all three vertices of the triangle
		cpuGeom.cols.push_back(colorOne);
		cpuGeom.cols.push_back(colorTwo);
		cpuGeom.cols.push_back(colorThree);
		};

	// Recursive function to subdivide the triangle
	std::function<void(glm::vec3, glm::vec3, int)> divideSnowflake = [&](glm::vec3 a, glm::vec3 b, int m) {
		if (m > 0) { // where m is the iteration level

			// Find points 1/3 and 2/3 the way along the side
			glm::vec3 P1 = a + (b - a) * (1.0f / 3.0f); // 1/3 of the way from A to B
			glm::vec3 P2 = a + (b - a) * (2.0f / 3.0f); // 2/3 of the way from A to B

			// Calculate the vector perpendicular to the side (to create the "bump")
			glm::vec3 dir = b - a;
			glm::vec3 perp(-dir.y, dir.x, 0); // Perpendicular vector in 2D

			// Normalize perpendicular vector and scale by the appropriate height
			perp = glm::normalize(perp) * (glm::length(dir) / 3.0f) * (sqrt(3.0f) / 6.0f);

			// Reverse the direction of the bump to go outward
			glm::vec3 P3 = (P1 + P2) / 2.0f - perp;

			// Recursively divide the new segments
			divideSnowflake(a, P1, m - 1); // Left segment
			divideSnowflake(P1, P3, m - 1); // Left to peak
			divideSnowflake(P3, P2, m - 1); // Peak to right
			divideSnowflake(P2, b, m - 1); // Right segment
		}
		else {
			// Draw the line when the recursion depth is zero
			cpuGeom.verts.push_back(a);
			cpuGeom.verts.push_back(b);

			// Add colors (can alternate colors or use the same for all points)
			cpuGeom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f)); // Red
			cpuGeom.cols.push_back(glm::vec3(0.0f, 1.0f, 0.0f)); // Green
		}
		};

	// Recursively divide each side of the triangle
	divideSnowflake(p1, p2, depth);
	divideSnowflake(p2, p3, depth);
	divideSnowflake(p3, p1, depth);

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
			else if (currentFractalType == 3) {
				cpuGeom = generateKochSnowflake(p.iterations);
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
