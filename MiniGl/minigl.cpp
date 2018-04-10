/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>
#include <stack>

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h

class Pixel;
class Matrix;
struct vertex;
class BoundingBox;
struct triangle;

MGLpoly_mode draw_mode;
vector<vertex> vectorlist;
vec3 current_color;
vector<triangle> trianglelist;
vector<Pixel> frameBuffer;
vector<Pixel> zBuffer;
bool failure = false;
MGLmatrix_mode matrix_mode;
stack <Matrix> modelMatrixS;
stack <Matrix> projMatrixS;

class Pixel
{
	public:
	int x,y;
	MGLpixel pColor;
	MGLfloat z;
	
	Pixel(int X, int Y, MGLpixel c, MGLfloat Z) : x(X), y(Y), pColor(c), z(Z) {}
};

class Matrix
{
	// variables
	public:
	MGLfloat matrix[4][4];
	
	Matrix()
	{
		clearMatrix();
		initMatrix (0, 0, 0);
	}
	
	Matrix(MGLfloat X, MGLfloat Y, MGLfloat Z)
	{
		clearMatrix();
		initMatrix (X, Y, Z);
	}
	
	Matrix& operator= (const Matrix& rhs)
	{
		if (this != &rhs) 
		{
			for (int row = 0; row < 4; ++row)
				for (int col = 0; col < 4; ++col)
					matrix[row][col] = rhs.matrix[row][col];
		}
		
		return *this;
	}
	
	// Currently unused. mglMultMatrix is used instead.
	Matrix operator* (const Matrix& rhs)
	{
		Matrix result;
		result.clearMatrix();
		
		MGLfloat sum = 0;
		
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				sum = 0;
				for (int k = 0; k < 4; ++k)
					sum += matrix[i][k] * rhs.matrix[k][j];
				result.matrix[i][j] = sum;
			}
		}
			
		return result;
	}
	
	friend ostream& operator<< (ostream& os, const Matrix& m)
	{
		for (int row = 0; row < 4; ++row)
		{
			for (int col = 0; col < 4; ++col)
				os << m.matrix[row][col] << " ";
			os << "\n";
		}
		return os;
		
	}
	
	void clearMatrix()
	{
		for (int row = 0; row < 4; ++row)
			for (int col = 0; col < 4; ++col)
				matrix[row][col] = 0;
	}
	
	void createScaler(float x, float y, float z)
	{
		clearMatrix();
		matrix[0][0] = x;
		matrix[1][1] = y;
		matrix[2][2] = z;
		matrix[3][3] = 1;
	}
	
	void createTranslater(MGLfloat X, MGLfloat Y, MGLfloat Z)
	{
		matrix[0][0] = 1;
		matrix[1][1] = 1;
		matrix[2][2] = 1;
		matrix[3][3] = 1;
		
		matrix[0][3] = X;
		matrix[1][3] = Y;
		matrix[2][3] = Z;
	}	
	
	private:
	void initMatrix (MGLfloat X, MGLfloat Y, MGLfloat Z)
	{	
		// set w
		matrix[3][3] = 1;
		
		// set the x, y, z coordinate
		matrix[0][3] = X; 
		matrix[1][3] = Y; 
		matrix[2][3] = Z;
	}
};

struct vertex {
  vec3 color;
  vec4 pos;
  vec4 xyzw_Screen;
  MGLpixel vColor[3];
  
  friend ostream& operator<< (ostream& os, const vertex& v)
	{
		os << "( " << v.pos[0] << ", " << v.pos[1] << ", " << v.pos[2] << ", " << v.pos[3] << " ) " << endl;
		return os;
		
	}
	
	vertex& operator= (const vertex& rhs)
	{
		if (this != &rhs) 
		{
			pos[0] = rhs.pos[0];
			pos[1] = rhs.pos[1];
			pos[2] = rhs.pos[2];
			pos[3] = rhs.pos[3];
			color[0] = rhs.color[0];
			color[1] = rhs.color[1];
			color[2] = rhs.color[2];
			
			xyzw_Screen[0] = rhs.xyzw_Screen[0];
			xyzw_Screen[1] = rhs.xyzw_Screen[1];
			xyzw_Screen[2] = rhs.xyzw_Screen[2];
			xyzw_Screen[3] = rhs.xyzw_Screen[3];
		}
		
		return *this;
	}
	
	/* This operation performs matrix and vector multiplication.
	 * The vector must be on the left of *, and the matrix must be on
	 * the right.
	 * 
	 * It follows the matrix * vector order of multiplication.
	 */
	vertex operator* (const Matrix& rhs)
	{
		MGLfloat currVertex[4] = {pos[0], pos[1], pos[2], pos[3]};
		MGLfloat newVertex[4] = {0, 0, 0, 0};
		
		MGLfloat sum, val = 0;
		
		// [ Vertex ] [ Matrix ]
		for (int row = 0; row < 4; ++row)
		{
			sum = 0;
			for (int col = 0; col < 4; ++col)
			{
				val = currVertex[col];
				sum += rhs.matrix[row][col] * val;	
			}
			
			newVertex[row] = sum;
			//cout << "row: " << row << " sum: " << sum << endl;
		}
		
		vertex v;
		v.pos =  vec4(newVertex[0], newVertex[1], newVertex[2], newVertex[3]);
		return v;
}

void ndclipping(MGLfloat & a)
	{
		if(a < -1)
			a = -1;
		if(a > 1)
			a = 1;
		
	}
  
  void vscaletoscreen(MGLsize width, MGLsize height)
  
  {
	  pos[0] = pos[0]/pos[3];
	  pos[1] = pos[1]/pos[3];
	  pos[2] = pos[2]/pos[3];
//	  pos[3] = pos[3]/pos[3];
	  
//	  ndclipping(pos[0]);
//	  ndclipping(pos[1]);
//	  ndclipping(pos[2]);
	

	 MGLfloat  i = ((pos[0] + 1) * width) /2 -0.5;
	 MGLfloat  j = ((pos[1] + 1) * height) /2 - 0.5;

	  if (i > width) 
		i = width;

	  xyzw_Screen[0] = i;
	  xyzw_Screen[1] = j;
	  xyzw_Screen[2] = pos[2];
	  xyzw_Screen[3] = pos[3];
	  cout << "TESTS:    " << xyzw_Screen[0] << " " << xyzw_Screen[1] << " " << xyzw_Screen[2] << " " << xyzw_Screen[3] << endl;
  }
  
  	void applyTransformations()
	{
		Matrix model = modelMatrixS.top();
		Matrix proj = projMatrixS.top();
		Matrix trans;
		trans.createTranslater(1, 1, 1);

		vertex v;
		v.pos = pos;		
		v = v * model;
		v = v * proj;	
		//v = v * trans;	
		
			
		// update the x, y, z, w values.
		pos[0] = v.pos[0];
		pos[1] = v.pos[1];
		pos[2] = v.pos[2];
		pos[3] = v.pos[3];
}		
};

struct triangle {
  vertex a,b,c;
  
  void scaletoscreen(MGLsize width, MGLsize height)
  {
	  a.vscaletoscreen(width, height);
	  b.vscaletoscreen(width, height);
	  c.vscaletoscreen(width, height);
	  
  }
};

class BoundingBox
{
public:
	float min_x, max_x, min_y, max_y;

	BoundingBox()
		: min_x(0), max_x(0), min_y(0), max_y(0)
	{}

	void initBB(const triangle& t)
	{
		min_x = getMin_X(t);
		min_y = getMin_Y(t);
		max_x = getMax_X(t);
		max_y = getMax_Y(t);
	}
	
	void initBB(const triangle& t1, const triangle& t2)
	{
		min_x = getMin_X(t1, t2.c);
		min_y = getMin_Y(t1, t2.c);
		max_x = getMax_X(t1, t2.c);
		max_y = getMax_Y(t1, t2.c);
	}

private:
	// returns the minimum value 
	float getMin_X(const triangle& t)
	{
		float min = t.a.xyzw_Screen[0];

		if (t.b.xyzw_Screen[0] < min) min = t.b.xyzw_Screen[0];
		if (t.c.xyzw_Screen[0] < min) min = t.c.xyzw_Screen[0];
		return min;
	}

	float getMin_Y(const triangle& t)
	{
		float min = t.a.xyzw_Screen[1];

		if (t.b.xyzw_Screen[1] < min) min = t.b.xyzw_Screen[1];
		if (t.c.xyzw_Screen[1] < min) min = t.c.xyzw_Screen[1];

		return min;
	}

	float getMax_X(const triangle& t)
	{
		float max = t.a.xyzw_Screen[0];

		if (t.b.xyzw_Screen[0] > max) max = t.b.xyzw_Screen[0];
		if (t.c.xyzw_Screen[0] > max) max = t.c.xyzw_Screen[0];

		return max;
	}

	float getMax_Y(const triangle& t)
	{
		float max = t.a.xyzw_Screen[1];

		if (t.b.xyzw_Screen[1] > max) max = t.b.xyzw_Screen[1];
		if (t.c.xyzw_Screen[1] > max) max = t.c.xyzw_Screen[1];

		return max;
	}
	
	float getMin_X(const triangle& t, const vertex& v)
	{
		float min = t.a.xyzw_Screen[0];

		if (t.b.xyzw_Screen[0] < min) min = t.b.xyzw_Screen[0];
		if (t.c.xyzw_Screen[0] < min) min = t.c.xyzw_Screen[0];
		if (v.xyzw_Screen[0] < min) min = v.xyzw_Screen[0];
		
		return min;
	}

	float getMin_Y(const triangle& t, const vertex& v)
	{
		float min = t.a.xyzw_Screen[1];

		if (t.b.xyzw_Screen[1] < min) min = t.b.xyzw_Screen[1];
		if (t.c.xyzw_Screen[1] < min) min = t.c.xyzw_Screen[1];
		if (v.xyzw_Screen[1] < min) min = v.xyzw_Screen[1];

		return min;
	}

	float getMax_X(const triangle& t, const vertex& v)
	{
		float max = t.a.xyzw_Screen[0];

		if (t.b.xyzw_Screen[0] > max) max = t.b.xyzw_Screen[0];
		if (t.c.xyzw_Screen[0] > max) max = t.c.xyzw_Screen[0];
		if (v.xyzw_Screen[0] > max) max = v.xyzw_Screen[0];

		return max;
	}

	float getMax_Y(const triangle& t , const vertex& v)
	{
		float max = t.a.xyzw_Screen[1];

		if (t.b.xyzw_Screen[1] > max) max = t.b.xyzw_Screen[1];
		if (t.c.xyzw_Screen[1] > max) max = t.c.xyzw_Screen[1];
		if (v.xyzw_Screen[1] > max) max = v.xyzw_Screen[1];

		return max;
	}

};



/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}

void set_pixel(int x, int y, MGLpixel c, MGLfloat z)
{
	if(z > -1 && z < 1)
	{
	Pixel pix(x,y,c,z);
	frameBuffer.push_back(pix);
	zBuffer.push_back(pix);
	}
}

void conv2screencoord(MGLsize width, MGLsize height, triangle& t)
{
	
		t.scaletoscreen(width, height);
	
}

float bary(float x, float y, float x_b, float y_b, float x_c, float y_c)
{
	return (y_b - y_c)*x + (x_c - x_b)*y + ((x_b*y_c) - (x_c*y_b));
}

bool sortByZ(Pixel a, Pixel b) { return a.z > b.z; }

MGLpixel mixColors(float alpha, float beta, float gamma, const vec3 vColor1, const vec3 vColor2, const vec3 vColor3)
{
	
		float newR = alpha*vColor1[0] + beta*vColor2[0] + gamma*vColor3[0];
		float newG = alpha*vColor1[1] + beta*vColor2[1] + gamma*vColor3[1];
		float newB = alpha*vColor1[2] + beta*vColor2[2] + gamma*vColor3[2];
		MGLpixel newColor;
		float r = 255* newR;
		float g = 255* newG;
		float b = 255* newB;
		newColor = Make_Pixel(r,g,b);
//		newColor = Make_Pixel(255, 255, 255);
		return newColor;
}

void drawTriangle(int x, int y, const vertex& a, const vertex& b, const vertex& c)
{
	// vertex A
	float x_a = a.xyzw_Screen[0];
	float y_a = a.xyzw_Screen[1];
	// vertex B
	float x_b = b.xyzw_Screen[0];
	float y_b = b.xyzw_Screen[1];

	// vertex C 
	float x_c = c.xyzw_Screen[0];
	float y_c = c.xyzw_Screen[1];
	

	
//	MGLfloat fabc = ((y_a - y_b)*x_c) + ((x_b - x_a)*y_c) + (x_a*y_b - x_b*y_a);
//  MGLfloat facb = ((y_c - y_a)*x_b) + ((x_a - x_c)*y_c) + (x_c*y_a - x_a*y_c);
//	MGLfloat fbca = ((y_b - y_c)*x_a) + ((x_c - x_b)*y_a) + (x_b*y_c - x_c*y_b);
	
	
//  MGLfloat fcar = ((y_c - y_a)*(80)) + ((x_a - x_c)*(60)) + (x_c*y_a - x_a*y_c);
//  MGLfloat fabr = ((y_a - y_b)*(80)) + ((x_b - x_a)*(60)) + (x_a*y_b - x_b*y_a);
//	MGLfloat fbcr = ((y_b - y_c)*(80)) + ((x_c - x_b)*(60)) + (x_b*y_c - x_c*y_b);

//	float alpha = bary(x, y, x_b, y_b, x_c, y_c) / bary(x_a, y_a, x_b, y_b, x_c, y_c);
//	float beta = bary(x, y, x_c, y_c, x_a, y_a) / bary(x_b, y_b, x_c, y_c, x_a, y_a);
//	float gamma = bary(x, y, x_a, y_a, x_b, y_b) / bary(x_c, y_c, x_a, y_a, x_b, y_b);


	
	float alpha = bary(x, y, x_b, y_b, x_c, y_c) / bary(x_a, y_a, x_b, y_b, x_c, y_c);
	float beta = bary(x_a, y_a, x, y, x_c, y_c) / bary(x_a, y_a, x_b, y_b, x_c, y_c);
	float gamma = bary(x_a, y_a, x_b, y_b, x, y) / bary(x_a, y_a, x_b, y_b, x_c, y_c);
	
//	bool one = (0 < alpha);
//    bool two = (0 < beta);
//   bool three = (0 < gamma);
//    bool four = (fbca*fbcr > 0);
//   bool five = (facb*fcar > 0);
//	bool six = (fabc*fabr > 0);

//	alpha = (alpha/a.xyzw_Screen[3])/((alpha/a.xyzw_Screen[3]) + (beta/b.xyzw_Screen[3]) + (gamma/c.xyzw_Screen[3]));
//	beta = (beta/b.xyzw_Screen[3])/((alpha/a.xyzw_Screen[3]) + (beta/b.xyzw_Screen[3]) + (gamma/c.xyzw_Screen[3]));
//	gamma = (gamma/c.xyzw_Screen[3])/((alpha/a.xyzw_Screen[3]) + (beta/b.xyzw_Screen[3]) + (gamma/c.xyzw_Screen[3]));
	
	// Inside the triangle
	if (alpha >= 0 && beta >= 0 && gamma >= 0)
	{

		
		MGLpixel pixelColor = mixColors(alpha, beta, gamma, a.color, b.color, c.color);
		MGLfloat newZ = alpha*a.xyzw_Screen[2] + beta*b.xyzw_Screen[2] + gamma*c.xyzw_Screen[2];
		
		set_pixel(x,y, pixelColor, newZ);
	}
}

void clipping(int & min, int & max, MGLsize upper)
{
	int u = upper;
	if (min < 0)
		min = 0;
	if (max > u)
		max = upper;
		
	cout << min << " " << max << endl;
}

void rasterizeT(MGLsize width, MGLsize height, const triangle& t)
{
	BoundingBox MGL_BoundingBox;

	MGL_BoundingBox.initBB(t);


	int x_min = static_cast<MGLpixel>(MGL_BoundingBox.min_x + 0.5);
	int x_max = static_cast<MGLpixel>(MGL_BoundingBox.max_x + 0.5);
	int y_min = static_cast<MGLpixel>(MGL_BoundingBox.min_y + 0.5);
	int y_max = static_cast<MGLpixel>(MGL_BoundingBox.max_y + 0.5);
	
	clipping(x_min, x_max, width);
	clipping(y_min, y_max, height);
	
	for (int i = x_min; i <= x_max; ++i)
		for (int j = y_min; j <= y_max; ++j)
			drawTriangle(i,j, t.a, t.b, t.c);
}

void rasterizeQ(MGLsize width, MGLsize height, const triangle& t1, const triangle& t2) 
{
	BoundingBox MGL_BoundingBox;
//	cout << t1.a.pos[0] << " " << t1.a.pos[1] << t1.b.pos[0] << " " << t1.b.pos[1] << t1.c.pos[0] << " " << t1.c.pos[1] << t2.c.pos[0] << " " << t2.c.pos[1] << endl;

	MGL_BoundingBox.initBB(t1, t2);


	int x_min = static_cast<MGLpixel>(MGL_BoundingBox.min_x + 0.5);
	int x_max = static_cast<MGLpixel>(MGL_BoundingBox.max_x + 0.5);
	int y_min = static_cast<MGLpixel>(MGL_BoundingBox.min_y + 0.5);
	int y_max = static_cast<MGLpixel>(MGL_BoundingBox.max_y + 0.5);
	
	clipping(x_min, x_max, width);
	clipping(y_min, y_max, height);
	
	for (int i = x_min; i <= x_max; ++i)
	{
			for (int j = y_min; j <= y_max; ++j)
		{
			drawTriangle(i,j, t1.a, t1.b, t1.c);
			drawTriangle(i,j, t2.a, t2.b, t2.c);
		}
	}
}
/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{	
	int numTriangles = trianglelist.size();
	for (int i = 0; i < numTriangles; i++)
	{
		if (draw_mode == MGL_TRIANGLES)
		{
			conv2screencoord(width, height, trianglelist.at(i));
			rasterizeT(width, height, trianglelist.at(i));
		}
		else if (draw_mode == MGL_QUADS)
		{
			conv2screencoord(width, height, trianglelist.at(i));
			conv2screencoord(width, height, trianglelist.at(i+1));
			rasterizeQ(width, height, trianglelist.at(i), trianglelist.at(i+1));
			i+=2;
		}
	}
	
		// sort the zBuffer by descending order.
	sort(zBuffer.begin(), zBuffer.end(), sortByZ);

	int size = zBuffer.size();
	 
	 for (int i = 0; i < size; ++i)
	 {
		int x = zBuffer[i].x;
		int y = zBuffer[i].y;
		
		MGLpixel color = zBuffer[i].pColor;
		*(data + y*width + x) = color;
		
	}
	



	trianglelist.clear();
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
  if (mode == MGL_TRIANGLES || mode == MGL_QUADS)
	 	draw_mode = mode;
  else
    failure = true;
}


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
  if (draw_mode == MGL_TRIANGLES && vectorlist.size()%3 == 0)
	{
      triangle t;
      for(unsigned int i = 0; i < vectorlist.size(); i = i + 3)
	{
     t.a = vectorlist.at(i);
	 t.b = vectorlist.at(i+1);
	 t.c = vectorlist.at(i+2);
	 trianglelist.push_back(t);
	}

    }
  if (draw_mode == MGL_QUADS && vectorlist.size()%4 == 0)
    {
      triangle t1;
      triangle t2;
	for (unsigned int i = 0; i < vectorlist.size(); i = i + 4)
	  {
	    t1.a = vectorlist.at(i);
	    t1.b = vectorlist.at(i+1);
	    t1.c = vectorlist.at(i+2);
	    t2.a = vectorlist.at(i);
	    t2.b = vectorlist.at(i+2);
	    t2.c = vectorlist.at(i+3);
	    
	    trianglelist.push_back(t1);
	    trianglelist.push_back(t2);
	  }
    }
  vectorlist.clear();
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x, MGLfloat y)
{
  mglVertex3(x,y,0);
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x, MGLfloat y, MGLfloat z)
{
	
	if (failure == true)
	{
		MGL_ERROR("Failure Aborting...Aborting! \n");
	}
	vertex v;
	v.color = current_color;
	v.pos = vec4(x,y,z,1);
//	cout << v.pos[0] << " " << v.pos[1] << " " << v.pos[2] << " " << v.pos[3] << endl;
	v.applyTransformations();
//	cout << v.pos[0] << " " << v.pos[1] << " " << v.pos[2] << " " << v.pos[3] << endl;
	vectorlist.push_back(v);
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
	
	if (mode == MGL_MODELVIEW || mode == MGL_PROJECTION)
		matrix_mode = mode;
	else
		failure = true;

}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
	if (matrix_mode == MGL_MODELVIEW && !modelMatrixS.empty())
		modelMatrixS.push(modelMatrixS.top());
	else if (matrix_mode == MGL_PROJECTION && !projMatrixS.empty());
		projMatrixS.push(projMatrixS.top());
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
	if (matrix_mode == MGL_MODELVIEW && !modelMatrixS.empty())
	{
		if (!modelMatrixS.empty())
			modelMatrixS.pop();
	}
	
	else if (matrix_mode == MGL_PROJECTION && !projMatrixS.empty())
	{
		if (!projMatrixS.empty())
			projMatrixS.pop();
	}
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
		
	Matrix Identity;
	Identity.matrix[0][0] = 1;
	Identity.matrix[1][1] = 1;
	Identity.matrix[2][2] = 1;
	
	
	if (matrix_mode == MGL_PROJECTION)
	{

		if (!projMatrixS.empty())
			projMatrixS.pop();
			
		projMatrixS.push(Identity);
		
	}
	else if (matrix_mode == MGL_MODELVIEW)
	{
		if (!modelMatrixS.empty())
			modelMatrixS.pop();
			
		modelMatrixS.push(Identity);
	}

}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
	Matrix load;
	for (int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
			{
				load.matrix[i][j] = matrix[j * 4 + i];
			}
	}
	
	if (matrix_mode == MGL_MODELVIEW)
	{
		if(!modelMatrixS.empty())
			modelMatrixS.pop();
		modelMatrixS.push(load);
	}
	else if (matrix_mode == MGL_PROJECTION)
	{
		if(!projMatrixS.empty())
			projMatrixS.pop();
		projMatrixS.push(load);
	}
	
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
	Matrix result;
	Matrix left;
	left.clearMatrix();
	result.clearMatrix();
	
	
	if (matrix_mode == MGL_MODELVIEW)
	{
		left = modelMatrixS.top();
				
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
		
				for (int k = 0; k < 4; ++k)
				{
				result.matrix[i][j] += left.matrix[i][k] * matrix[j * 4 + k];
				
				}
		
			}
		}
		if (!modelMatrixS.empty())
			modelMatrixS.pop();
		
		modelMatrixS.push(result);
	}
	else if (matrix_mode == MGL_PROJECTION)
	{
		left = projMatrixS.top();
		
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
		
				for (int k = 0; k < 4; ++k)
				{	
					result.matrix[i][j] += left.matrix[i][k] * matrix[j * 4 + k];
					
				}
		
			}
		}
		if(!projMatrixS.empty())
			projMatrixS.pop();
		
		projMatrixS.push(result);
	}
}


void mglMultMatrix(Matrix& left, const Matrix& m)
{
	Matrix result;
	result.clearMatrix();
	
	// [ currMatrix ] [ m ]
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
		
			for (int k = 0; k < 4; ++k)
			{
				result.matrix[i][j] += left.matrix[i][k] * m.matrix[k][j];
				
			}
		
		}
	}
	left = result;
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x, MGLfloat y, MGLfloat z)
{
	Matrix m;
	m.createTranslater(x,y,z);
	
	if (matrix_mode == MGL_PROJECTION)
		mglMultMatrix(projMatrixS.top(), m); 
		
	else if (matrix_mode == MGL_MODELVIEW)
		mglMultMatrix(modelMatrixS.top(), m); 
}

void normalize(MGLfloat &x, MGLfloat &y, MGLfloat &z)
{
	MGLfloat magnitude = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2) );
		
	if (magnitude == 0)
		return;
		
	x = x/magnitude;
	y = y/magnitude;
	z = z/magnitude; 
}
/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle, MGLfloat x, MGLfloat y, MGLfloat z)
{
	angle = angle * (M_PI/180);
	MGLfloat s = sin(angle);
	MGLfloat c = cos(angle);	
	
	normalize(x, y, z);
	
	MGLfloat a = (x*x) * (1-c) + c;
	MGLfloat b = (x*y) * (1-c) - z*s;
	MGLfloat cr = (x*z) * (1-c) + y*s;
	
	MGLfloat d = (y*x) * (1-c) + z*s;
	MGLfloat e = (y*y) * (1-c) + c;
	MGLfloat f = (y*z) * (1-c) - x*s;
	
	MGLfloat g = (x*z) * (1-c) - y*s;
	MGLfloat h = (y*z) * (1-c) + x*s;
	MGLfloat i = (z*z) * (1-c) + c;
	
	Matrix r(0, 0, 0);
	
	r.matrix[0][0] = a;
	r.matrix[0][1] = b;
	r.matrix[0][2] = cr;
	
	r.matrix[1][0] = d;
	r.matrix[1][1] = e;
	r.matrix[1][2] = f;
	
	r.matrix[2][0] = g;
	r.matrix[2][1] = h;
	r.matrix[2][2] = i;

	if (matrix_mode == MGL_PROJECTION)
		mglMultMatrix(projMatrixS.top(), r); 
		
	else if (matrix_mode == MGL_MODELVIEW)
		mglMultMatrix(modelMatrixS.top(), r); 
}


/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x, MGLfloat y, MGLfloat z)
{
	Matrix m;
	
	m.matrix[0][0] = x;
	m.matrix[1][1] = y;
	m.matrix[2][2] = z;
	
	if (matrix_mode == MGL_PROJECTION)
		mglMultMatrix(projMatrixS.top(), m); 
		
	else if (matrix_mode == MGL_MODELVIEW)
		mglMultMatrix(modelMatrixS.top(), m);  
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
	MGLfloat x,y,A,B,C,D;
	
	x = (near * 2.0f)/(right - left);
	y = (near * 2.0f)/(top - bottom);
	A = (right + left)/(right - left);
	B = (top + bottom)/(top - bottom);
	C = -(far + near)/(far - near);
	D = -(2.0f * far * near)/(far - near);
	
	Matrix frustrum(0, 0, D);
	
	frustrum.matrix[0][0] = x;
	frustrum.matrix[1][1] = y;
	
	frustrum.matrix[0][2] = A;
	frustrum.matrix[1][2] = B;
	frustrum.matrix[2][2] = C;
	frustrum.matrix[3][2] = -1;
	frustrum.matrix[3][3] = 0;
	
	if (matrix_mode == MGL_PROJECTION)
		mglMultMatrix(projMatrixS.top(), frustrum); 
		
	else if (matrix_mode == MGL_MODELVIEW)
		mglMultMatrix(modelMatrixS.top(), frustrum); 
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
 
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
	MGLfloat t_x, t_y, t_z, x, y, z;
	x = 2/(right - left);
	y = 2/(top - bottom);
	z = -2/(far - near);

	t_x = -(right + left)/(right - left);
	t_y = -(top + bottom)/(top-bottom);
	t_z = -(far + near)/(far - near);

	Matrix ortho;
	ortho.matrix[0][0] = x;
	ortho.matrix[1][1] = y;
	ortho.matrix[2][2] = z;
	
	ortho.matrix[0][3] = t_x;
	ortho.matrix[1][3] = t_y;
	ortho.matrix[2][3] = t_z;
	
	if (matrix_mode == MGL_PROJECTION)
		mglMultMatrix(projMatrixS.top(), ortho); 
		
	else if (matrix_mode == MGL_MODELVIEW)
		mglMultMatrix(modelMatrixS.top(), ortho); 
		
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
	current_color = vec3(red,green,blue);
}


