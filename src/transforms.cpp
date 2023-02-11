#include "transforms.h"

#include "CGL/matrix3x3.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"
#include <math.h>

namespace CGL {

Vector2D operator*(const Matrix3x3 &m, const Vector2D &v) {
	Vector3D mv = m * Vector3D(v.x, v.y, 1);
	return Vector2D(mv.x / mv.z, mv.y / mv.z);
}

Matrix3x3 translate(float dx, float dy) {
	// Part 3: Fill this in.
	Matrix3x3 A;

	A(0, 0) = 1;
	A(0, 1) = 0;
	A(0, 2) = dx;
	A(1, 0) = 0;
	A(1, 1) = 1;
	A(1, 2) = dy;
	A(2, 0) = 0;
	A(2, 1) = 0;
	A(2, 2) = 1;

	return A;
}

Matrix3x3 scale(float sx, float sy) {
	// Part 3: Fill this in.
	Matrix3x3 A;
	A(0, 0) = sx;
	A(0, 1) = 0;
	A(0, 2) = 0;
	A(1, 0) = 0;
	A(1, 1) = sy;
	A(1, 2) = 0;
	A(2, 0) = 0;
	A(2, 1) = 0;
	A(2, 2) = 1;

	return A;
}

// The input argument is in degrees counterclockwise
Matrix3x3 rotate(float deg) {
	// Part 3: Fill this in.
	Matrix3x3 A;

	A(0, 0) = cos(radians(deg));
	A(0, 1) = -sin(radians(deg));
	A(0, 2) = 0;
	A(1, 0) = sin(radians(deg));
	A(1, 1) = cos(radians(deg));
	A(1, 2) = 0;
	A(2, 0) = 0;
	A(2, 1) = 0;
	A(2, 2) = 1;

	return A;
	}

}
