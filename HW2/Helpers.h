#ifndef __HELPERS_H__
#define __HELPERS_H__

#define ABS(a) ((a) > 0 ? (a) : -1 * (a))
#define EPSILON 0.000000001
#define PI 3.14159265

#include "Matrix4.h"
#include "Vec3.h"
#include "Vec4.h"

/*
 * Calculate cross product of vec3 a, vec3 b and return resulting vec3.
 */
Vec3 crossProductVec3(Vec3 a, Vec3 b);

/*
 * Calculate dot product of vec3 a, vec3 b and return resulting value.
 */
double dotProductVec3(Vec3 a, Vec3 b);

/*
 * Find length (|v|) of vec3 v.
 */
double magnitudeOfVec3(Vec3 v);

/*
 * Normalize the vec3 to make it unit vec3.
 */
Vec3 normalizeVec3(Vec3 v);

/*
 * Return -v (inverse of vec3 v)
 */
Vec3 inverseVec3(Vec3 v);

/*
 * Add vec3 a to vec3 b and return resulting vec3 (a+b).
 */
Vec3 addVec3(Vec3 a, Vec3 b);

/*
 * Subtract vec3 b from vec3 a and return resulting vec3 (a-b).
 */
Vec3 subtractVec3(Vec3 a, Vec3 b);

/*
 * Multiply each element of vec3 with scalar.
 */
Vec3 multiplyVec3WithScalar(Vec3 v, double c);

/*
 * Prints elements in a vec3. Can be used for debugging purposes.
 */
void printVec3(Vec3 v);

/*
 * Check whether vec3 a and vec3 b are equal.
 * In case of equality, returns 1.
 * Otherwise, returns 0.
 */
int areEqualVec3(Vec3 a, Vec3 b);
/*
 * Returns an identity matrix (values on the diagonal are 1, others are 0).
*/
Matrix4 getIdentityMatrix();

/*
 * Multiply matrices m1 (Matrix4) and m2 (Matrix4) and return the result matrix r (Matrix4).
 */
Matrix4 multiplyMatrixWithMatrix(Matrix4 m1, Matrix4 m2);

/*
 * Multiply matrix m (Matrix4) with vector v (vec4) and store the result in vector r (vec4).
 */
Vec4 multiplyMatrixWithVec4(Matrix4 m, Vec4 v);

/*
* Returns a translation transformation matrix moving the point by x, y, z
*/
Matrix4 translationMatrix(double x, double y, double z);

/*
* Returns a rotation transformation matrix rotating the point around x axis by angle
*/
Matrix4 rotationxMatrix(double angle);

/*
* Returns a rotation transformation matrix rotating the point around y axis by angle
*/
Matrix4 rotationyMatrix(double angle);

/*
* Returns a rotation transformation matrix rotating the point around x axis by angle
*/
Matrix4 rotationzMatrix(double angle);

/*
* Returns a rotation transformation matrix rotating the point around (0, 0, 0), (x, y, z) axis by angle
*/
Matrix4 rotationMatrix(double angle, double x, double y, double z);

/*
* Returns a scaling transformation matrix scaling the point by x, y, z
*/
Matrix4 scalingMatrix(double x, double y, double z);

/*
* Return camera transformation matrix which sets camera pos as (0, 0, 0) and u, v, w ad major axises
*/
Matrix4 cameraTransformMatrix(Vec3 pos, Vec3 u, Vec3 v, Vec3 w);

/*
*  Retrurn orthographic projection matrix with given near and dar planes.
*/
Matrix4 orthographicProjectionMatrix(double left, double right, double bottom, double top, double near, double far);

/*
*  Return p2o matriz which turns a ortographic projection matrix into perspective projection matrix.
*/
Matrix4 p2OProjectionMatrix(double far, double near);

/*
*  Retrurn perspective projection matrix with given near and dar planes.
*/
Matrix4 perspectiveProjectionMatrix(double left, double right, double bottom, double top, double near, double far);

Matrix4 viewPortMatrix(int x, int y);

Vec3 perspectiveDivide(Vec4 vec);

Vec3 transformVertice(Vec3 vec, Matrix4 matrix);

Vec3 getNormal(Vec3 a, Vec3 b, Vec3 c);

Vec3 barycentricToCartesian(Vec3 ver1, Vec3 ver2, Vec3 ver3, double alpha, double beta);

#endif