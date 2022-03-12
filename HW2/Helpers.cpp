#include <iostream>
#include <cmath>
#include "Helpers.h"
#include "Matrix4.h"
#include "Vec3.h"
#include "Vec4.h"

using namespace std;

/*
 * Calculate cross product of vec3 a, vec3 b and return resulting vec3.
 */
Vec3 crossProductVec3(Vec3 a, Vec3 b)
{
    Vec3 result;

    result.x = a.y * b.z - b.y * a.z;
    result.y = b.x * a.z - a.x * b.z;
    result.z = a.x * b.y - b.x * a.y;

    return result;
}

/*
 * Calculate dot product of vec3 a, vec3 b and return resulting value.
 */
double dotProductVec3(Vec3 a, Vec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/*
 * Find length (|v|) of vec3 v.
 */
double magnitudeOfVec3(Vec3 v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

/*
 * Normalize the vec3 to make it unit vec3.
 */
Vec3 normalizeVec3(Vec3 v)
{
    Vec3 result;
    double d;

    d = magnitudeOfVec3(v);
    result.x = v.x / d;
    result.y = v.y / d;
    result.z = v.z / d;

    return result;
}

/*
 * Return -v (inverse of vec3 v)
 */
Vec3 inverseVec3(Vec3 v)
{
    Vec3 result;
    result.x = -v.x;
    result.y = -v.y;
    result.z = -v.z;

    return result;
}

/*
 * Add vec3 a to vec3 b and return resulting vec3 (a+b).
 */
Vec3 addVec3(Vec3 a, Vec3 b)
{
    Vec3 result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;

    return result;
}

/*
 * Subtract vec3 b from vec3 a and return resulting vec3 (a-b).
 */
Vec3 subtractVec3(Vec3 a, Vec3 b)
{
    Vec3 result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;

    return result;
}

/*
 * Multiply each element of vec3 with scalar.
 */
Vec3 multiplyVec3WithScalar(Vec3 v, double c)
{
    Vec3 result;
    result.x = v.x * c;
    result.y = v.y * c;
    result.z = v.z * c;

    return result;
}

/*
 * Prints elements in a vec3. Can be used for debugging purposes.
 */
void printVec3(Vec3 v)
{
    cout << "(" << v.x << "," << v.y << "," << v.z << ")" << endl;
}

/*
 * Check whether vec3 a and vec3 b are equal.
 * In case of equality, returns 1.
 * Otherwise, returns 0.
 */
int areEqualVec3(Vec3 a, Vec3 b)
{

    /* if x difference, y difference and z difference is smaller than threshold, then they are equal */
    if ((ABS((a.x - b.x)) < EPSILON) && (ABS((a.y - b.y)) < EPSILON) && (ABS((a.z - b.z)) < EPSILON))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

/*
 * Returns an identity matrix (values on the diagonal are 1, others are 0).
*/
Matrix4 getIdentityMatrix()
{
    Matrix4 result;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (i == j)
            {
                result.val[i][j] = 1.0;
            }
            else
            {
                result.val[i][j] = 0.0;
            }
        }
    }

    return result;
}

/*
 * Multiply matrices m1 (Matrix4) and m2 (Matrix4) and return the result matrix r (Matrix4).
 */
Matrix4 multiplyMatrixWithMatrix(Matrix4 m1, Matrix4 m2)
{
    Matrix4 result;
    double total;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            total = 0;
            for (int k = 0; k < 4; k++)
            {
                total += m1.val[i][k] * m2.val[k][j];
            }

            result.val[i][j] = total;
        }
    }

    return result;
}

/*
 * Multiply matrix m (Matrix4) with vector v (vec4) and store the result in vector r (vec4).
 */
Vec4 multiplyMatrixWithVec4(Matrix4 m, Vec4 v)
{
    double values[4];
    double total;

    for (int i = 0; i < 4; i++)
    {
        total = 0;
        for (int j = 0; j < 4; j++)
        {
            total += m.val[i][j] * v.getElementAt(j);
        }
        values[i] = total;
    }

    return Vec4(values[0], values[1], values[2], values[3], v.colorId);
}

/*
* Returns a translation transformation matrix moving the point by x, y, z
*/
Matrix4 translationMatrix(double x, double y, double z)
{
    Matrix4 returnVal = getIdentityMatrix();
    returnVal.val[0][3] = x;
    returnVal.val[1][3] = y;
    returnVal.val[2][3] = z;
    return returnVal;
}

/*
* Returns a rotation transformation matrix rotating the point around x axis by angle
*/
Matrix4 rotationxMatrix(double angle)
{
    Matrix4 returnVal = getIdentityMatrix();
    double angCos = cos(angle * PI / 180.0);
    double angSin = sin(angle * PI / 180.0);
    returnVal.val[1][1] = angCos;
    returnVal.val[1][2] = -angSin;
    returnVal.val[2][1] = angSin;
    returnVal.val[2][2] = angCos;
    return returnVal;
}

/*
* Returns a rotation transformation matrix rotating the point around y axis by angle
*/
Matrix4 rotationyMatrix(double angle)
{
    Matrix4 returnVal = getIdentityMatrix();
    double angCos = cos(angle * PI / 180.0);
    double angSin = sin(angle * PI / 180.0);
    returnVal.val[0][0] = angCos;
    returnVal.val[0][2] = angSin;
    returnVal.val[2][0] = -angSin;
    returnVal.val[2][2] = angCos;
    return returnVal;
}

/*
* Returns a rotation transformation matrix rotating the point around x axis by angle
*/
Matrix4 rotationzMatrix(double angle)
{
    Matrix4 returnVal = getIdentityMatrix();
    double angCos = cos(angle * PI / 180.0);
    double angSin = sin(angle * PI / 180.0);
    returnVal.val[0][0] = angCos;
    returnVal.val[0][1] = -angSin;
    returnVal.val[1][0] = angSin;
    returnVal.val[1][1] = angCos;
    return returnVal;
}

/*
* Returns a rotation transformation matrix rotating the point around (0, 0, 0), (x, y, z) axis by angle
*/
Matrix4 rotationMatrix(double angle, double x, double y, double z)
{
    Vec3 u = Vec3(x, y, z, 0);
    Vec3 v;
    Vec3 w;
    Matrix4 M;
    Matrix4 Mrev;
    u = normalizeVec3(u);
    double ux = ABS(u.x);
    double uy = ABS(u.y);
    double uz = ABS(u.z);

    if (ux <= uy && ux <= uz)
    {
        v.x = 0;
        v.y = -u.z;
        v.z = u.y;
    }
    else if (uy <= ux && uy <= uz)
    {
        v.x = -u.z;
        v.y = 0;
        v.z = u.x;
    }
    else
    {
        v.x = -u.y;
        v.y = u.x;
        v.z = 0;
    }

    w = crossProductVec3(u, v);
    v = normalizeVec3(v);
    w = normalizeVec3(w);

    M.val[0][0] = u.x;
    M.val[0][1] = u.y;
    M.val[0][2] = u.z;
    M.val[1][0] = v.x;
    M.val[1][1] = v.y;
    M.val[1][2] = v.z;
    M.val[2][0] = w.x;
    M.val[2][1] = w.y;
    M.val[2][2] = w.z;
    M.val[3][3] = 1.0;

    Mrev.val[0][0] = u.x;
    Mrev.val[1][0] = u.y;
    Mrev.val[2][0] = u.z;
    Mrev.val[0][1] = v.x;
    Mrev.val[1][1] = v.y;
    Mrev.val[2][1] = v.z;
    Mrev.val[0][2] = w.x;
    Mrev.val[1][2] = w.y;
    Mrev.val[2][2] = w.z;
    Mrev.val[3][3] = 1.0;

    return multiplyMatrixWithMatrix(Mrev, multiplyMatrixWithMatrix(rotationxMatrix(angle), M));
}

/*
* Returns a scaling transformation matrix scaling the point by x, y, z
*/
Matrix4 scalingMatrix(double x, double y, double z)
{
    Matrix4 returnVal;
    returnVal.val[0][0] = x;
    returnVal.val[1][1] = y;
    returnVal.val[2][2] = z;
    returnVal.val[3][3] = 1.0;
    return returnVal;
}

/*
* Return camera transformation matrix which sets camera pos as (0, 0, 0) and u, v, w ad major axises
*/
Matrix4 cameraTransformMatrix(Vec3 pos, Vec3 u, Vec3 v, Vec3 w)
{
    Matrix4 M;
    M.val[0][0] = u.x;
    M.val[0][1] = u.y;
    M.val[0][2] = u.z;
    M.val[1][0] = v.x;
    M.val[1][1] = v.y;
    M.val[1][2] = v.z;
    M.val[2][0] = w.x;
    M.val[2][1] = w.y;
    M.val[2][2] = w.z;
    M.val[3][3] = 1.0;

    return multiplyMatrixWithMatrix(M, translationMatrix(-pos.x, -pos.y, -pos.z));
}

/*
*  Retrurn orthographic projection matrix with given near and dar planes.
*/
Matrix4 orthographicProjectionMatrix(double left, double right, double bottom, double top, double near, double far)
{
    Matrix4 returnVal;
    returnVal.val[0][0] = 2.0 / (right - left);
    returnVal.val[0][3] = -(right + left) / (right - left);
    returnVal.val[1][1] = 2.0 / (top - bottom);
    returnVal.val[1][3] = -(top + bottom) / (top - bottom);
    returnVal.val[2][2] = -2.0 / (far - near);
    returnVal.val[2][3] = -(far + near) / (far - near);
    returnVal.val[3][3] = 1.0;
    return returnVal;
}

/*
*  Return p2o matriz which turns a ortographic projection matrix into perspective projection matrix.
*/
Matrix4 p2OProjectionMatrix(double far, double near)
{
    Matrix4 returnVal;
    returnVal.val[0][0] = near;
    returnVal.val[1][1] = near;
    returnVal.val[2][2] = far + near;
    returnVal.val[2][3] = far * near;
    returnVal.val[3][2] = -1;
    return returnVal;
}

/*
*  Retrurn perspective projection matrix with given near and dar planes.
*/
Matrix4 perspectiveProjectionMatrix(double left, double right, double bottom, double top, double near, double far)
{
    return multiplyMatrixWithMatrix(orthographicProjectionMatrix(left, right, bottom, top, near, far), p2OProjectionMatrix(far, near));
}

Matrix4 viewPortMatrix(int x, int y)
{
    double nx = (double)x;
    double ny = (double)y;
    Matrix4 returnVal;
    returnVal.val[0][0] = nx / 2.0;
    returnVal.val[0][3] = (nx - 1.0) / 2.0;
    returnVal.val[1][1] = ny / 2.0;
    returnVal.val[1][3] = (ny - 1.0) / 2.0;
    returnVal.val[2][2] = 0.5;
    returnVal.val[2][3] = 0.5;
    returnVal.val[3][3] = 1.0;
    
    return returnVal;
}

Vec3 perspectiveDivide(Vec4 vec)
{
    Vec3 returnVal;
    double w = vec.t;
    returnVal.x = vec.x / w;
    returnVal.y = vec.y / w;
    returnVal.z = vec.z / w;
    returnVal.colorId = vec.colorId;
    return returnVal;
}

Vec3 transformVertice(Vec3 vec, Matrix4 matrix)
{
    Vec4 returnVal = multiplyMatrixWithVec4(matrix, Vec4(vec.x, vec.y, vec.z, 1, vec.colorId));
    return perspectiveDivide(returnVal);
}

Vec3 getNormal(Vec3 a, Vec3 b, Vec3 c)
{
    return crossProductVec3(subtractVec3(c, b), subtractVec3(a, b));
}

Vec3 barycentricToCartesian(Vec3 ver1, Vec3 ver2, Vec3 ver3, double alpha, double beta)
{
    return addVec3(multiplyVec3WithScalar(ver1, alpha), addVec3(multiplyVec3WithScalar(ver2, beta), multiplyVec3WithScalar(ver3, 1 - alpha - beta)));
}
