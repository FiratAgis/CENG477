#ifndef __MESH_H__
#define __MESH_H__

#include <vector>
#include "Triangle.h"
#include "Vec4.h"
#include <iostream>

using namespace std;

class Mesh
{

public:
    int meshId;
    int type; // 0 for wireframe, 1 for solid
    int numberOfTransformations;
    vector< Vec3 > transformedVertices;
    vector<int> transformationIds;
    vector<char> transformationTypes;
    int numberOfTriangles;
    vector<Triangle> triangles;
    vector<Triangle> transformedTriangles;

    Mesh();
    Mesh(int meshId, int type, int numberOfTransformations,
          vector<int> transformationIds,
          vector<char> transformationTypes,
          int numberOfTriangles,
          vector<Triangle> triangles);

    Vec3 getVertex(int x);

    friend ostream &operator<<(ostream &os, const Mesh &m);
};

#endif
