#ifndef _SCENE_H_
#define _SCENE_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "Vec4.h"
#include "Matrix4.h"

using namespace std;

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	vector< Matrix4 > translationMatrixes;
	vector< Matrix4 > rotationMatrixes;
	vector< Matrix4 > scalingMatrixes;
	vector< Matrix4 > meshModelingMatrix;
	vector< Matrix4 > cameraMatrixes;

	vector< vector<Color> > image;
	vector< Camera* > cameras;
	vector< Vec3* > vertices;
	vector< Color* > colorsOfVertices;
	vector< Scaling* > scalings;
	vector< Rotation* > rotations;
	vector< Translation* > translations;
	vector< Mesh* > meshes;

	Scene(const char *xmlPath);

	void initilizeModelingTransformations();
	void initilizaCameraTransformations();
	void transformMesh(Mesh* mesh, Matrix4 transformation);
	void transformMeshes(Camera* camera);
	double f(Vec3 ver1, Vec3 ver2, double x, double y);
	void viewPort(Matrix4 viewPortTransformationMatrix, Mesh* mesh);
	void renderTriangle(Vec3 ver1, Vec3 ver2, Vec3 ver3, int camMaxx, int camMaxy);
	bool backCull(Vec3 ver1, Vec3 ver2, Vec3 ver3);
	void clearMeshes();
	void renderLine(Vec3 ver1, Vec3 ver2, int camMaxx, int camMaxy);
	void renderWireframeTriangle(Vec3 ver1, Vec3 ver2, Vec3 ver3, int camMaxx, int camMaxy);

	void clipTriangles(Mesh* mesh);
	bool triangleCulling(Triangle triangle, Mesh* mesh);
	bool volumeCulling(Mesh* mesh);
	bool isVisible(double den, double num, double &te, double &tl);
	vector<Vec3> clipper(Vec3 first, Vec3 second, bool & changed_1, bool & changed_2);

	void initializeImage(Camera* camera);
	void forwardRenderingPipeline(Camera* camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera* camera);
	void convertPPMToPNG(string ppmFileName, int osType);
};

#endif
