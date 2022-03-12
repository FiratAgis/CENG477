#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>

#include "Scene.h"
#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "tinyxml2.h"
#include "Helpers.h"
#include <chrono>


using namespace tinyxml2;
using namespace std;


void Scene::initilizeModelingTransformations()
{
	int transSize = this->translations.size();
	int rotSize = this->rotations.size();
	int scaleSize = this->scalings.size();
	int meshSize = this->meshes.size();
	Rotation* rot;
	Translation* trans;
	Scaling* scale;

	for (int i = 0; i < transSize; i++)
	{
		trans = translations[i];
		translationMatrixes.insert(translationMatrixes.begin() + i, translationMatrix(trans->tx, trans->ty, trans->tz));
	}

	for (int i = 0; i < rotSize; i++)
	{
		rot = rotations[i];
		rotationMatrixes.insert(rotationMatrixes.begin() + i, rotationMatrix(rot->angle, rot->ux, rot->uy, rot->uz));
	}

	for (int i = 0; i < scaleSize; i++)
	{
		scale = scalings[i];
		scalingMatrixes.insert(scalingMatrixes.begin() + i, scalingMatrix(scale->sx, scale->sy, scale->sz));
	}

	for (int i = 0; i < meshSize; i++)
	{
		Mesh* mesh = meshes[i];
		Matrix4 meshMatrix = getIdentityMatrix();
		for (int j = 0; j < mesh->numberOfTransformations; j++)
		{
			if (mesh->transformationTypes[j] == 't')
			{
				meshMatrix = multiplyMatrixWithMatrix(translationMatrixes[mesh->transformationIds[j] - 1], meshMatrix);
			}
			else if (mesh->transformationTypes[j] == 's')
			{
				meshMatrix = multiplyMatrixWithMatrix(scalingMatrixes[mesh->transformationIds[j] - 1], meshMatrix);
			}
			else
			{
				meshMatrix = multiplyMatrixWithMatrix(rotationMatrixes[mesh->transformationIds[j] - 1], meshMatrix);
			}
		}
		meshModelingMatrix.insert(meshModelingMatrix.begin() + i, meshMatrix);
	}
}

void Scene::initilizaCameraTransformations()
{
	int camSize = cameras.size();
	for (int i = 0; i < camSize; i++)
	{
		Camera* cam = cameras[i];
		Matrix4 cameraTransformation = cameraTransformMatrix(cam->pos, cam->u, cam->v, cam->w);
		if (cam->projectionType == 0)
		{
			cameraMatrixes.insert(cameraMatrixes.begin() + i, multiplyMatrixWithMatrix(orthographicProjectionMatrix(cam->left, cam->right, cam->bottom, cam->top, cam->near, cam->far), cameraTransformation));
		}
		else
		{
			cameraMatrixes.insert(cameraMatrixes.begin() + i, multiplyMatrixWithMatrix(perspectiveProjectionMatrix(cam->left, cam->right, cam->bottom, cam->top, cam->near, cam->far), cameraTransformation));
		}
	}
}

void Scene::transformMesh(Mesh* mesh, Matrix4 transformation)
{
	int verticeAmount = vertices.size();
	for (int i = 0; i < verticeAmount; i++)
	{
		mesh->transformedVertices.insert(mesh->transformedVertices.begin() + i, transformVertice(*(vertices[i]), transformation));
	}
}

void Scene::transformMeshes(Camera* camera)
{
	int meshCount = meshes.size();
	Matrix4 cameraMatrix = cameraMatrixes[camera->cameraId - 1];
	for (int i = 0; i < meshCount; i++)
	{
		transformMesh(meshes[i], multiplyMatrixWithMatrix(cameraMatrix, meshModelingMatrix[i]));
	}
}

void Scene::clearMeshes()
{
	for (Mesh* mesh : meshes)
	{
		mesh->transformedVertices.clear();
	}
}

void Scene::viewPort(Matrix4 viewPortTransformationMatrix, Mesh* mesh)
{
	int verticeAmount = mesh->transformedVertices.size();
	for (int i = 0; i < verticeAmount; i++)
	{
		Vec3 temp2 = mesh->transformedVertices[i];
		Vec4 temp = multiplyMatrixWithVec4(viewPortTransformationMatrix, Vec4(temp2.x, temp2.y, temp2.z, 1, temp2.colorId));
		mesh->transformedVertices[i] = Vec3(temp.x, temp.y, temp.z, temp.colorId);

	}
}

double Scene::f(Vec3 ver1, Vec3 ver2, double i, double j)
{
	return i * (ver1.y - ver2.y) + j * (ver2.x - ver1.x) + ver1.x * ver2.y - ver1.y * ver2.x;
}

void Scene::renderLine(Vec3 ver1, Vec3 ver2, int camMaxx, int camMaxy)
{
	if(ABS(ver1.x - ver2.x) > ABS(ver1.y - ver2.y))
	{
		double x0 = min(ver1.x, ver2.x);
		double x1 = max(ver1.x, ver2.x);
		double y0 = ABS(x0 - ver1.x) < EPSILON ? ver1.y: ver2.y;
		double y1 = ABS(x1 - ver1.x) < EPSILON ? ver1.y: ver2.y;

		double xdif = x1 - x0;
		double ydif = y1 - y0;
		for(int x = (int)max(0.0, x0); x <= min((double)camMaxx - 1, x1); x++)
		{
			double alpha = (x - x0) / xdif;
			int y = (int)(alpha  * ydif) + y0;
			if (y >= 0.0 && y < camMaxy)
			{
				Color* color1 = colorsOfVertices[ABS(x0 - ver1.x) < EPSILON ? ver1.colorId - 1 : ver2.colorId - 1];
				Color* color2 = colorsOfVertices[ABS(x1 - ver1.x) < EPSILON ? ver1.colorId - 1 : ver2.colorId - 1];

				image[x][y].r = makeBetweenZeroAnd255((color1->r * (1 - alpha)) + (alpha * color2->r));
				image[x][y].g = makeBetweenZeroAnd255((color1->g * (1 - alpha)) + (alpha * color2->g));
				image[x][y].b = makeBetweenZeroAnd255((color1->b * (1 - alpha)) + (alpha * color2->b));
			}

		}
	}
	else
	{
		double y0 = min(ver1.y, ver2.y);
		double y1 = max(ver1.y, ver2.y);
		double x0 = ABS(y0 - ver1.y) < EPSILON ? ver1.x: ver2.x;
		double x1 = ABS(y1 - ver1.y) < EPSILON ? ver1.x: ver2.x;

		double xdif = x1 - x0;
		double ydif = y1 - y0;
		for(int y = (int)max(0.0, y0); y <= min((double)camMaxy - 1, y1); y++)
		{
			double alpha = (y - y0) / ydif;
			int x = (int)(alpha  * xdif) + x0;
			if (x >= 0.0 && x < camMaxx)
			{
				Color* color1 = colorsOfVertices[ABS(y0 - ver1.y) < EPSILON ? ver1.colorId - 1 : ver2.colorId -1];
				Color* color2 = colorsOfVertices[ABS(y1 - ver1.y) < EPSILON ? ver1.colorId - 1 : ver2.colorId - 1];

				image[x][y].r = makeBetweenZeroAnd255((color1->r * (1 - alpha)) + (alpha * color2->r));
				image[x][y].g = makeBetweenZeroAnd255((color1->g * (1 - alpha)) + (alpha * color2->g));
				image[x][y].b = makeBetweenZeroAnd255((color1->b * (1 - alpha)) + (alpha * color2->b));
			}
		}
	}
}

void Scene::renderWireframeTriangle(Vec3 ver1, Vec3 ver2, Vec3 ver3, int camMaxx, int camMaxy)
{
	renderLine(ver1, ver2, camMaxx, camMaxy);
	renderLine(ver1, ver3, camMaxx, camMaxy);
	renderLine(ver2, ver3, camMaxx, camMaxy);
}

void Scene::renderTriangle(Vec3 ver1, Vec3 ver2, Vec3 ver3, int camMaxx, int camMaxy)
{
	int xmax, ymax, xmin, ymin;
	xmax = min((double)camMaxx - 1, max(ver1.x, max(ver2.x, ver3.x)));
	ymax = min((double)camMaxy - 1, max(ver1.y, max(ver2.y, ver3.y)));
	xmin = max(0.0, min(ver1.x, min(ver2.x, ver3.x)));
	ymin = max(0.0, min(ver1.y, min(ver2.y, ver3.y)));
	for (int i = xmin; i <= xmax; i++)
	{
		for (int j = ymin; j <= ymax; j++)
		{
			double alpha, beta, gama, f01, f12, f20;
			alpha = f(ver2, ver3, i, j) / f(ver2, ver3, ver1.x, ver1.y);
			beta = f(ver3, ver1, i, j) / f(ver3, ver1, ver2.x, ver2.y);
			gama = f(ver1, ver2, i, j) / f(ver1, ver2, ver3.x, ver3.y);
			if (alpha >= 0 && beta >= 0 && gama >= 0)
			{
				Color* color1 = colorsOfVertices[ver1.colorId - 1];
				Color* color2 = colorsOfVertices[ver2.colorId - 1];
				Color* color3 = colorsOfVertices[ver3.colorId - 1];

				
				this->image[i][j].r = makeBetweenZeroAnd255(color1->r * alpha + color2->r * beta + color3->r * gama);
				this->image[i][j].g = makeBetweenZeroAnd255(color1->g * alpha + color2->g * beta + color3->g * gama);
				this->image[i][j].b = makeBetweenZeroAnd255(color1->b * alpha + color2->b * beta + color3->b * gama);
			}
		}
	}
}

bool Scene::backCull(Vec3 ver1, Vec3 ver2, Vec3 ver3) 
{
	Vec3 normal = getNormal(ver1, ver2, ver3);

	if (dotProductVec3(normal, barycentricToCartesian(ver1, ver2, ver3, 1.0 / 3.0, 1.0 / 3.0)) > 0.0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//Slayttaki clip kodu, bunu da aynen aldım. Return değeri olarak Vertice1 ve Vertice2 çeviriyor, clipplendi veya cliplenmedi farketmeden.
vector<Vec3> Scene::clipper(Vec3 first, Vec3 second, bool & changed_1, bool & changed_2) {
		Vec3 minVertex = Vec3(-1.0,-1.0,-1.0,-1);
		Vec3 maxVertex = Vec3(1.0,1.0,1.0,-1);
		vector<Vec3> ret_vec ;
		Vec3 clipped_1, clipped_2 ;
		//Cliplenmemişse de hala vertex döndürmek gerektiği için önce gelen değerlere eşitliyorum
		clipped_1.x = first.x ;
		clipped_1.y = first.y ;
		clipped_1.z = first.z ;
		clipped_2.x = second.x ;
		clipped_2.y = second.y ;
		clipped_2.z = second.z ;
		double te = 0;
		double tl = 1;
		double dx = second.x - first.x ;
		double dy = second.y - first.y ;
		double dz = second.z - first.z ;
		//Buradan problem çıkmaz sanırım, ama çıkarsa iç içe devasa ifception a da çevrilebilir,
		//slaytta öyle duruyodu çünkü.
		if ( isVisible(dx, minVertex.x - first.x, te, tl ) &&
			 isVisible(-dx,first.x - maxVertex.x, te, tl ) &&
			 isVisible(dy, minVertex.y - first.y, te, tl ) &&
			 isVisible(-dy,first.y - maxVertex.y, te, tl ) &&
			 isVisible(dz, minVertex.z - first.z, te, tl ) &&
			 isVisible(-dz,first.z - maxVertex.z, te, tl )  )
		{
			if (tl < 1) {
				clipped_1.x = first.x + dx * tl ;
				clipped_1.y = first.y + dy * tl ;
				clipped_1.z = first.z + dz * tl ;
				changed_1 = true;
			}
			if (te > 0) {
				clipped_2.x = first.x + dx * te ;
				clipped_2.y = first.y + dy * te ;
				clipped_2.z = first.z + dy * te ;
				changed_2 = true;
			}
		}
		ret_vec.push_back(clipped_1) ;
		ret_vec.push_back(clipped_2) ;
		return ret_vec ;
}

void Scene::clipTriangles(Mesh* mesh) {
	//Hangi vertexler bölündü kontrol etmek için boolean values
	bool changed_1, changed_2, changed_3 ;
	//bölünürse diye vertex id'leri tutmak için integerlar. v_11 : 1. vertice 1. clip, v_12 : 1. vertice 2. clip ...
	int v_11, v_12, v_21, v_22, v_31, v_32 ;
	for (Triangle triangle : mesh->triangles) {
		vector<int> new_vertice_ids ; //yeni oluşan köşeler için vector
		vector<Vec3> coords ; // cliplendikten sonra koordinatları tutmak için vector
		Vec3 ver1 = mesh->getVertex(triangle.getFirstVertexId()); //vertice 1
		Vec3 ver2 = mesh->getVertex(triangle.getSecondVertexId()); // vertice 2
		Vec3 ver3 = mesh->getVertex(triangle.getThirdVertexId()); // vertice 3

		v_11 = triangle.getFirstVertexId() ; //değişmezse default id'si önceki lokasyonu köşenin
		v_12 = v_11 ;
		v_21 = triangle.getSecondVertexId() ;
		v_22 = v_21 ;
		v_31 = triangle.getThirdVertexId() ;
		v_32 = triangle.getThirdVertexId() ;
		//hiçbiri değişmedi'ye set
		changed_1 = false ;
		changed_2 = false ;
		changed_3 = false ;

		//1,2 - 1,3 - 2,3 köşelerini cliplemeye yollayıp, sonra geri gelen iki köşeyi vector'e pushla.
		//clip fonksiyonu hangisi cliplendi ise changed' değerini True'ya setliyor.
		vector<Vec3> clip1, clip2, clip3;
		clip1 = clipper(ver1, ver2, changed_1, changed_2);
		clip2 = clipper(ver1, ver2, changed_1, changed_3);
		clip3 = clipper(ver1, ver2, changed_2, changed_3); 
		for (int i = 0; i < 2; i++ ) {
		coords.push_back( clip1[i] );
		coords.push_back( clip2[i] );
		coords.push_back( clip3[i] );
		}


		//eğer 1. değişti ise, coords'daki 0 ve 2. konumdaki elemanlar cliplenmiş versiyonları demektir.
		//onları transformed Vertice'e pushlayıp, konumlarını v_11 ve v_12'ye tutuyorum, onu da new_vertice_ids'e atıyorum
		//sonrasında üçgen oluşturmak için.
		if (changed_1)
		{
			mesh->transformedVertices.push_back(coords[0]) ;
			v_11 = mesh->transformedVertices.size() ;
			mesh->transformedVertices.push_back(coords[2]) ;
			v_12 = mesh->transformedVertices.size() ;
			new_vertice_ids.push_back(v_11) ;
			new_vertice_ids.push_back(v_12);
		}
		//bir şey değişmemişse sadece v_11'i atıyorum, eski id o çünkü.
		else
		{
			new_vertice_ids.push_back(v_11) ;
		}
		//Same goes here
		if (changed_2)
		{
			mesh->transformedVertices.push_back(coords[1]) ;
			v_21 = mesh->transformedVertices.size() ;
			mesh->transformedVertices.push_back(coords[4]) ;
			v_22 = mesh->transformedVertices.size() ;
			new_vertice_ids.push_back(v_21) ;
			new_vertice_ids.push_back(v_22) ;

		}
		else
		{
			new_vertice_ids.push_back(v_21) ;

		}
		//Same goes here
		if (changed_3)
		{
			mesh->transformedVertices.push_back(coords[3]) ;
			v_31 = mesh->transformedVertices.size() ;
			mesh->transformedVertices.push_back(coords[5]) ;
			v_32 = mesh->transformedVertices.size() ;
			new_vertice_ids.push_back(v_31) ;
			new_vertice_ids.push_back(v_32) ;

		}
		else
		{
			new_vertice_ids.push_back(v_31) ;
		}
		//new_vertice_ids içinde şimdi elimizdeki bütün yeni köşeler var.
		//İlk üçü ile üçgen oluşturuyorum. Sonra ikinci köşeyi listeden siliyorum ki overlapping üçgenler oluşturmayalım.
		//Bunu iki köşe kalıncaya kadar yapıyorum. Böylece atıyorum
		// [A,B,C,D,E] köşeleri varsa, ABC - ACD - ADE üçgenleri oluşuyor.
		// Yalnız, Triangle ile oluşturup atıyorum ama, umarım silinmiyordur sonrasında. (new ile allocate etmem gerekli mi bilemedim.)
		while ( new_vertice_ids.size() > 2) {
			mesh->transformedTriangles.push_back( Triangle(new_vertice_ids[0], new_vertice_ids[1], new_vertice_ids[2] )) ;
			new_vertice_ids.erase(new_vertice_ids.begin() + 1) ;
		}
	}
}

//Slayytaki isVisible kodu, te ve tl'yi editleyerek gidiyordu o da öyle yaptım.
bool Scene::isVisible(double den, double num, double &te, double &tl) {
	int t ;
	if (den > 0) {
		t = num / den ;
		if (t > tl) {
			return false;
		}
		if (t > te) {
			te = t ;
		}
	}
	else if ( den < 0) {
		t = num / den ;
		if ( t > te) {
			return false ;
		}
		if (t < tl) {
			tl = t ;
		}
	}
	else if (num > 0) {
		return false ;
	}
	return true;
}

bool Scene::volumeCulling(Mesh* mesh) {
	//min'i en yukarıdan, max'i en aşağıdan
	Vec3 minVertex = Vec3(1.0,1.0,1.0,-1);
	Vec3 maxVertex = Vec3(-1.0,-1.0,-1.0,-1);
	for (Triangle triangle : mesh->triangles) {
		Vec3 ver1 = mesh->getVertex(triangle.getFirstVertexId());
		Vec3 ver2 = mesh->getVertex(triangle.getSecondVertexId());
		Vec3 ver3 = mesh->getVertex(triangle.getThirdVertexId());
		minVertex.x = fmin( fmin( fmin(minVertex.x, ver1.x), ver2.x )  , ver3.x ) ;
		minVertex.y = fmin( fmin( fmin(minVertex.y, ver1.y), ver2.y )  , ver3.y ) ;
		minVertex.z = fmin( fmin( fmin(minVertex.z, ver1.z), ver2.z )  , ver3.z ) ;

		maxVertex.x = fmax( fmax( fmax(maxVertex.x, ver1.x), ver2.x )  , ver3.x ) ;
		maxVertex.x = fmax( fmax( fmax(maxVertex.y, ver1.y), ver2.y )  , ver3.y ) ;
		maxVertex.x = fmax( fmax( fmax(maxVertex.z, ver1.z), ver2.z )  , ver3.z ) ;

	}
	if (minVertex.x > 1.0 || minVertex.y > 1.0 || minVertex.z > 1.0 ||
	    maxVertex.x < -1.0 || maxVertex.y < -1.00 || maxVertex.z < -1.00 ) {
			return false;
	}
	else
	{
		return true;
	}
}

bool Scene::triangleCulling(Triangle triangle, Mesh* mesh) {
	//min'i en yukarıdan, max'i en aşağıdan
	Vec3 minVertex = Vec3(1.0,1.0,1.0,-1);
	Vec3 maxVertex = Vec3(-1.0,-1.0,-1.0,-1);

	Vec3 ver1 = mesh->getVertex(triangle.getFirstVertexId());
	Vec3 ver2 = mesh->getVertex(triangle.getSecondVertexId());
	Vec3 ver3 = mesh->getVertex(triangle.getThirdVertexId());

	minVertex.x = fmin( fmin( fmin(minVertex.x, ver1.x), ver2.x )  , ver3.x ) ;
	minVertex.y = fmin( fmin( fmin(minVertex.y, ver1.y), ver2.y )  , ver3.y ) ;
	minVertex.z = fmin( fmin( fmin(minVertex.z, ver1.z), ver2.z )  , ver3.z ) ;

	maxVertex.x = fmax( fmax( fmax(maxVertex.x, ver1.x), ver2.x )  , ver3.x ) ;
	maxVertex.x = fmax( fmax( fmax(maxVertex.y, ver1.y), ver2.y )  , ver3.y ) ;
	maxVertex.x = fmax( fmax( fmax(maxVertex.z, ver1.z), ver2.z )  , ver3.z ) ;

	if (minVertex.x > 1.0 || minVertex.y || 1.0 or minVertex.z > 1.0 ||
	    maxVertex.x < -1.0 || maxVertex.y < -1.00 || maxVertex.z < -1.00 ) {
			return false;
	}
	else {
		return true;
	}
}

void Scene::forwardRenderingPipeline(Camera* camera)
{
	// TODO: Implement this function.
	transformMeshes(camera);

	int meshCount = meshes.size();
	Matrix4 view = viewPortMatrix(camera->horRes, camera->verRes);
	for (int i = 0; i < meshCount; i++)
	{
		Mesh* mesh = meshes[i];
		vector <Triangle> willBeRendered;
		for (Triangle triangle : mesh->triangles)
		{
			Vec3 ver1 = mesh->getVertex(triangle.getFirstVertexId());
			Vec3 ver2 = mesh->getVertex(triangle.getSecondVertexId());
			Vec3 ver3 = mesh->getVertex(triangle.getThirdVertexId());
			if (!cullingEnabled || (cullingEnabled && backCull(ver1, ver2, ver3)))
			{
				willBeRendered.push_back(triangle);
			}
		}

		viewPort(view, mesh);
		for (Triangle triangle : willBeRendered)
		{
			Vec3 ver1 = mesh->getVertex(triangle.getFirstVertexId());
			Vec3 ver2 = mesh->getVertex(triangle.getSecondVertexId());
			Vec3 ver3 = mesh->getVertex(triangle.getThirdVertexId());
			if (!cullingEnabled || (cullingEnabled && backCull(ver1, ver2, ver3)))
			{
				if (mesh->type == 0)
				{
					renderWireframeTriangle(ver1, ver2, ver3, camera->horRes, camera->verRes);
				}
				else
				{
					renderTriangle(ver1, ver2, ver3, camera->horRes, camera->verRes);
				}
			}
		}
		willBeRendered.clear();
	}
	clearMeshes();

}

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *pElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *pRoot = xmlDoc.FirstChild();

	// read background color
	pElement = pRoot->FirstChildElement("BackgroundColor");
	str = pElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	pElement = pRoot->FirstChildElement("Culling");
	if (pElement != NULL) {
		str = pElement->GetText();

		if (strcmp(str, "enabled") == 0) {
			cullingEnabled = true;
		}
		else {
			cullingEnabled = false;
		}
	}

	// read cameras
	pElement = pRoot->FirstChildElement("Cameras");
	XMLElement *pCamera = pElement->FirstChildElement("Camera");
	XMLElement *camElement;
	while (pCamera != NULL)
	{
		Camera *cam = new Camera();

		pCamera->QueryIntAttribute("id", &cam->cameraId);

		// read projection type
		str = pCamera->Attribute("type");

		if (strcmp(str, "orthographic") == 0) {
			cam->projectionType = 0;
		}
		else {
			cam->projectionType = 1;
		}

		camElement = pCamera->FirstChildElement("Position");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->pos.x, &cam->pos.y, &cam->pos.z);

		camElement = pCamera->FirstChildElement("Gaze");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->gaze.x, &cam->gaze.y, &cam->gaze.z);

		camElement = pCamera->FirstChildElement("Up");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->v.x, &cam->v.y, &cam->v.z);

		cam->gaze = normalizeVec3(cam->gaze);
		cam->u = crossProductVec3(cam->gaze, cam->v);
		cam->u = normalizeVec3(cam->u);

		cam->w = inverseVec3(cam->gaze);
		cam->v = crossProductVec3(cam->u, cam->gaze);
		cam->v = normalizeVec3(cam->v);

		camElement = pCamera->FirstChildElement("ImagePlane");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &cam->left, &cam->right, &cam->bottom, &cam->top,
			   &cam->near, &cam->far, &cam->horRes, &cam->verRes);

		camElement = pCamera->FirstChildElement("OutputName");
		str = camElement->GetText();
		cam->outputFileName = string(str);

		cameras.push_back(cam);

		pCamera = pCamera->NextSiblingElement("Camera");
	}

	// read vertices
	pElement = pRoot->FirstChildElement("Vertices");
	XMLElement *pVertex = pElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (pVertex != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = pVertex->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = pVertex->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		vertices.push_back(vertex);
		colorsOfVertices.push_back(color);

		pVertex = pVertex->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	pElement = pRoot->FirstChildElement("Translations");
	XMLElement *pTranslation = pElement->FirstChildElement("Translation");
	while (pTranslation != NULL)
	{
		Translation *translation = new Translation();

		pTranslation->QueryIntAttribute("id", &translation->translationId);

		str = pTranslation->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		translations.push_back(translation);

		pTranslation = pTranslation->NextSiblingElement("Translation");
	}

	// read scalings
	pElement = pRoot->FirstChildElement("Scalings");
	XMLElement *pScaling = pElement->FirstChildElement("Scaling");
	while (pScaling != NULL)
	{
		Scaling *scaling = new Scaling();

		pScaling->QueryIntAttribute("id", &scaling->scalingId);
		str = pScaling->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		scalings.push_back(scaling);

		pScaling = pScaling->NextSiblingElement("Scaling");
	}

	// read rotations
	pElement = pRoot->FirstChildElement("Rotations");
	XMLElement *pRotation = pElement->FirstChildElement("Rotation");
	while (pRotation != NULL)
	{
		Rotation *rotation = new Rotation();

		pRotation->QueryIntAttribute("id", &rotation->rotationId);
		str = pRotation->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		rotations.push_back(rotation);

		pRotation = pRotation->NextSiblingElement("Rotation");
	}

	// read meshes
	pElement = pRoot->FirstChildElement("Meshes");

	XMLElement *pMesh = pElement->FirstChildElement("Mesh");
	XMLElement *meshElement;
	while (pMesh != NULL)
	{
		Mesh *mesh = new Mesh();

		pMesh->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = pMesh->Attribute("type");

		if (strcmp(str, "wireframe") == 0) {
			mesh->type = 0;
		}
		else {
			mesh->type = 1;
		}

		// read mesh transformations
		XMLElement *pTransformations = pMesh->FirstChildElement("Transformations");
		XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

		while (pTransformation != NULL)
		{
			char transformationType;
			int transformationId;

			str = pTransformation->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			pTransformation = pTransformation->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *clone_str;
		int v1, v2, v3;
		XMLElement *pFaces = pMesh->FirstChildElement("Faces");
        str = pFaces->GetText();
		clone_str = strdup(str);

		row = strtok(clone_str, "\n");
		while (row != NULL)
		{
			sscanf(row, "%d %d %d", &v1, &v2, &v3);
			mesh->triangles.push_back(Triangle(v1, v2, v3));
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		meshes.push_back(mesh);

		pMesh = pMesh->NextSiblingElement("Mesh");
	}

	//Added initilization
	initilizeModelingTransformations();
	initilizaCameraTransformations();
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFileName.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFileName << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}
