#include <iostream>
#include <cmath>
#include <limits>
#include "parser.h"
#include "ppm.h"

typedef unsigned char RGB[3];


struct CameraInstance
{
    float pixelWidth;
    float pixelHeight;
    parser::Camera camera;
    parser::Vec3f u;
    parser::Vec3f v;
    parser::Vec3f w;
    parser::Vec3f m;
    parser::Vec3f q;
};

struct Rayf
{
    parser::Vec3f origin;
    parser::Vec3f direction;
};

struct TriangleCamera
{
    parser::Vec3f vec1;
    parser::Vec3f vec2;
    parser::Vec3f vec3;
    float det;
};

struct Surface
{
    parser::Triangle triangle;
    parser::Vec3f normal;
    std::vector<TriangleCamera> relationship;
};

parser::Scene scene;
std::vector<Surface> surfaces;
std::vector<CameraInstance> cameras;
enum close {none, surf, sph};

parser::Vec3f addVector(parser::Vec3f vec1, parser::Vec3f vec2)
{
    parser::Vec3f result;
    result.x = vec1.x + vec2.x;
    result.y = vec1.y + vec2.y;
    result.z = vec1.z + vec2.z;
    return result;
}

parser::Vec3f scaleVector(float scale, parser::Vec3f vec)
{
    parser::Vec3f result;
    result.x = scale * vec.x;
    result.y = scale * vec.y;
    result.z = scale * vec.z;
    return result;
}

float dotProduct(parser::Vec3f vec1, parser::Vec3f vec2)
{
    return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
}

parser::Vec3f crossProduct(parser::Vec3f vec1, parser::Vec3f vec2)
{
    parser::Vec3f result;
    result.x = (vec1.y * vec2.z) - (vec1.z * vec2.y);
    result.y = (vec1.z * vec2.x) - (vec1.x * vec2.z);
    result.z = (vec1.x * vec2.y) - (vec1.y * vec2.x);
    return result;
}

float vectorLength(parser::Vec3f vec)
{
    return sqrtf(powf(vec.x, 2) + powf(vec.y, 2) + powf(vec.z, 2));
}

float pointDistance(parser::Vec3f vec1, parser::Vec3f vec2)
{
    return vectorLength(addVector(vec1, scaleVector(-1, vec2)));
}

parser::Vec3f normalize(parser::Vec3f vec)
{
    parser::Vec3f result;
    float length = vectorLength(vec);
    result = scaleVector(1.0f / length, vec);
    return result;
}

parser::Vec3i castToInt(parser::Vec3f vec)
{
    parser::Vec3i result;
    result.x = (int)(vec.x + 0.5f);
    result.y = (int)(vec.y + 0.5f);
    result.z = (int)(vec.z + 0.5f);
    return result;
}

int clamp(int min, int val, int max)
{
    return std::max(min, std::min(val, max));
}

parser::Vec3i colorNormalize(parser::Vec3i vec)
{
    parser::Vec3i result;
    result.x = clamp(0, vec.x, 255);
    result.y = clamp(0, vec.y, 255);
    result.z = clamp(0, vec.z, 255);
    return result;
}

parser::Vec3i colorNormalize(parser::Vec3f vec)
{
    return colorNormalize(castToInt(vec));
}

parser::Vec3f combineVector(parser::Vec3f vec1, parser::Vec3f vec2)
{
    parser::Vec3f result;
    result.x = vec1.x * vec2.x;
    result.y = vec1.y * vec2.y;
    result.z = vec1.z * vec2.z;
    return result;
}

parser::Vec3f combineVector(parser::Vec3i vec1, parser::Vec3f vec2)
{
    parser::Vec3f result;
    result.x = vec1.x * vec2.x;
    result.y = vec1.y * vec2.y;
    result.z = vec1.z * vec2.z;
    return result;
}

parser::Vec3f combineVector(parser::Vec3f vec1, parser::Vec3i vec2)
{
    return combineVector(vec2, vec1);
}

float determinant(parser::Vec3f vec1, parser::Vec3f vec2, parser::Vec3f vec3)
{
    return (vec1.x * ((vec2.y * vec3.z) - (vec3.y * vec2.z)))
        + (vec1.y * ((vec3.x * vec2.z) - (vec2.x * vec3.z)))
        + (vec1.z * ((vec2.x * vec3.y) - (vec2.y * vec3.x)));
}

parser::Vec3f getVertex(int vertID)
{
    return scene.vertex_data.at(vertID - 1);
}

parser::Material getMaterial(int matID){
    return scene.materials.at(matID - 1);
}

TriangleCamera getRelationship(parser::Camera camera, parser::Face face)
{
    TriangleCamera result;

    parser::Vec3f a = getVertex(face.v0_id);
    parser::Vec3f b = getVertex(face.v1_id);
    parser::Vec3f c = getVertex(face.v2_id);

    result.vec1 = addVector(a, scaleVector(-1, b));
    result.vec2 = addVector(a, scaleVector(-1, c));
    result.vec3 = addVector(a, scaleVector(-1, camera.position));

    result.det = determinant(result.vec1, result.vec2, result.vec3);

    return result;
}

parser::Vec3f getNormal(parser::Face face)
{
    return normalize(crossProduct(addVector(getVertex(face.v2_id), scaleVector(-1, getVertex(face.v1_id))), addVector(getVertex(face.v0_id), scaleVector(-1, getVertex(face.v1_id)))));
}

float triangleIntersection(parser::Face face, Rayf ray)
{
    parser::Vec3f a = getVertex(face.v0_id);
    parser::Vec3f b = getVertex(face.v1_id);
    parser::Vec3f c = getVertex(face.v2_id);

    parser::Vec3f vec1 = addVector(a, scaleVector(-1, b));
    parser::Vec3f vec2 = addVector(a, scaleVector(-1, c));
    parser::Vec3f vec3 = addVector(a, scaleVector(-1, ray.origin));

    float detA = determinant(vec1, vec2, ray.direction);
    float beta = determinant(vec3, vec2, ray.direction) / detA;
    float gamma = determinant(vec1, vec3, ray.direction) / detA;
    if(beta + gamma <= 1 && beta >= 0 && gamma >= 0)
    {
        return determinant(vec1, vec2, vec3) / detA;
    }
    else
    {
        return -1;
    }
}

float triangleIntersection(TriangleCamera relationship, Rayf ray)
{
    float detA = determinant(relationship.vec1, relationship.vec2, ray.direction);
    float beta = determinant(relationship.vec3, relationship.vec2, ray.direction) / detA;
    float gamma = determinant(relationship.vec1, relationship.vec3, ray.direction) / detA;
    if(beta + gamma <= 1 && beta >= 0 && gamma >= 0)
    {
        return relationship.det / detA;
    }
    else
    {
        return -1;
    }
}

float sphereIntersection(parser::Sphere sphere, Rayf ray)
{
    float A,B,C, delta, t, d;
    parser::Vec3f c = getVertex(sphere.center_vertex_id);
    parser::Vec3f o_c = addVector(ray.origin, scaleVector(-1, c));
    delta = powf(dotProduct(ray.direction, o_c),2) - dotProduct(ray.direction, ray.direction)*(dotProduct(o_c, o_c) - powf(sphere.radius, 2));
    if (delta < 0) {
        t = -1;
    }
    else
    {
        d = dotProduct(ray.direction, ray.direction);
        delta = sqrtf(delta);
        t = -dotProduct(ray.direction, o_c);
        t = std::min((t + delta)/d , (t - delta)/d );
    }
    return t;
}

Rayf getEyeRay(CameraInstance instance, int x, int y)
{
    Rayf result;
    result.origin = instance.camera.position;
    float xDistance = (((float)x) + 0.5f) * instance.pixelWidth;
    float yDistance = (((float)y) + 0.5f) * instance.pixelHeight;
    result.direction = addVector(addVector(instance.q, addVector(scaleVector(xDistance, instance.u), scaleVector(-yDistance, instance.v))), scaleVector(-1, result.origin));

    return result;
}

parser::Vec3f getRecievedIrrediance(parser::Vec3f origin, parser::PointLight light)
{
    return scaleVector((1.0f / powf(pointDistance(origin, light.position), 2)), light.intensity);
}

parser::Vec3f getAmbient(Surface surface)
{
    parser::Vec3f light = scene.ambient_light;
    parser::Vec3f material = getMaterial(surface.triangle.material_id).ambient;
    return combineVector(light, material);
}

parser::Vec3f getAmbient(parser::Sphere sphere)
{
    parser::Vec3f light = scene.ambient_light;
    parser::Vec3f material = getMaterial(sphere.material_id).ambient;
    return combineVector(light, material);
}

parser::Vec3f getPoint(Rayf ray, float t)
{
    return addVector(ray.origin, scaleVector(t, ray.direction));
}

parser::Vec3f getSpecularDiffused(parser::Vec3f x, parser::Vec3f n, parser::Vec3f w0dir, parser::Material material, int count)
{
    parser::Vec3f specular;
    parser::Vec3f diffused;
    specular.x = 0;
    specular.y = 0;
    specular.z = 0;
    diffused.x = 0;
    diffused.y = 0;
    diffused.z = 0;
    for (parser::PointLight light : scene.point_lights)
    {
        parser::Vec3f widir = normalize(addVector(light.position, scaleVector(-1, x)));
        Rayf lightRay;
        lightRay.origin = x;
        lightRay.direction = widir;
        float t;
        float minT = pointDistance(x, light.position);
        bool unblocked = true;
        for(Surface& surface : surfaces)
        {
            t = triangleIntersection(surface.triangle.indices, lightRay);
            if(t > scene.shadow_ray_epsilon && t < minT)
            {
                unblocked = false;
                break;
            }
        }
        if(unblocked)
        {
            for (parser::Sphere& sphere : scene.spheres)
            {
                t = sphereIntersection(sphere, lightRay);
                if (t > scene.shadow_ray_epsilon && t < minT)
                {
                    unblocked = false;
                    break;
                }
            }
        }
        if(unblocked)
        {
            parser::Vec3f h = normalize(addVector(w0dir, widir));
            float costheta = dotProduct(widir, n);
            if (costheta > 0)
            {
                float cosalpha = powf(std::max(0.0f, dotProduct(n, h)), material.phong_exponent);
                parser::Vec3f radiance = getRecievedIrrediance(x, light);
                specular = addVector(specular, scaleVector(cosalpha, radiance));
                diffused = addVector(diffused, scaleVector(costheta, radiance));
            }
        }

    }
    specular = combineVector(specular, material.specular);
    diffused = combineVector(diffused, material.diffuse);
    if (material.is_mirror && count <= scene.max_recursion_depth) {
        Rayf reflect;
        reflect.origin = x;
        reflect.direction = normalize(addVector(scaleVector(-1,w0dir), scaleVector(2*dotProduct(n, w0dir), n)));
        Surface closestSurface;
        parser::Sphere closestSphere;
        close closePoint = none;
        float t;
        float minT = std::numeric_limits<float>::max();
        for(Surface& surface : surfaces)
        {
            t = triangleIntersection(surface.triangle.indices, reflect);
            if(t > scene.shadow_ray_epsilon && t < minT)
            {
                minT = t;
                closestSurface = surface;
                closePoint = surf;
            }
        }
        for (parser::Sphere& sphere : scene.spheres)
        {
            t = sphereIntersection(sphere, reflect);
            if (t > scene.shadow_ray_epsilon && t < minT)
            {
                minT = t;
                closestSphere = sphere;
                closePoint = sph;
            }
        }

        parser::Vec3f colorf;
        switch(closePoint)
        {
            case none:
                colorf.x = 0;
                colorf.y = 0;
                colorf.z = 0;
                break;
            case surf:
                colorf = addVector(getAmbient(closestSurface), getSpecularDiffused(getPoint(reflect, minT), closestSurface.normal, normalize(scaleVector(-1, reflect.direction)), getMaterial(closestSurface.triangle.material_id), count+1));
                colorf = combineVector(colorf, material.mirror);
                break;
            case sph:
                colorf = addVector(getAmbient(closestSphere), getSpecularDiffused(getPoint(reflect, minT), normalize(addVector(getPoint(reflect, minT), scaleVector(-1, getVertex(closestSphere.center_vertex_id)))), scaleVector(-1, reflect.direction), getMaterial(closestSphere.material_id), count+1));
                colorf = combineVector(colorf, material.mirror);
                break;
        }
        specular = addVector(specular, colorf);
    }
    return addVector(specular, diffused);
}

void computeImage(int camera_id)
{
    int imageRef = 0;
    CameraInstance instance = cameras.at(camera_id);
    unsigned char* image = new unsigned char [instance.camera.image_height *  instance.camera.image_width * 3];
    for(int i = 0; i < instance.camera.image_height; i++)
    {
        for(int j = 0; j < instance.camera.image_width; j++)
        {
            Rayf eyeray = getEyeRay(instance, j, i);
            close closePoint = none;
            float t;
            float minT = std::numeric_limits<float>::max();
            Surface closestSurface;
            parser::Sphere closestSphere;
            for(Surface& surface : surfaces)
            {
                t = triangleIntersection(surface.relationship[camera_id], eyeray);
                if(t > 0 && t < minT)
                {
                    minT = t;
                    closestSurface = surface;
                    closePoint = surf;
                }
            }
            for (parser::Sphere& sphere : scene.spheres)
            {
                t = sphereIntersection(sphere, eyeray);
                if (t > 0 && t < minT)
                {
                    minT = t;
                    closestSphere = sphere;
                    closePoint = sph;
                }
            }
            parser::Vec3f colorf;
            parser::Vec3i colori;
            switch(closePoint)
            {
                case none:
                    image[imageRef++] = scene.background_color.x;
                    image[imageRef++] = scene.background_color.y;
                    image[imageRef++] = scene.background_color.z;
                    break;
                case surf:
                    colorf = addVector(getAmbient(closestSurface), getSpecularDiffused(getPoint(eyeray, minT), closestSurface.normal, normalize(scaleVector(-1, eyeray.direction)), getMaterial(closestSurface.triangle.material_id), 0));

                    colori = colorNormalize(colorf);
                    image[imageRef++] = colori.x;
                    image[imageRef++] = colori.y;
                    image[imageRef++] = colori.z;
                    break;
                case sph:
                    colorf = addVector(getAmbient(closestSphere), getSpecularDiffused(getPoint(eyeray, minT), normalize(addVector(getPoint(eyeray, minT), scaleVector(-1, getVertex(closestSphere.center_vertex_id)))), scaleVector(-1, eyeray.direction), getMaterial(closestSphere.material_id), 0));

                    colori = colorNormalize(colorf);
                    image[imageRef++] = colori.x;
                    image[imageRef++] = colori.y;
                    image[imageRef++] = colori.z;
                    break;
            }
        }
    }
    write_ppm(instance.camera.image_name.c_str(), image, instance.camera.image_width, instance.camera.image_height);
}

void setup()
{
    for (parser::Triangle triangle : scene.triangles)
    {
        Surface surface;
        surface.triangle = triangle;
        surface.normal = getNormal(triangle.indices);
        surfaces.push_back(surface);
    }

    for (parser::Mesh mesh : scene.meshes)
    {
        for(parser::Face face : mesh.faces)
        {
            Surface surface;
            surface.triangle.indices = face;
            surface.triangle.material_id = mesh.material_id;
            surface.normal = getNormal(surface.triangle.indices);
            surfaces.push_back(surface);
        }
    }

    for(parser::Camera camera : scene.cameras)
    {
        CameraInstance instance;
        instance.camera = camera;
        parser::Vec4f viewPlane = camera.near_plane;
        instance.pixelWidth = (viewPlane.y - viewPlane.x) / camera.image_width;
        instance.pixelHeight = (viewPlane.w - viewPlane.z) / camera.image_height;
        instance.w = scaleVector(-1, camera.gaze);
        instance.v.x = camera.up.x;
        instance.v.y = camera.up.y;
        instance.v.z = camera.up.z;
        instance.u = crossProduct(instance.v, instance.w);
        instance.w = normalize(instance.w);
        instance.v = normalize(instance.v);
        instance.u = normalize(instance.u);
        instance.m = addVector(camera.position, scaleVector(-camera.near_distance, instance.w));
        instance.q = addVector(instance.m, addVector(scaleVector(viewPlane.x, instance.u), scaleVector(viewPlane.w, instance.v)));
        cameras.push_back(instance);

        for (Surface& surface : surfaces)
        {
            surface.relationship.push_back(getRelationship(camera, surface.triangle.indices));
        }
    }
}

int main(int argc, char* argv[])
{
    scene.loadFromXml(argv[1]);

    setup();

    for (int i = 0; i < cameras.size(); i++)
    {
        computeImage(i);
    }
}
