#include "primitive.h"
#include <math.h>

// Ray-sphere intersection
bool sphere_hit(const Sphere* sphere, const Ray* ray, float t_min, float t_max,
                HitRecord* rec) {
    // Vector from ray origin to sphere center
    Vec3 oc = vec3_sub(ray->origin, sphere->center);
    
    // Quadratic equation coefficients: at² + bt + c = 0
    float a = vec3_dot(ray->direction, ray->direction);
    float half_b = vec3_dot(oc, ray->direction);
    float c = vec3_dot(oc, oc) - sphere->radius * sphere->radius;
    
    // Calculate discriminant
    float discriminant = half_b * half_b - a * c;
    
    // No intersection if discriminant is negative
    if (discriminant < 0) {
        return false;
    }
    
    float sqrtd = sqrtf(discriminant);
    
    // Find the nearest root that lies in the valid range [t_min, t_max]
    float root = (-half_b - sqrtd) / a;
    if (root < t_min || root > t_max) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || root > t_max) {
            return false;
        }
    }
    
    // Fill hit record
    rec->t = root;
    rec->point = ray_at(*ray, rec->t);
    Vec3 outward_normal = vec3_div(vec3_sub(rec->point, sphere->center), sphere->radius);
    
    // Determine if ray hit front or back face
    rec->front_face = vec3_dot(ray->direction, outward_normal) < 0;
    rec->normal = rec->front_face ? outward_normal : vec3_scale(outward_normal, -1.0f);
    
    return true;
}

// Ray-triangle intersection (Möller-Trumbore algorithm)
bool triangle_hit(const Triangle* triangle, const Ray* ray, float t_min, float t_max,
                  HitRecord* rec) {
    const float EPSILON = 0.0000001f;

    // Compute edge vectors
    Vec3 edge1 = vec3_sub(triangle->v1, triangle->v0);
    Vec3 edge2 = vec3_sub(triangle->v2, triangle->v0);
    
    // Compute h = ray.direction × edge2
    Vec3 h = vec3_cross(ray->direction, edge2);
    float a = vec3_dot(edge1, h);
    
    // Check if ray is parallel to triangle
    if (fabsf(a) < EPSILON) {
        return false;
    }
    
    float f = 1.0f / a;
    Vec3 s = vec3_sub(ray->origin, triangle->v0);
    
    // Compute barycentric coordinate u
    float u = f * vec3_dot(s, h);
    if (u < 0.0f || u > 1.0f) {
        return false;
    }
    
    // Compute q = s × edge1
    Vec3 q = vec3_cross(s, edge1);
    
    // Compute barycentric coordinate v
    float v = f * vec3_dot(ray->direction, q);
    if (v < 0.0f || u + v > 1.0f) {
        return false;
    }
    
    // Compute t parameter
    float t = f * vec3_dot(edge2, q);
    
    // Check if t is within valid range
    if (t < t_min || t > t_max) {
        return false;
    }
    
    // Fill hit record
    rec->t = t;
    rec->point = ray_at(*ray, rec->t);
    
    // Determine front/back face
    Vec3 outward_normal = triangle->normal;
    rec->front_face = vec3_dot(ray->direction, outward_normal) < 0;
    rec->normal = rec->front_face ? outward_normal : vec3_scale(outward_normal, -1.0f);
    
    return true;
}

// Generic primitive hit test
bool primitive_hit(const Primitive* prim, const Ray* ray, float t_min, float t_max,
                   HitRecord* rec) {
    bool hit = false;

    switch (prim->type) {
        case PRIMITIVE_SPHERE:
            hit = sphere_hit(&prim->sphere, ray, t_min, t_max, rec);
            break;
        case PRIMITIVE_TRIANGLE:
            hit = triangle_hit(&prim->triangle, ray, t_min, t_max, rec);
            break;
        default:
            return false;
    }

    if (hit) {
        rec->material = &prim->material;
    }

    return hit;
}