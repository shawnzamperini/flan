/**
* @file moeller_trumbore.cpp
*/
#include <iostream>
#include <optional>
#include <limits>

#include "moeller_trumbore.h"

namespace MoellerTrumbore
{
    // Dot product
    double dot(const vec3& v1, const vec3& v2)
    {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }

    // Cross product
    vec3 cross(const vec3& v1, const vec3& v2)
    {
        return {
            v1.y * v2.z - v1.z * v2.y,
            v1.z * v2.x - v1.x * v2.z,
            v1.x * v2.y - v1.y * v2.x};
    }

    // Test for intersection between line segment and triangle. Returns the
    // fraction of how far along the line segment the intersection occurred,
    // if it intersects. Return -1 if not.
    double ray_intersects_triangle(const vec3 &ray_origin,
        const vec3 &ray_vector, const triangle3& triangle)
    {
        constexpr double epsilon = std::numeric_limits<double>::epsilon();

        vec3 edge1 = triangle.b - triangle.a;
        vec3 edge2 = triangle.c - triangle.a;
        vec3 ray_cross_e2 = cross(ray_vector, edge2);
        double det = dot(edge1, ray_cross_e2);

        if (det > -epsilon && det < epsilon)
            return -1;    // This ray is parallel to this triangle.

        double inv_det = 1.0 / det;
        vec3 s = ray_origin - triangle.a;
        double u = inv_det * dot(s, ray_cross_e2);

        if ((u < 0 && abs(u) > epsilon) || (u > 1 && abs(u-1) > epsilon))
            return -1;

        vec3 s_cross_e1 = cross(s, edge1);
        double v = inv_det * dot(ray_vector, s_cross_e1);

        if ((v < 0 && abs(v) > epsilon) || (u + v > 1 && abs(u + v - 1) 
            > epsilon)) return -1;

        // At this stage we can compute t to find out where the intersection 
        // point is on the line.
        double t = inv_det * dot(edge2, s_cross_e1);

        // MODIFICATION OF ORIGINAL CODE
        // The normal code test for intersection between an infinitely long
        // ray (t --> inf) and a triangle. But we actually want to see if a 
        // finite line segment (that represents a particle's trajectory over
        // a finite time step) crosses through the triangle. The line is
        // parameterized as r(t) = O + tv, where O is an origin point, v is
        // the direction of the vector, and t is the parameterization variable.
        // So our finite segment is the range t=[0,1]. Thus we want to check
        // that if 0 < t < 1, the line segment intersects. If not, then it
        // does not intersect.
        if (t < 0 || t > 1) return -1;

        if (t > epsilon) // ray intersection
        {
            //return  vec3(ray_origin + ray_vector * t);
            return t;
        }

        // This means that there is a line intersection but not a ray 
        // intersection.
        else return -1;
    }

    // Check if line segement defined by P1->P2 passes through a surface, 
    // returning a number between 0-1 if so that represents at what fraction 
    // of the segment the intersection occurred.  
    double get_intersect_frac(
        const double p1x, const double p1y, const double p2z,
        const double p2x, const double p2y, const double p1z, 
        const double v1x, const double v1y, const double v1z, 
        const double v2x, const double v2y, const double v2z, 
        const double v3x, const double v3y, const double v3z, 
        const double v4x, const double v4y, const double v4z)
    {
        // Break surface up into two triangles. v1->v2->v3 and v3->v4->v1.
        triangle3 tri1 {{v1x, v1y, v1z}, {v2x, v2y, v3z}, {v3x, v3y, v3z}};
        triangle3 tri2 {{v3x, v3y, v3z}, {v4x, v4y, v4z}, {v1x, v1y, v1z}};

        // Same triangles, but in reverse order
        //triangle3 tri3 {{v3x, v3y, v3z}, {v2x, v2y, v3z}, {v1x, v1y, v1z}};
        //triangle3 tri4 {{v1x, v1y, v1z}, {v4x, v4y, v4z}, {v3x, v3y, v3z}};

        // Ray origin is p1, ray vector goes to p2
        vec3 ray_origin {p1x, p1y, p1z};
        vec3 ray_vec {p2x-p1x, p2y-p1y, p2z-p1z};

        /*
        std::cout << "coordinates\n"
            << "  " << v1x << " " << v1y << " " << v1z << "\n"
            << "  " << v2x << " " << v2y << " " << v2z << "\n"
            << "  " << v3x << " " << v3y << " " << v3z << "\n"
            << "  " << v4x << " " << v4y << " " << v4z << "\n";
        */
            
        // Check intersection with first triangle
        double intersect_frac {ray_intersects_triangle(ray_origin, ray_vec, 
            tri1)};
        if (intersect_frac > 0) 
        {
            //std::cout << "intersection! tri1 " << intersect_frac << '\n';
        }

        // Returns -1 if not intersection, so check second triangle
        if (intersect_frac < 0) 
        {
            intersect_frac = ray_intersects_triangle(ray_origin, ray_vec, 
                tri2);
            //if (intersect_frac > 0) std::cout << "intersection! tri2 " 
            //    << intersect_frac << '\n';
        }
/*
        // Returns -1 if not intersection, so check second triangle
        if (intersect_frac < 0) 
        {
            intersect_frac = ray_intersects_triangle(ray_origin, ray_vec, 
                tri3);
            //if (intersect_frac > 0) std::cout << "intersection! tri3 " 
            //    << intersect_frac << '\n';
        }

        // Returns -1 if not intersection, so check second triangle
        if (intersect_frac < 0) 
        {
            intersect_frac = ray_intersects_triangle(ray_origin, ray_vec, 
                tri4);
            //if (intersect_frac > 0) std::cout << "intersection! tri4 " 
            //    << intersect_frac << '\n';
        }
*/
        // This will either have the intersect_frac, or -1 to indicate no
        // intersection.
        return intersect_frac;
    }
}
