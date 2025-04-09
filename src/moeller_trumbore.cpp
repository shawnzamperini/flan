/**
* @file moeller_trumbore.cpp
*/
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>

#include "moeller_trumbore.h"

namespace MoellerTrumbore
{
	// Overload the operator<< for vec3
	std::ostream& operator<<(std::ostream& os, const vec3& v) 
	{
		os << "(" << v.x << ", " << v.y << ", " << v.z << ")"
			<< "  R = " << std::sqrt(v.x*v.x + v.y*v.y);
		return os;
	}

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

	// Test for intersection between a ray and a triangle
    bool ray_intersects_triangle(const vec3 &ray_origin,
        const vec3 &ray_vector, const triangle3& triangle, const bool debug)
    {
        constexpr double epsilon = std::numeric_limits<double>::epsilon();

		if (debug) std::cout << "ray_origin = " << ray_origin << '\n';
		if (debug) std::cout << "ray_vector = " << ray_vector << '\n';
		if (debug) std::cout << "triangle.a = " << triangle.a << '\n';
		if (debug) std::cout << "triangle.b = " << triangle.b << '\n';
		if (debug) std::cout << "triangle.c = " << triangle.c << '\n';
        vec3 edge1 = triangle.b - triangle.a;    // E1
		if (debug) std::cout << "edge1 = " << edge1 << '\n';
        vec3 edge2 = triangle.c - triangle.a;    // E2
		if (debug) std::cout << "edge2 = " << edge2 << '\n';
        vec3 ray_cross_e2 = cross(ray_vector, edge2);  // P
		if (debug) std::cout << "ray_cross_e2 = " << ray_cross_e2 << '\n';
        //double det = dot(edge1, ray_cross_e2);  // E1 dot P
        double det = std::abs(dot(edge1, ray_cross_e2));  // E1 dot P

		if (debug) std::cout << "det = " << det << '\n';
        if (det > -epsilon && det < epsilon)
		{
			if (debug) std::cout << "parallel\n";
            return false;    // This ray is parallel to this triangle.
		}

        double inv_det = 1.0 / det;
        vec3 s = ray_origin - triangle.a;  // T
		if (debug) std::cout << "s = " << s << '\n';
        double u = inv_det * dot(s, ray_cross_e2);  // T dot P
		if (debug) std::cout << "u = " << u << '\n';

        if ((u < 0 && std::abs(u) > epsilon) || 
			(u > 1 && std::abs(u-1) > epsilon))
		{
			if (debug) std::cout << "barycentric fail #1...\n";
            return false;
		}

        vec3 s_cross_e1 = cross(s, edge1);  // Q = T x E1
		if (debug) std::cout << "s_cross_e1 = " << s_cross_e1 << '\n';
        double v = inv_det * dot(ray_vector, s_cross_e1);  // Q dot D
		if (debug) std::cout << "v = " << v << '\n';

		// I note a warning here of a particularly hard bug to find. This
		// code was copied from Wikipedia, but it had just abs and not
		// std::abs. abs is a compiler command and casts the number to an 
		// integer, depending on the compiler. This messes everything up,
		// but I eventually caught it after much pain.
        if ((v < 0 && std::abs(v) > epsilon) || 
			(u + v > 1 && std::abs(u + v - 1) > epsilon)) 
		{
			if (debug) std::cout << "barycentric fail #2...\n";
			return false;
		}

        // At this stage we can compute t to find out where the intersection 
        // point is on the line.
        double t = inv_det * dot(edge2, s_cross_e1);  // E2 dot Q
		if (debug) std::cout << "t = " << t << '\n';

		if (t < 0 || t > 1)
		{
			if (debug) std::cout << "t not between [0,1]: " << t << '\n';
			return false;
		}
        if (t > epsilon) // ray intersection
        {
            //return  vec3(ray_origin + ray_vector * t);
			if (debug) std::cout << "(u,v,w) = " << "(" << u << ", " << v 
				<< ", "	<< (1-u-v) << ")\n";

			// I don't know why this is happeneing, but I am adding an
			// additional check. The algorithm will return a positive t for
			// surfaces that are in the opposite direction of the ray, when
			// it should return a negative value to indicate this. But since 
			// we now know the coordinate at which an intersection occurs,
			// we want to see if it is past the point along the ray. 
			//vec3 int_pt {triangle.a * (1 - u - v) + triangle.b * u 
			//	+ triangle.c * v};
			//std::cout << "int_pt = " << int_pt << '\n';

			// If it's magnitude is greater than that where to origin is, the
			// intersection is past and thus not really an intersection since
			// the line stops at the origin. Can just compare d^2.
			//double origin_mag {ray_origin.x * ray_origin.x
			//	+ ray_origin.y * ray_origin.y + ray_origin.z * ray_origin.z};
			//double int_mag {int_pt.x * int_pt.x + int_pt.y * int_pt.y 
			//	+ int_pt.z * int_pt.z};
			//std::cout << "origin_mag = " << origin_mag << "  int_mag = " 
			//	<< int_mag << '\n';
			//if (int_mag > origin_mag)
			//{
			//	std::cout << "Not really an intersection.\n";
			//	return false;
			//}

            return true;
        }

        // This means that there is a line intersection but not a ray 
        // intersection.
        else return false;
    }

    // Test for intersection between line segment and triangle. Returns the
    // fraction of how far along the line segment the intersection occurred,
    // if it intersects. Return -1 if not.
    double seg_intersects_triangle(const vec3 &ray_origin,
        const vec3 &ray_vector, const triangle3& triangle)
    {
        constexpr double epsilon = std::numeric_limits<double>::epsilon();

        vec3 edge1 = triangle.b - triangle.a;
        vec3 edge2 = triangle.c - triangle.a;
        vec3 ray_cross_e2 = cross(ray_vector, edge2);
        double det = dot(edge1, ray_cross_e2);

        if (det > -epsilon && det < epsilon)
		{
			//std::cout << "parallel\n";
            return -1;    // This ray is parallel to this triangle.
		}

        double inv_det = 1.0 / det;
        vec3 s = ray_origin - triangle.a;
        double u = inv_det * dot(s, ray_cross_e2);
		//std::cout << "u = " << u << '\n';

        if ((u < 0 && std::abs(u) > epsilon) || 
			(u > 1 && std::abs(u-1) > epsilon))
		{
			//std::cout << "barycentric fail #1...\n";
            return -1;
		}

        vec3 s_cross_e1 = cross(s, edge1);
        double v = inv_det * dot(ray_vector, s_cross_e1);
		//std::cout << "v = " << v << '\n';

        if ((v < 0 && std::abs(v) > epsilon) || 
			(u + v > 1 && std::abs(u + v - 1) > epsilon)) 
		{
			//std::cout << "barycentric fail #2...\n";
			return -1;
		}


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
		//std::cout << "t = " << t << '\n';
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

    bool check_intersect(
        const double p1x, const double p1y, const double p1z,
        const double p2x, const double p2y, const double p2z, 
        const double v1x, const double v1y, const double v1z, 
        const double v2x, const double v2y, const double v2z, 
        const double v3x, const double v3y, const double v3z, 
        const double v4x, const double v4y, const double v4z,
		const bool debug)
    {
        // Break surface up into two triangles. v1->v2->v3 and v3->v4->v1.
        //triangle3 tri1 {{v1x, v1y, v1z}, {v2x, v2y, v2z}, {v3x, v3y, v3z}};
        //triangle3 tri2 {{v3x, v3y, v3z}, {v4x, v4y, v4z}, {v1x, v1y, v1z}};
        triangle3 tri1 {{v3x, v3y, v3z}, {v2x, v2y, v2z}, {v1x, v1y, v1z}};
        triangle3 tri2 {{v1x, v1y, v1z}, {v4x, v4y, v4z}, {v3x, v3y, v3z}};

        // Same triangles, but in reverse order
        //triangle3 tri3 {{v3x, v3y, v3z}, {v2x, v2y, v2z}, {v1x, v1y, v1z}};
        //triangle3 tri4 {{v1x, v1y, v1z}, {v4x, v4y, v4z}, {v3x, v3y, v3z}};

        // Ray origin is p1, ray vector goes to p2
        vec3 ray_origin {p1x, p1y, p1z};
        vec3 ray_vec {p2x-p1x, p2y-p1y, p2z-p1z};
		if (debug) std::cout << "ray_origin = " << ray_origin << '\n';
		if (debug) std::cout << "ray_vec = " << ray_vec << '\n';

		// Must normalize!!!
		ray_vec = ray_vec.normalize();

        // Check intersection with first triangle
        bool intersect {ray_intersects_triangle(ray_origin, ray_vec, 
            tri1, debug)};

        // Returns false if not intersection, so check second triangle
        if (!intersect) 
        {
            intersect = ray_intersects_triangle(ray_origin, ray_vec, 
                tri2, debug);
        }

		/*
        if (intersect_frac < 0) 
        {
            intersect_frac = seg_intersects_triangle(ray_origin, ray_vec, 
                tri3);
            //if (intersect_frac > 0) std::cout << "intersection! tri3 " 
            //    << intersect_frac << '\n';
        }

        if (intersect_frac < 0) 
        {
            intersect_frac = seg_intersects_triangle(ray_origin, ray_vec, 
                tri4);
            //if (intersect_frac > 0) std::cout << "intersection! tri4 " 
            //    << intersect_frac << '\n';
        }
		*/

		// This will return false if we didn't intersect with either of the
		// triangles that make up the surface, true if so
        return intersect;
    }

    // Check if line segement defined by P1->P2 passes through a surface, 
    // returning a number between 0-1 if so that represents at what fraction 
    // of the segment the intersection occurred.  
    double get_intersect_frac(
        const double p1x, const double p1y, const double p1z,
        const double p2x, const double p2y, const double p2z, 
        const double v1x, const double v1y, const double v1z, 
        const double v2x, const double v2y, const double v2z, 
        const double v3x, const double v3y, const double v3z, 
        const double v4x, const double v4y, const double v4z)
    {
        // Break surface up into two triangles. v1->v2->v3 and v3->v4->v1.
        //triangle3 tri1 {{v1x, v1y, v1z}, {v2x, v2y, v2z}, {v3x, v3y, v3z}};
        //triangle3 tri2 {{v3x, v3y, v3z}, {v4x, v4y, v4z}, {v1x, v1y, v1z}};
        triangle3 tri1 {{v3x, v3y, v3z}, {v2x, v2y, v2z}, {v1x, v1y, v1z}};
        triangle3 tri2 {{v1x, v1y, v1z}, {v4x, v4y, v4z}, {v3x, v3y, v3z}};

        // Same triangles, but in reverse order
        //triangle3 tri3 {{v3x, v3y, v3z}, {v2x, v2y, v2z}, {v1x, v1y, v1z}};
        //triangle3 tri4 {{v1x, v1y, v1z}, {v4x, v4y, v4z}, {v3x, v3y, v3z}};

        // Ray origin is p1, ray vector goes to p2
        vec3 ray_origin {p1x, p1y, p1z};
        vec3 ray_vec {p2x-p1x, p2y-p1y, p2z-p1z};

		// Must normalize!!!
		ray_vec = ray_vec.normalize();

        /*
        std::cout << "coordinates\n"
            << "  " << v1x << " " << v1y << " " << v1z << "\n"
            << "  " << v2x << " " << v2y << " " << v2z << "\n"
            << "  " << v3x << " " << v3y << " " << v3z << "\n"
            << "  " << v4x << " " << v4y << " " << v4z << "\n";
        */
            
        // Check intersection with first triangle
        double intersect_frac {seg_intersects_triangle(ray_origin, ray_vec, 
            tri1)};
        if (intersect_frac > 0) 
        {
            //std::cout << "intersection! tri1 " << intersect_frac << '\n';
        }

        // Returns -1 if not intersection, so check second triangle
        if (intersect_frac < 0) 
        {
            intersect_frac = seg_intersects_triangle(ray_origin, ray_vec, 
                tri2);
            //if (intersect_frac > 0) std::cout << "intersection! tri2 " 
            //    << intersect_frac << '\n';
        }
		/*
        // Returns -1 if not intersection, so check second triangle
        if (intersect_frac < 0) 
        {
            intersect_frac = seg_intersects_triangle(ray_origin, ray_vec, 
                tri3);
            //if (intersect_frac > 0) std::cout << "intersection! tri3 " 
            //    << intersect_frac << '\n';
        }

        // Returns -1 if not intersection, so check second triangle
        if (intersect_frac < 0) 
        {
            intersect_frac = seg_intersects_triangle(ray_origin, ray_vec, 
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
