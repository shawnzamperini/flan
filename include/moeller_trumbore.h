#ifndef MOELLER_TRUMBORE_H
#define MOELLER_TRUMBORE_H

#include <cmath>

/**
* @file moeller_trumbore.h
*
* @brief Implementation of the Moeller-Trumbore algorithm for testing
* intersection between a vector and a 3D triangular surface. 
*
* The main function ray_intersects_triangle comes from the Wikipedia page
* for Moeller-Trumbore intersection algorithm. I just had to define vec3 and
* triangle3 since those were not included.
*/

namespace MoellerTrumbore
{
    /**
    * @brief Simple implementation of a 3D vector defined by three doubles
    */
    struct vec3
    {
        // Coordinates
        double x, y, z;

        // Addition operator to subtract two vectors
        vec3 operator+(const vec3& other) const
        {
            return {x + other.x, y + other.y, z + other.z};
        }

        // Subtraction operator to subtract two vectors
        vec3 operator-(const vec3& other) const
        {
            return {x - other.x, y - other.y, z - other.z};
        }

        // Multiplication to mutiply by a scalar
        vec3 operator*(const double scalar) const
        {
            return {x * scalar, y * scalar, z * scalar};
        }
			
		// Normalize function
		vec3 normalize() const 
		{
			double magnitude = std::sqrt(x * x + y * y + z * z);
			return {x / magnitude, y / magnitude, z / magnitude};
		}
    };

    /**
    * @brief A 3D triangle defined by a trio of vec3 (x,y,z coordinates)
    */
    struct triangle3
    {
        // The three coordinates that make up the triangle
        vec3 a, b, c;
    };
    
    /**
    * @brief Calculate dot product for two vec3.
    */
    double dot(const vec3& v1, const vec3& v2);

    /**
    * @brief Calculate cross product for two vec3.
    */
    vec3 cross(const vec3& v1, const vec3& v2);

	/**
	* @brief Check if ray intersects with the quadrilateral surface defined 
	* by v1, v2, v3, v4.
	*
	* @return Retuns true if it does, false if not
	*/
    bool check_intersect(
        const double p1x, const double p1y, const double p1z,
        const double p2x, const double p2y, const double p2z, 
        const double v1x, const double v1y, const double v1z, 
        const double v2x, const double v2y, const double v2z, 
        const double v3x, const double v3y, const double v3z, 
        const double v4x, const double v4y, const double v4z,
		const bool debug=false);

    /**
    * @brief Get the fraction along the line segment from p1-->p2 at which it
    * intersects with the quadrilateral surface defined by v1, v2, v3, v4.
    *
    * @return Returns fraction between 0 and 1 if intersection, -1 otherwise.
    */
    double get_intersect_frac(
        const double p1x, const double p1y, const double p1z,
        const double p2x, const double p2y, const double p2z, 
        const double v1x, const double v1y, const double v1z, 
        const double v2x, const double v2y, const double v2z, 
        const double v3x, const double v3y, const double v3z, 
        const double v4x, const double v4y, const double v4z);
}

#endif
