/**
* @file impurity.h
*
* @brief Header file for impurity.cpp
*/
#ifndef IMPURITY_H
#define IMPURITY_H

namespace Impurity
{
	/**
	* @brief Straightforward implementation of an impurity particle 
	*
	* Impurity particle class. Contains:
	*   t,x,y,z: Coordinates of the particle in space-time
	*   vx,vy,vz: Velocity components of the particle
	*   weight: The Monte Carlo weight of the particle
	*   charge: Current charge of the particle
	*   mass: The mass of the particle (kg)
	*/
	class Impurity
	{
	private:
		double m_t {};
		double m_x {};
		double m_y {};
		double m_z {};
		double m_vx {};
		double m_vy {};
		double m_vz {};
		double m_weight {};
		int m_charge {};
		double m_mass {};
		int m_atom_num {};
	
	public:

		/**
		* @brief Constructor: Start impurity with everything equal to zero
		*
		* This probably shouldn't be used?
		*/
		Impurity();
	
		/**
		* @brief Constructor: Fully specify each member variable
		*/
		Impurity(const double t, const double x, 
			const double y, const double z, const double vx, const double vy, 
			const double vz, const double weight, const int charge, 
			const double mass, const int atom_num);
		
		/**
		* @brief Accessor for particle weight (units depend, see imp_stats.cpp)
		* @return Returns particle weight as double
		*/
		double get_weight() const;
		/**
		* @brief Accessor for particle t (s)
		* @return Returns particle t as double
		*/
		double get_t() const;
		/**
		* @brief Accessor for particle x position (m)
		* @return Returns particle x position as double
		*/
		double get_x() const;
		/**
		* @brief Accessor for particle y position (m)
		* @return Returns particle y position as double
		*/
		double get_y() const;
		/**
		* @brief Accessor for particle z position (m)
		* @return Returns particle z position as double
		*/
		double get_z() const;
		/**
		* @brief Accessor for particle x velocity (m/s)
		* @return Returns particle x velocity as double
		*/
		double get_vx() const;
		/**
		* @brief Accessor for particle y velocity (m/s)
		* @return Returns particle y velocity as double
		*/
		double get_vy() const;
		/**
		* @brief Accessor for particle z velocity (m/s)
		* @return Returns particle z velocity as double
		*/
		double get_vz() const;
		/**
		* @brief Accessor for particle mass (kg)
		* @return Returns particle mass as double
		*/
		double get_mass() const;
		/**
		* @brief Accessor for particle charge (unitless)
		* @return Returns particle charge as int
		*/
		int get_charge() const;
		/**
		* @brief Accessor for atomic number (unitless)
		* @return Returns atomic number as int
		*/
		int get_atom_num() const;

		/**
		* @brief Setter for particle weight
		*/
		void set_weight(double w);
		/**
		* @brief Setter for particle t
		*/
		void set_t(double t);
		/**
		* @brief Setter for particle x position
		*/
		void set_x(double x);
		/**
		* @brief Setter for particle y position
		*/
		void set_y(double y);
		/**
		* @brief Setter for particle z position
		*/
		void set_z(double z);
		/**
		* @brief Setter for particle x velocity
		*/
		void set_vx(double vx);
		/**
		* @brief Setter for particle y velocity
		*/
		void set_vy(double vy);
		/**
		* @brief Setter for particle z velocity
		*/
		void set_vz(double vz);
		/**
		* @brief Setter for particle charge
		*/
		void set_charge(int charge);
	};
}

#endif
