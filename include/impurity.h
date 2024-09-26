/**
* @file impurity.h
* @brief Header file for impurity.cpp
*/
#ifndef IMPURITY_H
#define IMPURITY_H

namespace Impurity
{
	// Impurity particle class. Contains:
	// t,x,y,z: Coordinates of the particle in space-time
	// vx,vy,vz: Velocity components of the particle
	// weight: The Monte Carlo weight of the particle
	// charge: Current charge of the particle
	// mass: The mass of the particle (kg)
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
	
	public:

		// Constructor: Start impurity at rest and at (0,0,0)
		Impurity();
	
		// Constructor: Fully specify each member variable
		Impurity(const double t, const double x, 
			const double y, const double z, const double vx, const double vy, 
			const double vz, const double weight, const int charge, 
			const double mass);
		
		// Accessors
		double get_weight() const;
		double get_t() const;
		double get_x() const;
		double get_y() const;
		double get_z() const;
		double get_vx() const;
		double get_vy() const;
		double get_vz() const;
		double get_mass() const;
		int get_charge() const;

		// Setters
		void set_weight(double w);
		void set_t(double t);
		void set_x(double x);
		void set_y(double y);
		void set_z(double z);
		void set_vx(double vx);
		void set_vy(double vy);
		void set_vz(double vz);
		void set_charge(int charge);
	};
}

#endif
