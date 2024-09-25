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
		Impurity()
		{}
	
		// Constructor: Start impurity at rest
		//Impurity(double t, double x, double y, double z)
		//	: m_t {t}
		//	, m_y {y}
		//	, m_y {y}
		//	, m_z {z}
		//{}

		// Constructor: Fully specify each member variable
		Impurity(const double t, const double x, 
			const double y, const double z, const double vx, const double vy, 
			const double vz, const double weight, const int charge, 
			const double mass)
			: m_t {t}
			, m_x {x}
			, m_y {y}
			, m_z {z}
			, m_vx {vx}
			, m_vy {vy}
			, m_vz {vz} 
			, m_weight {weight}
			, m_charge {charge}
			, m_mass {mass}
		{}
		
		// Accessors
		double get_weight() const {return m_weight;}
		double get_t() const {return m_t;}
		double get_x() const {return m_x;}
		double get_y() const {return m_y;}
		double get_z() const {return m_z;}
		double get_vx() const {return m_vx;}
		double get_vy() const {return m_vy;}
		double get_vz() const {return m_vz;}
		double get_mass() const {return m_mass;}
		int get_charge() const {return m_charge;}

		// Setters
		void set_weight(double w) {m_weight = w;}
		void set_t(double t) {m_t = t;}
		void set_x(double x) {m_x = x;}
		void set_y(double y) {m_y = y;}
		void set_z(double z) {m_z = z;}
		void set_vx(double vx) {m_vx = vx;}
		void set_vy(double vy) {m_vy = vy;}
		void set_vz(double vz) {m_vz = vz;}
		void set_charge(int charge) {m_charge = charge;}
	};
}

#endif
