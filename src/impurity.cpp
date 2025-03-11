/**
* @file impurity.cpp
*
* @brief Contains a straightforward implementation of an impurity particle
*/
#include "impurity.h"

namespace Impurity
{
	Impurity::Impurity()
	{}

	Impurity::Impurity(const double t, const double x, 
		const double y, const double z, const double X, const double Y, 
		const double Z, const double vX, const double vY, 
		const double vZ, const double weight, const int charge, 
		const double mass, const int atom_num)
		: m_t {t}
		, m_x {x}
		, m_y {y}
		, m_z {z}
		, m_X {X}
		, m_Y {Y}
		, m_Z {Z}
		, m_vX {vX}
		, m_vY {vY}
		, m_vZ {vZ} 
		, m_weight {weight}
		, m_charge {charge}
		, m_mass {mass}
		, m_atom_num {atom_num}
	{}

	double Impurity::get_weight() const {return m_weight;}
	double Impurity::get_t() const {return m_t;}
	double Impurity::get_x() const {return m_x;}
	double Impurity::get_y() const {return m_y;}
	double Impurity::get_z() const {return m_z;}
	double Impurity::get_X() const {return m_X;}
	double Impurity::get_Y() const {return m_Y;}
	double Impurity::get_Z() const {return m_Z;}
	double Impurity::get_vX() const {return m_vX;}
	double Impurity::get_vY() const {return m_vY;}
	double Impurity::get_vZ() const {return m_vZ;}
	double Impurity::get_mass() const {return m_mass;}
	int Impurity::get_charge() const {return m_charge;}
	int Impurity::get_atom_num() const {return m_atom_num;}
		
	void Impurity::set_weight(double w) {m_weight = w;}
	void Impurity::set_t(double t) {m_t = t;}
	void Impurity::set_x(double x) {m_x = x;}
	void Impurity::set_y(double y) {m_y = y;}
	void Impurity::set_z(double z) {m_z = z;}
	void Impurity::set_X(double X) {m_X = X;} void Impurity::set_Y(double Y) {m_Y = Y;}
	void Impurity::set_Z(double Z) {m_Z = Z;}
	void Impurity::set_vX(double vX) {m_vX = vX;}
	void Impurity::set_vY(double vY) {m_vY = vY;}
	void Impurity::set_vZ(double vZ) {m_vZ = vZ;}
	void Impurity::set_charge(int charge) {m_charge = charge;}
}
