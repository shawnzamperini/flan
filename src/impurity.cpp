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
		const double y, const double z, const double vx, const double vy, 
		const double vz, const double weight, const int charge, 
		const double mass, const int atom_num)
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
		, m_atom_num {atom_num}
	{}

	double Impurity::get_weight() const {return m_weight;}
	double Impurity::get_t() const {return m_t;}
	double Impurity::get_x() const {return m_x;}
	double Impurity::get_y() const {return m_y;}
	double Impurity::get_z() const {return m_z;}
	double Impurity::get_vx() const {return m_vx;}
	double Impurity::get_vy() const {return m_vy;}
	double Impurity::get_vz() const {return m_vz;}
	double Impurity::get_mass() const {return m_mass;}
	int Impurity::get_charge() const {return m_charge;}
	int Impurity::get_atom_num() const {return m_atom_num;}
		
	void Impurity::set_weight(double w) {m_weight = w;}
	void Impurity::set_t(double t) {m_t = t;}
	void Impurity::set_x(double x) {m_x = x;}
	void Impurity::set_y(double y) {m_y = y;}
	void Impurity::set_z(double z) {m_z = z;}
	void Impurity::set_vx(double vx) {m_vx = vx;}
	void Impurity::set_vy(double vy) {m_vy = vy;}
	void Impurity::set_vz(double vz) {m_vz = vz;}
	void Impurity::set_charge(int charge) {m_charge = charge;}
}
