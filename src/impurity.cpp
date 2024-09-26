/**
* @file impurity.cpp
*
* @brief Contains a straightforward implementation of an impurity particle
*/
#include "impurity.h"

namespace Impurity
{
	/**
	* @class Impurity
	* @brief Straightforward implementation of an impurity particle
	*/

	/**
	* @brief Constructor: Start impurity with everything equal to zero
	*
	* This probably shouldn't be used?
	*/
	Impurity::Impurity()
	{}

	/**
	* @brief Constructor: Fully specify each member variable
	*/
	Impurity::Impurity(const double t, const double x, 
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

	/**
	* @brief Accessor for particle weight (units depend, see imp_stats.cpp)
	* @return Returns particle weight as double
	*/
	double Impurity::get_weight() const {return m_weight;}
	/**
	* @brief Accessor for particle t (s)
	* @return Returns particle t as double
	*/
	double Impurity::get_t() const {return m_t;}
	/**
	* @brief Accessor for particle x position (m)
	* @return Returns particle x position as double
	*/
	double Impurity::get_x() const {return m_x;}
	/**
	* @brief Accessor for particle y position (m)
	* @return Returns particle y position as double
	*/
	double Impurity::get_y() const {return m_y;}
	/**
	* @brief Accessor for particle z position (m)
	* @return Returns particle z position as double
	*/
	double Impurity::get_z() const {return m_z;}
	/**
	* @brief Accessor for particle x velocity (m/s)
	* @return Returns particle x velocity as double
	*/
	double Impurity::get_vx() const {return m_vx;}
	/**
	* @brief Accessor for particle y velocity (m/s)
	* @return Returns particle y velocity as double
	*/
	double Impurity::get_vy() const {return m_vy;}
	/**
	* @brief Accessor for particle z velocity (m/s)
	* @return Returns particle z velocity as double
	*/
	double Impurity::get_vz() const {return m_vz;}
	/**
	* @brief Accessor for particle mass (kg)
	* @return Returns particle mass as double
	*/
	double Impurity::get_mass() const {return m_mass;}
	/**
	* @brief Accessor for particle charge (unitless)
	* @return Returns particle charge as int
	*/
	int Impurity::get_charge() const {return m_charge;}
		
	/**
	* @brief Setter for particle weight
	*/
	void Impurity::set_weight(double w) {m_weight = w;}
	/**
	* @brief Setter for particle t
	*/
	void Impurity::set_t(double t) {m_t = t;}
	/**
	* @brief Setter for particle x position
	*/
	void Impurity::set_x(double x) {m_x = x;}
	/**
	* @brief Setter for particle y position
	*/
	void Impurity::set_y(double y) {m_y = y;}
	/**
	* @brief Setter for particle z position
	*/
	void Impurity::set_z(double z) {m_z = z;}
	/**
	* @brief Setter for particle x velocity
	*/
	void Impurity::set_vx(double vx) {m_vx = vx;}
	/**
	* @brief Setter for particle y velocity
	*/
	void Impurity::set_vy(double vy) {m_vy = vy;}
	/**
	* @brief Setter for particle z velocity
	*/
	void Impurity::set_vz(double vz) {m_vz = vz;}
	/**
	* @brief Setter for particle charge
	*/
	void Impurity::set_charge(int charge) {m_charge = charge;}
}
