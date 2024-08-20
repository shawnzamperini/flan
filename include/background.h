#ifndef BACKGROUND_H
#define BACKGROUND_H

#include <vector>
#include "vectors.h"

namespace Background
{
	// Encapsulating class containing the Gkeyll background.
	class Background
	{
	private:
		std::vector<double> m_times {};
		Vectors::Vector4D m_ne {};
		Vectors::Vector4D m_te {};
		Vectors::Vector4D m_ti {};
		Vectors::Vector4D m_vp {};
		Vectors::Vector4D m_b {};

	public:

		// Accessors
		Vectors::Vector4D& get_ne() {return m_ne;}
		Vectors::Vector4D& get_te() {return m_te;}
		Vectors::Vector4D& get_ti() {return m_ti;}
		Vectors::Vector4D& get_vp() {return m_vp;}
		Vectors::Vector4D& get_b() {return m_b;}

		// Setters using move semantics. No idea if this is the proper
		// way to do it but it works.
		void move_into_ne(Vectors::Vector4D& ne) 
		{
			m_ne.move_into_data(ne);
		}
		void move_into_te(Vectors::Vector4D& te) 
		{
			m_te.move_into_data(te);
		}
		void move_into_ti(Vectors::Vector4D& ti) 
		{
			m_ti.move_into_data(ti);
		}
		void move_into_vp(Vectors::Vector4D& vp) 
		{
			m_vp.move_into_data(vp);
		}
		void move_into_b(Vectors::Vector4D& b) 
		{
			m_b.move_into_data(b);
		}
	};
}

#endif
