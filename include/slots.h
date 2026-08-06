#pragma once

#include <vector>

#include "slots_device.h"

/*
#ifdef USE_CUDA
#include "slots_device.cuh"
#endif
*/

namespace Slots
{

	// Structure to hold initialized particle
	struct ParticleInit {
		double t, x, y, z;
		double vx, vy, vz;
		double vX, vY, vZ;
		double weight;
		int tidx, xidx, yidx, zidx, q, state;
	};

	/**
	* Impurity particle class. Contains:
	*   t: Time of the particle
	*   x,y,z: Coordinates of the particle in computational space
	*   X,Y,Z: Coordinates of the particle in Cartesian space
	*   vX,vY,vZ: Velocity components of the particle in Cartesian space
	*   weight: The Monte Carlo weight of the particle
	*   q: Current charge of the particle
	*   mass: The mass of the particle (kg)
	*/
    class Slots
    {
    private:
        int m_N;
		double m_mass;

        std::vector<double> m_t;
        std::vector<double> m_x;
        std::vector<double> m_y;
        std::vector<double> m_z;
        std::vector<double> m_vx;
        std::vector<double> m_vy;
        std::vector<double> m_vz;
        std::vector<double> m_vX;
        std::vector<double> m_vY;
        std::vector<double> m_vZ;

        std::vector<int> m_tidx;
        std::vector<int> m_xidx;
        std::vector<int> m_yidx;
        std::vector<int> m_zidx;

        std::vector<double> m_weight;
        std::vector<int> m_q;
        std::vector<int> m_state;

    public:
        // Constructor / Destructor
        explicit Slots(int N);
        ~Slots();

        // Accessors
        int N() const noexcept;
		double mass() const noexcept;
        const std::vector<double>& t() const noexcept;
        const std::vector<double>& x() const noexcept;
        const std::vector<double>& y() const noexcept;
        const std::vector<double>& z() const noexcept;
        const std::vector<double>& vx() const noexcept;
        const std::vector<double>& vy() const noexcept;
        const std::vector<double>& vz() const noexcept;
        const std::vector<double>& vX() const noexcept;
        const std::vector<double>& vY() const noexcept;
        const std::vector<double>& vZ() const noexcept;
        const std::vector<int>& tidx() const noexcept;
        const std::vector<int>& xidx() const noexcept;
        const std::vector<int>& yidx() const noexcept;
        const std::vector<int>& zidx() const noexcept;
        const std::vector<double>& weight() const noexcept;
        const std::vector<int>& q() const noexcept;
        const std::vector<int>& state() const noexcept;

		// Assume all particles share the same mass
		void set_mass(double mass);

        // Element-level setters
        void set_t(int i, double val);
        void set_x(int i, double val);
        void set_y(int i, double val);
        void set_z(int i, double val);
        void set_vx(int i, double val);
        void set_vy(int i, double val);
        void set_vz(int i, double val);
        void set_vX(int i, double val);
        void set_vY(int i, double val);
        void set_vZ(int i, double val);
        void set_tidx(int i, int val);
        void set_xidx(int i, int val);
        void set_yidx(int i, int val);
        void set_zidx(int i, int val);
        void set_weight(int i, double val);
        void set_q(int i, int val);
        void set_state(int i, int val);

        // GPU transfer
        SlotsDevice to_device(int device_id = 0);
        Slots to_host(const SlotsDevice& slots_d);
    };

	// GPU cleanup
	void free_slots(SlotsDevice& slots_d, int device_id = 0);

	// Replace dead particles with alive ones
	void fill_slots_cpu(Slots& slots, int& rem_parts);
}

