
namespace Impurity
{
	class Impurity
	{
	private:
		double m_weight {};
		double m_x {};
		double m_y {};
		double m_z {};
		double m_vx {};
		double m_vy {};
		double m_vz {};
	
	public:

		// Constructor: Start impurity at rest
		Impurity(double x, double y, double z)
			: m_x {x}
			, m_y {y}
			, m_z {z}
		{}
		
		// Accessors
		double get_weight() {return m_weight;}
		double get_x() {return m_x;}
		double get_y() {return m_y;}
		double get_z() {return m_z;}
		double get_vx() {return m_vx;}
		double get_vy() {return m_vy;}
		double get_vz() {return m_vz;}

		// Setters
		void set_weight(double w) {m_weight = w;}
		void set_x(double x) {m_x = x;}
		void set_y(double y) {m_y = y;}
		void set_z(double z) {m_z = z;}
		void set_vx(double vx) {m_vx = vx;}
		void set_vy(double vy) {m_vy = vy;}
		void set_vz(double vz) {m_vz = vz;}
	};
}
