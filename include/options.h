/**
* @file options.h
*/
#ifndef options_h
#define options_h

namespace Options
{
	class Options
	{
	private:

		// Input options and their defaults are defined here
		int m_imp_num {1};

	public:

		Options();

		/**
		* @brief Setter for imp_num
		*/
		void set_imp_num(int imp_num);

		/**
		* @brief Accessor for imp_num
		*/
		int imp_num();
	};
}

#endif
