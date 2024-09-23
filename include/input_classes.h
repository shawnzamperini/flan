/**
* @file input_classes.h
* @brief input_classes.cpp header file
*/

#ifndef INPUT_CLASSES_H
#define INPUT_CLASSES_H


namespace Input
{
	
	// Base class to derive different input option types (int, double. string, ...)
	// from.  
	class OptionBase
	{
	private:
		std::string m_name {};
	
	public:
		OptionBase(std::string name);

		virtual ~OptionBase() = default;
		
		std::string_view get_name() const;

	};

	// Integer input option class.
	class OptionInt : public OptionBase
	{
	private:
		int m_value {0};
	public:

		// Constructor
		OptionInt(std::string name, int value);

		// Getter and putter
		int get_value() const;

		// Overloading to allow string input
		void put_value(int value);
		void put_value(std::string value);
		void put_value(std::string_view sv);
	};

	// Integer input option class.
	class OptionDbl : public OptionBase
	{
	private:
		double m_value {0};
	public:

		// Constructor
		OptionDbl(std::string name, double value);

		// Getter and putter
		double get_value() const;

		// Overloading to allow string input
		void put_value(double value);
		void put_value(std::string value);
		void put_value(std::string_view sv);
	};

	// String input option class.
	class OptionStr : public OptionBase
	{
	private:
		std::string m_value {"undefined"};
	public:
		
		// Constructor
		OptionStr(std::string name, std::string value);

		// Getter and putter
		std::string_view get_value() const;
		void put_value(std::string_view value);
	};

	// Overload << so we can easily print an option out in a nice format
	std::ostream& operator<<(std::ostream& out, const OptionStr opt);
	std::ostream& operator<<(std::ostream& out, const OptionInt opt);
	std::ostream& operator<<(std::ostream& out, const OptionDbl opt);
	
	template <typename T>
	T sv_to_value(std::string_view sv);

}

#endif
