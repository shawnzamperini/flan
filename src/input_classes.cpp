#include <string>
#include <iomanip>
#include <charconv>
#include <iostream>

#include "input_classes.h"

namespace Input
{
	// Base class to derive different input option types (int, double. string, ...)
	// from.  
	OptionBase::OptionBase(std::string name)
		: m_name {name}
	{}

	std::string_view OptionBase::get_name() const {return m_name;}
	
	// Integer option class
	OptionInt::OptionInt(std::string name, int value)
		: OptionBase(name)
		, m_value {value}
	{}

	int OptionInt::get_value() const {return m_value;}
	void OptionInt::put_value(int value) {m_value = value;}
	void OptionInt::put_value(std::string value) {m_value = std::stoi(value);}
	void OptionInt::put_value(std::string_view sv)
	{
		m_value = sv_to_value<int>(sv);
	}
	
	// Double option class
	OptionDbl::OptionDbl(std::string name, double value)
		: OptionBase(name)
		, m_value {value}
	{}

	double OptionDbl::get_value() const {return m_value;}
	void OptionDbl::put_value(double value) {m_value = value;}
	void OptionDbl::put_value(std::string value) {m_value = std::stod(value);}
	void OptionDbl::put_value(std::string_view sv)
	{
		m_value = sv_to_value<double>(sv);
	}

	// String input option class.
	OptionStr::OptionStr(std::string name, std::string value)
		: OptionBase(name)
		, m_value {value}
	{}

	// Getter and putter
	std::string_view OptionStr::get_value() const {return m_value;}
	void OptionStr::put_value(std::string_view value) {m_value = value;}

	// Overload << so we can easily print an option out in a nice format
	std::ostream& operator<<(std::ostream& out, const OptionStr opt)
	{
		out << std::setw(15) << opt.get_name() << "  " << opt.get_value();
		return out;
	}

	std::ostream& operator<<(std::ostream& out, const OptionInt opt)
	{
		out << std::setw(15) << opt.get_name() << "  " << opt.get_value();
		return out;
	}

	std::ostream& operator<<(std::ostream& out, const OptionDbl opt)
	{
		out << std::setw(15) << opt.get_name() << "  " << opt.get_value();
		return out;
	}

	// Function needed to convert a string_view to either an int or double. Used
	// in OptionInt and OptionDouble. 
	template <typename T>
	T sv_to_value(std::string_view sv)
	{
		T value {};
		auto result = std::from_chars(sv.data(), sv.data() + sv.size(), value);
		if (result.ec == std::errc::invalid_argument)
		{
			std::cerr << "Could not convert to int/double: " << sv << '\n';
		}
		return value;
	}

}
