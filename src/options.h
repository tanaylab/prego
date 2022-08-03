#ifndef util_options_h
#define util_options_h 1

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>



class options {

private:

	std::map<std::string, std::string> scalars;

	std::map<std::string, std::vector<int> *> int_arrays;
	std::map<std::string, std::vector<float> *> float_arrays;
	std::map<std::string, std::vector<std::string> *> string_arrays;

	std::vector<int> m_empty_int_vec;
	std::vector<float> m_empty_float_vec;
	std::vector<std::string> m_empty_str_vec;
//auxilaries

	std::string cur_scope;
	std::string cur_name;

public:

	bool m_logit;

	options() {}

	options(std::istream &in) { read(in); }

	void load_defaults(const char *vals[], int size);
	void parse_argv(int &argc, char *argv[]);

	void read(std::istream &in, bool logit = false);

//Access methods

	int get_int(const std::string &sc, const std::string &nm, bool must = true, int def = 0);
	int get_g_int(const std::string &nm, bool must = true, int def = 0) {
		return(get_int("", nm, must, def));
	}
	float get_float(const std::string &sc, const std::string &nm, bool must = true, float def = 0);
	float get_g_float(const std::string &nm, bool must = true, int def = 0) {
		return(get_float("", nm, must, def));
	}
	const std::string &get_str(const std::string &sc, const std::string &nm, bool must = true, const std::string *def = 0);

	const std::string &get_g_str(const std::string &nm, bool must = true, const std::string *def = 0) {
		return(get_str("", nm, must, def));
	}

	const std::vector<int> &get_ints(const std::string &sc, const std::string &nm,
							bool must = true);
	const std::vector<float> &get_floats(const std::string &sc, const std::string &nm,
							bool must = true);
	const std::vector<std::string> &get_strs(const std::string &sc, const std::string &nm,
							bool must = true);

	bool have_scalar(const std::string &sc, const std::string &nm) {
		return(scalars.count(sc+"::"+nm));
	}
	bool have_g_scalar(const std::string &nm) {
		return(scalars.count("::"+nm));
	}

	void set_scalar(const std::string &sc, const std::string &nm, const std::string &val) {
		scalars[sc+"::"+nm] = val;
	}
	void set_g_scalar(const std::string &nm, const std::string &val) {
		scalars["::"+nm] = val;
	}

private :

//Parsing methods

	bool is_comment_line(std::string &line);
	bool is_scope_line(std::string &line);
	bool get_option_name(std::string &line);
	bool is_option_line(std::string &line);
	bool is_array_line(std::string &line, std::istream &in);

};

#endif //util_options_h
