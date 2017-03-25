#include <cstdio>

#include "./tuple.hpp"

namespace lib_math {
	static const t_tuple4f zero_tuple_4f = {  0.0f,   0.0f,   0.0f,   0.0f};
	static const t_tuple4f ones_tuple_4f = {  1.0f,   1.0f,   1.0f,   1.0f};
	static const t_tuple4f  eps_tuple_4f = {M_FEPS, M_FEPS, M_FEPS, M_FEPS};

	static const t_tuple4d zero_tuple_4d = {   0.0,    0.0,    0.0,    0.0};
	static const t_tuple4d ones_tuple_4d = {   1.0,    1.0,    1.0,    1.0};
	static const t_tuple4d  eps_tuple_4d = {M_DEPS, M_DEPS, M_DEPS, M_DEPS};

	static const t_tuple4i zero_tuple_4i = {0,  0,  0,  0};
	static const t_tuple4i ones_tuple_4i = {1,  1,  1,  1};
	static const t_tuple4i  eps_tuple_4i = {0,  0,  0,  0};

	static const t_tuple4ui zero_tuple_4ui = {0,  0,  0,  0};
	static const t_tuple4ui ones_tuple_4ui = {1,  1,  1,  1};
	static const t_tuple4ui  eps_tuple_4ui = {0,  0,  0,  0};

	template<> const t_tuple4d& t_tuple4d::zero_tuple() { return zero_tuple_4d; }
	template<> const t_tuple4d& t_tuple4d::ones_tuple() { return ones_tuple_4d; }
	template<> const t_tuple4d& t_tuple4d:: eps_tuple() { return  eps_tuple_4d; }

	template<> const t_tuple4f& t_tuple4f::zero_tuple() { return zero_tuple_4f; }
	template<> const t_tuple4f& t_tuple4f::ones_tuple() { return ones_tuple_4f; }
	template<> const t_tuple4f& t_tuple4f:: eps_tuple() { return  eps_tuple_4f; }

	template<> const t_tuple4ui& t_tuple4ui::zero_tuple() { return zero_tuple_4ui; }
	template<> const t_tuple4ui& t_tuple4ui::ones_tuple() { return ones_tuple_4ui; }
	template<> const t_tuple4ui& t_tuple4ui:: eps_tuple() { return  eps_tuple_4ui; }

	template<> void t_tuple4d::print() const { std::printf("<x=%.3g, y=%.3g, z=%.3g, w=%.3g>\n", x(), y(), z(), w()); }
	template<> void t_tuple4f::print() const { std::printf("<x=%.3f, y=%.3f, z=%.3f, w=%.3f>\n", x(), y(), z(), w()); }
	// template<> void t_tuple4d::print() const { std::printf("<x=%.3g, y=%.3g, z=%.3g, w=%.3g h=%u>\n", x(), y(), z(), w(), hash()); }
	// template<> void t_tuple4f::print() const { std::printf("<x=%.3f, y=%.3f, z=%.3f, w=%.3f h=%u>\n", x(), y(), z(), w(), hash()); }
};

