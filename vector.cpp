#include "./vector.hpp"

namespace lib_math {
	static const t_vector4f    pi_vector_4f = {M_PI  , M_PI  , M_PI  ,  0.0f};
	static const t_vector4f   eps_vector_4f = {M_FEPS, M_FEPS, M_FEPS,  0.0f};
	static const t_vector4f   inf_vector_4f = {M_FINF, M_FINF, M_FINF,  0.0f};
	static const t_vector4f  zero_vector_4f = {  0.0f,   0.0f,   0.0f,  0.0f};
	static const t_vector4f  ones_vector_4f = {  1.0f,   1.0f,   1.0f,  0.0f};
	static const t_vector4f error_vector_4f = { -1.0f,  -1.0f,  -1.0f,  0.0f};
	static const t_vector4f     x_vector_4f = {  1.0f,   0.0f,   0.0f,  0.0f};
	static const t_vector4f     y_vector_4f = {  0.0f,   1.0f,   0.0f,  0.0f};
	static const t_vector4f     z_vector_4f = {  0.0f,   0.0f,   1.0f,  0.0f};
	static const t_vector4f     w_vector_4f = {  0.0f,   0.0f,   0.0f,  1.0f};
	static const t_vector4f    xz_vector_4f = x_vector_4f + z_vector_4f;
	static const t_vector4f    xy_vector_4f = x_vector_4f + y_vector_4f;
	static const t_vector4f    yz_vector_4f = y_vector_4f + z_vector_4f;
	static const t_vector4f   xyz_vector_4f = x_vector_4f + y_vector_4f + z_vector_4f;
	static const t_vector4f  xyzw_vector_4f = x_vector_4f + y_vector_4f + z_vector_4f + w_vector_4f;

	template<> const t_vector4f& t_vector4f::  eps_vector() { return   eps_vector_4f; }
	template<> const t_vector4f& t_vector4f::  inf_vector() { return   inf_vector_4f; }
	template<> const t_vector4f& t_vector4f:: zero_vector() { return  zero_vector_4f; }
	template<> const t_vector4f& t_vector4f:: ones_vector() { return  ones_vector_4f; }
	template<> const t_vector4f& t_vector4f::error_vector() { return error_vector_4f; }
	template<> const t_vector4f& t_vector4f::    x_vector() { return     x_vector_4f; }
	template<> const t_vector4f& t_vector4f::    y_vector() { return     y_vector_4f; }
	template<> const t_vector4f& t_vector4f::    z_vector() { return     z_vector_4f; }
	template<> const t_vector4f& t_vector4f::    w_vector() { return     w_vector_4f; }
	template<> const t_vector4f& t_vector4f::   xz_vector() { return    xz_vector_4f; }
	template<> const t_vector4f& t_vector4f::   xy_vector() { return    xy_vector_4f; }
	template<> const t_vector4f& t_vector4f::   yz_vector() { return    yz_vector_4f; }
	template<> const t_vector4f& t_vector4f::  xyz_vector() { return   xyz_vector_4f; }
	template<> const t_vector4f& t_vector4f:: xyzw_vector() { return  xyzw_vector_4f; }
};

