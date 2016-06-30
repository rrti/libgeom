#include <cstdio>

#include "./matrix.hpp"
#include "./point.hpp"
#include "./vector.hpp"

namespace lib_math {
	static const float null_values_44f[4 * 4] = {
		0.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 0.0f,
	};
	static const float unit_values_44f[4 * 4] = {
		1.0f, 0.0f, 0.0f, 0.0f, // 1st column
		0.0f, 1.0f, 0.0f, 0.0f, // 2nd column
		0.0f, 0.0f, 1.0f, 0.0f, // 3rd column
		0.0f, 0.0f, 0.0f, 1.0f, // 4th column
	};


	#if 0
	template<> t_matrix44f t_matrix44f::invert_general(const float eps) const {
		t_matrix44f mat;
		t_matrix44f& cofac = mat;

		for (unsigned int i = 0; i < LIBGEOM_MATRIX_SIZE; i++) {
			for (unsigned int j = 0; j < LIBGEOM_MATRIX_SIZE; j++) {
				cofac.m_values[i][j] = CalculateCofactor(m_values, i, j);
			}
		}

		const float det =
			(m_values[0][0] * cofac.m_values[0][0]) +
			(m_values[0][1] * cofac.m_values[0][1]) +
			(m_values[0][2] * cofac.m_values[0][2]) +
			(m_values[0][3] * cofac.m_values[0][3]);

		if (det < eps) {
			// <this> is a singular matrix
			return mat;
		}

		mat *= (1.0f / det);
		mat.transpose_ref();
		return mat;
	}
	#endif

	template<> void t_matrix44f::set_unit_values() {
		#if 0
		for (unsigned int n = 0; n < 4 * 4; n += 4) {
			m_values[n + 0] = float((n + 0) ==  0);
			m_values[n + 1] = float((n + 1) ==  5);
			m_values[n + 2] = float((n + 2) == 10);
			m_values[n + 3] = float((n + 3) == 15);
		}
		#else
		for (unsigned int n = 0; n < 4 * 4; n += 4) {
			m_values[n + 0] = float(n ==  0);
			m_values[n + 1] = float(n ==  4);
			m_values[n + 2] = float(n ==  8);
			m_values[n + 3] = float(n == 12);
		}
		#endif
	}
	template<> void t_matrix44f::print() {
		printf("\n X      Y      Z      T\n");
		for (unsigned int n = 0; n < LIBGEOM_MATRIX_SIZE; n += 4) { printf(" %.3f ", m_values[n]); } printf("\n"); // 1st row
		for (unsigned int n = 1; n < LIBGEOM_MATRIX_SIZE; n += 4) { printf(" %.3f ", m_values[n]); } printf("\n"); // 2nd row
		for (unsigned int n = 2; n < LIBGEOM_MATRIX_SIZE; n += 4) { printf(" %.3f ", m_values[n]); } printf("\n"); // 3rd row
		for (unsigned int n = 3; n < LIBGEOM_MATRIX_SIZE; n += 4) { printf(" %.3f ", m_values[n]); } printf("\n"); // 4th row
		printf("\n");
	}


	template<> unsigned int t_matrix44f::is_orthonormal(const float eps) const {
		const m_vector_type& xv = get_x_vector();
		const m_vector_type& yv = get_y_vector();
		const m_vector_type& zv = get_z_vector();

		// test orthogonality
		if (std::fabs(xv.inner_product(yv)) >= eps) { return 1; }
		if (std::fabs(yv.inner_product(zv)) >= eps) { return 2; }
		if (std::fabs(xv.inner_product(zv)) >= eps) { return 3; }
		// test normality
		if (std::fabs(1.0f - xv.sq_magnit()) >= eps) { return 4; }
		if (std::fabs(1.0f - yv.sq_magnit()) >= eps) { return 5; }
		if (std::fabs(1.0f - zv.sq_magnit()) >= eps) { return 6; }

		return 0;
	}

	template<> unsigned int t_matrix44f::is_identity(const float eps) const {
		const m_vector_type& xv = get_x_vector();
		const m_vector_type& yv = get_y_vector();
		const m_vector_type& zv = get_z_vector();
		const m_vector_type& tv = get_t_vector();

		if (!xv.equals(m_vector_type::x_axis_vector(), m_vector_type::x_axis_vector() * eps)) return 1;
		if (!yv.equals(m_vector_type::y_axis_vector(), m_vector_type::y_axis_vector() * eps)) return 2;
		if (!zv.equals(m_vector_type::z_axis_vector(), m_vector_type::z_axis_vector() * eps)) return 3;
		if (!tv.equals(m_vector_type::w_axis_vector(), m_vector_type::w_axis_vector() * eps)) return 4;

		return 0;
	}


	template<> const float* t_matrix44f::unit_values() { return &unit_values_44f[0]; }
	template<> const float* t_matrix44f::null_values() { return &null_values_44f[0]; }
};

