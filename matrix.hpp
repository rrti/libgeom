#ifndef LIBGEOM_MATRIX_HDR
#define LIBGEOM_MATRIX_HDR

#include <cmath>
#include <cassert>
#include <cstddef>
#include <algorithm>

#include "./math_defs.hpp"
#include "./tuple.hpp"

// forward-declarations
namespace lib_math {
	template<typename type> struct t_point;
	template<typename type> struct t_vector;

	template<typename type> struct t_matrix {
	typedef t_point<type> m_point_type;
	typedef t_vector<type> m_vector_type;
	public:
		t_matrix<type>(const type* values) { set_values(values); }
		t_matrix<type>(const t_matrix<type>& m) { *this = m; }
		t_matrix<type>(
			const m_vector_type& x_vec = m_vector_type::x_axis_vector(),
			const m_vector_type& y_vec = m_vector_type::y_axis_vector(),
			const m_vector_type& z_vec = m_vector_type::z_axis_vector(),
			const m_vector_type& t_vec = m_vector_type::w_axis_vector()
		) {
			set_x_vector(x_vec);
			set_y_vector(y_vec);
			set_z_vector(z_vec);
			set_t_vector(t_vec);
		}

		t_matrix<type>& operator = (const t_matrix<type>& m) {
			set_values(m.get_values()); return *this;
		}



		// matrix * point = point; 3x3 rotation + 1x3 translation
		// homogeneous w-coordinate of a point is 1 --> T included
		m_point_type operator * (const m_point_type& p) const {
			m_point_type tp;
			tp.x() = (m_values[0] * p.x()) + (m_values[4] * p.y()) + (m_values[ 8] * p.z()) + (m_values[12] * p.w());
			tp.y() = (m_values[1] * p.x()) + (m_values[5] * p.y()) + (m_values[ 9] * p.z()) + (m_values[13] * p.w());
			tp.z() = (m_values[2] * p.x()) + (m_values[6] * p.y()) + (m_values[10] * p.z()) + (m_values[14] * p.w());
			tp.w() = p.w();
			return tp;
		}

		// matrix * vector; interpret vector as point on unit-sphere
		// homogeneous w-coordinate of a vector is 0 --> T excluded
		//
		// only valid if matrix is orthonormal and vector has unit-length
		// in that case INVERSE(TRANSPOSE(R33)) == R33 and the vector can
		// be rotated directly (without requiring re-normalization due to
		// being translated as a point would be)
		//
		// aka "transform_vector"
		m_vector_type operator * (const m_vector_type& v) const {
			m_vector_type tv;
			tv.x() = (m_values[0] * v.x()) + (m_values[4] * v.y()) + (m_values[ 8] * v.z()); // equals n.inner_product(row[0])
			tv.y() = (m_values[1] * v.x()) + (m_values[5] * v.y()) + (m_values[ 9] * v.z()); // equals n.inner_product(row[1])
			tv.z() = (m_values[2] * v.x()) + (m_values[6] * v.y()) + (m_values[10] * v.z()); // equals n.inner_product(row[2])
			tv.w() = v.w();
			return tv;
		}



		// matrix * matrix = matrix
		t_matrix<type> operator * (const t_matrix<type>& m) const {
			t_matrix<type> r;

			// rows 0-2 with column 0 of m ( 0- 3)
			r[ 0] = (m_values[0] * m[ 0]) + (m_values[4] * m[ 1]) + (m_values[ 8] * m[ 2]) + (m_values[12] * m[ 3]);
			r[ 1] = (m_values[1] * m[ 0]) + (m_values[5] * m[ 1]) + (m_values[ 9] * m[ 2]) + (m_values[13] * m[ 3]);
			r[ 2] = (m_values[2] * m[ 0]) + (m_values[6] * m[ 1]) + (m_values[10] * m[ 2]) + (m_values[14] * m[ 3]);

			// rows 0-2 with column 1 of m ( 4- 7)
			r[ 4] = (m_values[0] * m[ 4]) + (m_values[4] * m[ 5]) + (m_values[ 8] * m[ 6]) + (m_values[12] * m[ 7]);
			r[ 5] = (m_values[1] * m[ 4]) + (m_values[5] * m[ 5]) + (m_values[ 9] * m[ 6]) + (m_values[13] * m[ 7]);
			r[ 6] = (m_values[2] * m[ 4]) + (m_values[6] * m[ 5]) + (m_values[10] * m[ 6]) + (m_values[14] * m[ 7]);

			// rows 0-2 with column 2 of m ( 8-11)
			r[ 8] = (m_values[0] * m[ 8]) + (m_values[4] * m[ 9]) + (m_values[ 8] * m[10]) + (m_values[12] * m[11]);
			r[ 9] = (m_values[1] * m[ 8]) + (m_values[5] * m[ 9]) + (m_values[ 9] * m[10]) + (m_values[13] * m[11]);
			r[10] = (m_values[2] * m[ 8]) + (m_values[6] * m[ 9]) + (m_values[10] * m[10]) + (m_values[14] * m[11]);

			// rows 0-2 with column 3 of m (12-15)
			r[12] = (m_values[0] * m[12]) + (m_values[4] * m[13]) + (m_values[ 8] * m[14]) + (m_values[12] * m[15]);
			r[13] = (m_values[1] * m[12]) + (m_values[5] * m[13]) + (m_values[ 9] * m[14]) + (m_values[13] * m[15]);
			r[14] = (m_values[2] * m[12]) + (m_values[6] * m[13]) + (m_values[10] * m[14]) + (m_values[14] * m[15]);

			return r;
		}

		t_matrix<type> operator * (const type values[LIBGEOM_MATRIX_SIZE]) const {
			return ((*this) * t_matrix<type>(&values[0]));
		}


		// matrix + matrix = matrix
		t_matrix<type> operator + (const t_matrix<type>& m) const {
			t_matrix<type> m_r;

			for (unsigned int n = 0; n < LIBGEOM_MATRIX_SIZE; n += 4) {
				m_r[n + 0] = m_values[n + 0] + m[n + 0];
				m_r[n + 1] = m_values[n + 1] + m[n + 1];
				m_r[n + 2] = m_values[n + 2] + m[n + 2];
				m_r[n + 3] = m_values[n + 3] + m[n + 3];
			}

			return m_r;
		}

		/*
		t_matrix<type>  operator +  (const t_matrix<type>& m) const { t_matrix<type> mr(false); mr.add_values(m_values); mr.add_values(m.get_values()); return    mr; }
		t_matrix<type>  operator *  (const t_matrix<type>& m) const { t_matrix<type> mr(false); mr.add_values(m_values); mr.mul_values(m.get_values()); return    mr; }
		t_matrix<type>& operator += (const t_matrix<type>& m) const {                                                       add_values(m.get_values()); return *this; }
		t_matrix<type>& operator *= (const t_matrix<type>& m) const {                                                       mul_values(m.get_values()); return *this; }
		*/

		// matrix {*,+} matrix = matrix
		// slower than direct in-place updates but avoids duplication
		t_matrix<type>& operator *= (const t_matrix<type>& m) { return ((*this) = (*this) * m); }
		t_matrix<type>& operator += (const t_matrix<type>& m) { return ((*this) = (*this) + m); }


		// matrix * scalar = matrix
		t_matrix<type> operator * (const type s) const {
			t_matrix<type> m_r;

			for (unsigned int n = 0; n < LIBGEOM_MATRIX_SIZE; n += 4) {
				m_r[n + 0] = m_values[n + 0] * s;
				m_r[n + 1] = m_values[n + 1] * s;
				m_r[n + 2] = m_values[n + 2] * s;
				m_r[n + 3] = m_values[n + 3] * s;
			}

			return m_r;
		}

		t_matrix<type>& operator *= (const type s) { return ((*this) = (*this) * s); }


		type  operator [] (unsigned int idx) const { return m_values[idx]; }
		type& operator [] (unsigned int idx)       { return m_values[idx]; }

		// <idx> should be one of {0,4,8,12}
		const type* get_raw_vector(unsigned int idx) const { return &m_values[idx]; }
			  type* get_raw_vector(unsigned int idx)       { return &m_values[idx]; }
		const type* get_values() const { return &m_values[0]; }
		      type* get_values()       { return &m_values[0]; }


		// column accessors (NOTE: the t-vector is really a point!)
		m_vector_type get_x_vector() const { m_vector_type v; v.x() = m_values[ 0]; v.y() = m_values[ 1]; v.z() = m_values[ 2]; v.w() = m_values[ 3]; return v; }
		m_vector_type get_y_vector() const { m_vector_type v; v.x() = m_values[ 4]; v.y() = m_values[ 5]; v.z() = m_values[ 6]; v.w() = m_values[ 7]; return v; }
		m_vector_type get_z_vector() const { m_vector_type v; v.x() = m_values[ 8]; v.y() = m_values[ 9]; v.z() = m_values[10]; v.w() = m_values[11]; return v; }
		m_vector_type get_t_vector() const { m_vector_type v; v.x() = m_values[12]; v.y() = m_values[13]; v.z() = m_values[14]; v.w() = m_values[15]; return v; }

		// random accessors (<idx> should be one of {0,1,2,3})
		m_vector_type get_row_vector(unsigned int idx) const { return (m_vector_type(m_values[idx], m_values[idx + 4], m_values[idx + 8], m_values[idx + 12])); }
		m_vector_type get_col_vector(unsigned int idx) const { return (m_vector_type(m_values[(idx * 4) + 0], m_values[(idx * 4) + 1], m_values[(idx * 4) + 2], m_values[(idx * 4) + 3])); }


		t_matrix<type> translate(const m_vector_type& v) const {
			t_matrix<type> r = *this;
			r[12] += ((v.x() * r[0]) + (v.y() * r[4]) + (v.z() * r[ 8])); // equals tv.inner_product(rows[0])
			r[13] += ((v.x() * r[1]) + (v.y() * r[5]) + (v.z() * r[ 9])); // equals tv.inner_product(rows[1])
			r[14] += ((v.x() * r[2]) + (v.y() * r[6]) + (v.z() * r[10])); // equals tv.inner_product(rows[2])
			r[15] += ((v.x() * r[3]) + (v.y() * r[7]) + (v.z() * r[11])); // equals tv.inner_product(rows[3])
			return r;
		}
		t_matrix<type> scale(const m_vector_type& sv) const {
			t_matrix r = *this;

			r[ 0] *= sv.x(); r[ 4] *= sv.y();
			r[ 1] *= sv.x(); r[ 5] *= sv.y();
			r[ 2] *= sv.x(); r[ 6] *= sv.y();

			r[ 8] *= sv.z(); // r[12] *= sv.w();
			r[ 9] *= sv.z(); // r[13] *= sv.w();
			r[10] *= sv.z(); // r[14] *= sv.w();

			return r;
		}

		t_matrix<type> transpose_rotation() const {
			t_matrix<type> r = *this;
			std::swap(r[1], r[4]);
			std::swap(r[2], r[8]);
			std::swap(r[6], r[9]);
			return r;
		}
		t_matrix<type> transpose_translation() const {
			t_matrix<type> r = *this;
			std::swap(r[ 3], r[12]);
			std::swap(r[ 7], r[13]);
			std::swap(r[11], r[14]);
			return r;
		}
		t_matrix<type> transpose() const {
			t_matrix<type> r = *this;
			r.transpose_rotation_ref();
			r.transpose_translation_ref();
			return r;
		}

		// "affine" assumes matrix only does translation and rotation
		t_matrix<type> invert_affine() const {
			t_matrix<type> r = *this;

			// transpose the rotation
			r.transpose_rotation_ref();

			// get the inverse translation
			const m_vector_type t_pre = -r.get_t_vector();
			const m_vector_type t_inv = r * t_pre;

			// do the positional inversion
			r.set_t_vector(t_inv);
			return r;
		}
		// generalized inverse for non-orthonormal 4x4 matrices
		// A^-1 = (1 / det(A)) (C^T)_{ij} = (1 / det(A)) C_{ji}
		// where C is the matrix of cofactors
		//
		t_matrix<type> invert_general(const type eps = t_tuple<type>::eps_scalar()) const;

		t_matrix& translate_ref(const m_vector_type& tv) { return ((*this) = translate(tv)); }
		t_matrix& scale_ref(const m_vector_type& sv) { return ((*this) = scale(sv)); }

		t_matrix<type>& transpose_ref() { return ((*this) = transpose()); }
		t_matrix<type>& transpose_rotation_ref() { return ((*this) = transpose_rotation()); }
		t_matrix<type>& transpose_translation_ref() { return ((*this) = transpose_translation()); }
		t_matrix<type>& invert_affine_ref() { return ((*this) = invert_affine()); }


		// NOTE:
		//   these are intrinsic: y-coordinate of Z-axis
		//   drops after rotate_x and again after rotate_y
		//
		// rotate_x on an identity-vector (v) gives v' = <x' = x,  y' = ca*y - sa*z,  z' = sa*y + ca*z>  (CCW)
		// rotate_x on an identity-matrix (R) gives R' = {X' = X,  Y' = <0, ca, -sa>,  Z' = <0, sa, ca>}  (CW?)
		//
		// assume v = Z = <0,0,1>, then v' = Z' = <0, ca*0 - sa*1, sa*0 + ca*1> = <0, -sa, ca>
		// which is the OPPOSITE (clockwise instead of counter-clockwise) of the Z-column of R'
		// (R * v != v')
		//
		// rotate in yz-plane
		t_matrix<type>& rotate_x_int(const type angle_radians) {
			const type ca = std::cos(angle_radians);
			const type sa = std::sin(angle_radians);

			t_matrix<type> r;
			r[ 5] = +ca; r[ 9] = +sa; // y.y(), z.y()
			r[ 6] = -sa; r[10] = +ca; // y.z(), z.z()

			// make the rotation counter-clockwise
			r.transpose_rotation_ref();

			return ((*this) *= r);
		}
		// rotate in xz-plane
		t_matrix<type>& rotate_y_int(const type angle_radians) {
			const type ca = std::cos(angle_radians);
			const type sa = std::sin(angle_radians);

			t_matrix<type> r;
			r[ 0] = +ca; r[ 8] = -sa; // x.x(), z.x()
			r[ 2] = +sa; r[10] = +ca; // x.z(), z.z()

			// make the rotation counter-clockwise
			r.transpose_rotation_ref();

			return ((*this) *= r);
		}
		// rotate in xy-plane
		t_matrix<type>& rotate_z_int(const type angle_radians) {
			const type ca = std::cos(angle_radians);
			const type sa = std::sin(angle_radians);

			t_matrix<type> r;
			r[0] = +ca; r[4] = +sa; // x.x(), y.x()
			r[1] = -sa; r[5] = +ca; // x.y(), y.y()

			// make the rotation counter-clockwise
			r.transpose_rotation_ref();

			return ((*this) *= r);
		}

		// rotation about arbitrary axis, ie. in the
		// plane crossing the origin whose normal is
		// <axis.x, axis.y, axis.z>
		// any such rotation can be decomposed into
		// rotations about the three component axes
		t_matrix<type>& rotate_ref(const m_vector_type& rot_axis, const type angle_radians) {
			const m_vector_type vx = (rot_axis.outer_product(t_vector<type>::y_axis_vector())).normalize();
			const m_vector_type vy = (rot_axis.outer_product(                           vx  )).normalize();

			const t_matrix<type> m_fwd = t_matrix<type>(vx, vy, rot_axis);
			const t_matrix<type> m_inv = m_fwd.invert_affine();
			const t_matrix<type> m_rot = t_matrix<type>().rotate_z_int(angle_radians);

			return ((*this) = m_fwd * m_rot * m_inv);
		}
		t_matrix<type>& rotate_ref_tmp(const m_vector_type& rot_axis, const type angle_radians) {
			return (rotate_ref(rot_axis, angle_radians));

			/*
			// does not preserve orthonormality
			// for non-principal rotation axes!
			const type ca = std::cos(angle_radians);
			const type sa = std::sin(angle_radians);

			const m_vector_type x_axis = get_x_vector();
			const m_vector_type y_axis = get_y_vector();
			const m_vector_type z_axis = get_z_vector();

			// project rotation axis onto each principal axis
			const m_vector_type fwd_axis_x = rot_axis * x_axis.inner_product(rot_axis);
			const m_vector_type fwd_axis_y = rot_axis * y_axis.inner_product(rot_axis);
			const m_vector_type fwd_axis_z = rot_axis * z_axis.inner_product(rot_axis);

			// NOTE: rgt and up define the rotational plane
			const m_vector_type rgt_axis_x = x_axis - fwd_axis_x;
			const m_vector_type rgt_axis_y = y_axis - fwd_axis_y;
			const m_vector_type rgt_axis_z = z_axis - fwd_axis_z;

			const m_vector_type up_axis_x = (rot_axis.outer_product(rgt_axis_x)).normalize();
			const m_vector_type up_axis_y = (rot_axis.outer_product(rgt_axis_y)).normalize();
			const m_vector_type up_axis_z = (rot_axis.outer_product(rgt_axis_z)).normalize();

			set_x_vector((fwd_axis_x + (rgt_axis_x * ca + up_axis_x * sa)));
			set_y_vector((fwd_axis_y + (rgt_axis_y * ca + up_axis_y * sa)));
			set_z_vector((fwd_axis_z + (rgt_axis_z * ca + up_axis_z * sa)));

			return *this;
			*/
		}



		void set_x_vector(const m_vector_type& v) { m_values[ 0] = v.x(); m_values[ 1] = v.y(); m_values[ 2] = v.z(); m_values[ 3] = type(0); } // v.w();
		void set_y_vector(const m_vector_type& v) { m_values[ 4] = v.x(); m_values[ 5] = v.y(); m_values[ 6] = v.z(); m_values[ 7] = type(0); } // v.w();
		void set_z_vector(const m_vector_type& v) { m_values[ 8] = v.x(); m_values[ 9] = v.y(); m_values[10] = v.z(); m_values[11] = type(0); } // v.w();
		void set_t_vector(const m_vector_type& v) { m_values[12] = v.x(); m_values[13] = v.y(); m_values[14] = v.z(); m_values[15] = type(1); } // v.w();

		void set_row_vector(unsigned int idx, const m_vector_type& v) { m_values[(idx    ) + 0] = v.x(); m_values[(idx    ) + 4] = v.y(); m_values[(idx    ) + 8] = v.z(); m_values[(idx    ) + 12] = v.w(); }
		void set_col_vector(unsigned int idx, const m_vector_type& v) { m_values[(idx * 4) + 0] = v.x(); m_values[(idx * 4) + 1] = v.y(); m_values[(idx * 4) + 2] = v.z(); m_values[(idx * 4) +  3] = v.w(); }


		t_matrix<type>& set_xz_vectors_from_y_vector(const m_vector_type& vy_new) {
			set_y_vector(vy_new);

			// if y=<0,1,0> and z=<0,0,1>, sets x=<1,0,0> and z=<0,0,1>
			const m_vector_type& vz_old = get_z_vector();
			const m_vector_type  vx_new = (vy_new.outer_product(vz_old)).normalize();
			const m_vector_type  vz_new = (vx_new.outer_product(vy_new)).normalize();

			// avoid degeneracy if abs(inner_product(vy_new, vz_old)) ~= 1
			assert(std::fabs(vy_new.inner_product(vz_old)) < (type(1) - t_tuple<type>::eps_scalar()));

			set_x_vector(vx_new);
			set_z_vector(vz_new);
			return *this;
		}


		t_matrix<type>& add_values(const float values[LIBGEOM_MATRIX_SIZE]) {
			assert(values != NULL);
			for (unsigned int i = 0; i < LIBGEOM_MATRIX_SIZE; i += 4) {
				m_values[i + 0] += values[i + 0];
				m_values[i + 1] += values[i + 1];
				m_values[i + 2] += values[i + 2];
				m_values[i + 3] += values[i + 3];
			}
			return *this;
		}

		t_matrix<type>& set_values(const float values[LIBGEOM_MATRIX_SIZE]) {
			assert(values != NULL);
			for (unsigned int i = 0; i < LIBGEOM_MATRIX_SIZE; i += 4) {
				m_values[i + 0] = values[i + 0];
				m_values[i + 1] = values[i + 1];
				m_values[i + 2] = values[i + 2];
				m_values[i + 3] = values[i + 3];
			}
			return *this;
		}

		void set_unit_values();
		void print();

		unsigned int is_identity(const type eps = t_tuple<type>::eps_scalar()) const;
		unsigned int is_orthonormal(const type eps = t_tuple<type>::eps_scalar()) const;

		//	matrix_outer(v, w) =
		//	  ( 0  -v3  v2)   (w1)
		//	  ( v3  0  -v1) * (w2)
		//	  (-v2  v1  0 )   (w3)
		//
		// NOTE: useful to generalize this to higher dimensions
		static m_vector_type outer_product(const m_vector_type& v, const m_vector_type& w) {
			t_matrix<type> m;
			m.set_x_vector(m_vector_type(type(0),   v.z(),  -v.y()));
			m.set_y_vector(m_vector_type( -v.z(), type(0),   v.x()));
			m.set_z_vector(m_vector_type(  v.y(),  -v.x(), type(0)));
			return (m * w);
		}

		static const type* unit_values();
		static const type* null_values();

	private:
		// stores elements in column-major order such that
		// m_data[0...3] represents first column (x), etc.
		type m_values[LIBGEOM_MATRIX_SIZE];
	};

	typedef t_matrix< float> t_matrix44f;
	typedef t_matrix<double> t_matrix44d;
};

#endif

