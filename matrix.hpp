#ifndef LIBGEOM_MATRIX_HDR
#define LIBGEOM_MATRIX_HDR

#include <cmath>
#include <cassert>
#include <cstddef>

#include "./math_defs.hpp"
#include "./quaternion.hpp"
#include "./tuple.hpp"

#define MATRIX_SIZE 16

// forward-declarations
namespace lib_math {
	template<typename type> struct t_point4t;
	template<typename type> struct t_vector4t;

	template<typename type> struct t_matrix44t {
	typedef t_point4t<type> m_pos_type;
	typedef t_vector4t<type> m_vec_type;
	typedef t_matrix44t<type> m_mat_type;
	public:
		t_matrix44t<type>(const type* m) { set_raw_matrix(m); }
		t_matrix44t<type>(const m_mat_type& m) { *this = m; }
		t_matrix44t<type>(
			const m_vec_type& x_vec = m_vec_type::x_vector(),
			const m_vec_type& y_vec = m_vec_type::y_vector(),
			const m_vec_type& z_vec = m_vec_type::z_vector(),
			const m_vec_type& t_vec = m_vec_type::w_vector()
		) {
			set_x_vec(x_vec);
			set_y_vec(y_vec);
			set_z_vec(z_vec);
			set_t_vec(t_vec);
		}

		m_mat_type& operator = (const m_mat_type& m) {
			set_raw_matrix(m.xyzt()); return *this;
		}


		bool operator == (const m_mat_type& m) const {
			unsigned int r = 0;

			for (unsigned int n = 0; n < MATRIX_SIZE; n += 4) {
				r += fp_eq(m_xyzt[n + 0], m[n + 0], M_FEPS);
				r += fp_eq(m_xyzt[n + 1], m[n + 1], M_FEPS);
				r += fp_eq(m_xyzt[n + 2], m[n + 2], M_FEPS);
				r += fp_eq(m_xyzt[n + 3], m[n + 3], M_FEPS);
			}

			return (r == MATRIX_SIZE);
		}

		bool operator != (const m_mat_type& m) const {
			return (!((*this) == m));
		}


		// matrix * point = point; 3x3 rotation + 1x3 translation
		// homogeneous w-coordinate of a point is 1 --> T included
		m_pos_type operator * (const m_pos_type& p) const {
			m_pos_type tp;
			tp.x() = (m_xyzt[0] * p.x()) + (m_xyzt[4] * p.y()) + (m_xyzt[ 8] * p.z()) + (m_xyzt[12] * p.w());
			tp.y() = (m_xyzt[1] * p.x()) + (m_xyzt[5] * p.y()) + (m_xyzt[ 9] * p.z()) + (m_xyzt[13] * p.w());
			tp.z() = (m_xyzt[2] * p.x()) + (m_xyzt[6] * p.y()) + (m_xyzt[10] * p.z()) + (m_xyzt[14] * p.w());
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
		m_vec_type operator * (const m_vec_type& v) const {
			m_vec_type tv;
			tv.x() = (m_xyzt[0] * v.x()) + (m_xyzt[4] * v.y()) + (m_xyzt[ 8] * v.z()); // same as n.inner(row[0])
			tv.y() = (m_xyzt[1] * v.x()) + (m_xyzt[5] * v.y()) + (m_xyzt[ 9] * v.z()); // same as n.inner(row[1])
			tv.z() = (m_xyzt[2] * v.x()) + (m_xyzt[6] * v.y()) + (m_xyzt[10] * v.z()); // same as n.inner(row[2])
			tv.w() = v.w();
			return tv;
		}



		// matrix * matrix = matrix
		m_mat_type operator * (const m_mat_type& m) const {
			m_mat_type r;

			// rows 0-2 with column 0 of m ( 0- 3)
			r[ 0] = (m_xyzt[0] * m[ 0]) + (m_xyzt[4] * m[ 1]) + (m_xyzt[ 8] * m[ 2]) + (m_xyzt[12] * m[ 3]);
			r[ 1] = (m_xyzt[1] * m[ 0]) + (m_xyzt[5] * m[ 1]) + (m_xyzt[ 9] * m[ 2]) + (m_xyzt[13] * m[ 3]);
			r[ 2] = (m_xyzt[2] * m[ 0]) + (m_xyzt[6] * m[ 1]) + (m_xyzt[10] * m[ 2]) + (m_xyzt[14] * m[ 3]);

			// rows 0-2 with column 1 of m ( 4- 7)
			r[ 4] = (m_xyzt[0] * m[ 4]) + (m_xyzt[4] * m[ 5]) + (m_xyzt[ 8] * m[ 6]) + (m_xyzt[12] * m[ 7]);
			r[ 5] = (m_xyzt[1] * m[ 4]) + (m_xyzt[5] * m[ 5]) + (m_xyzt[ 9] * m[ 6]) + (m_xyzt[13] * m[ 7]);
			r[ 6] = (m_xyzt[2] * m[ 4]) + (m_xyzt[6] * m[ 5]) + (m_xyzt[10] * m[ 6]) + (m_xyzt[14] * m[ 7]);

			// rows 0-2 with column 2 of m ( 8-11)
			r[ 8] = (m_xyzt[0] * m[ 8]) + (m_xyzt[4] * m[ 9]) + (m_xyzt[ 8] * m[10]) + (m_xyzt[12] * m[11]);
			r[ 9] = (m_xyzt[1] * m[ 8]) + (m_xyzt[5] * m[ 9]) + (m_xyzt[ 9] * m[10]) + (m_xyzt[13] * m[11]);
			r[10] = (m_xyzt[2] * m[ 8]) + (m_xyzt[6] * m[ 9]) + (m_xyzt[10] * m[10]) + (m_xyzt[14] * m[11]);

			// rows 0-2 with column 3 of m (12-15)
			r[12] = (m_xyzt[0] * m[12]) + (m_xyzt[4] * m[13]) + (m_xyzt[ 8] * m[14]) + (m_xyzt[12] * m[15]);
			r[13] = (m_xyzt[1] * m[12]) + (m_xyzt[5] * m[13]) + (m_xyzt[ 9] * m[14]) + (m_xyzt[13] * m[15]);
			r[14] = (m_xyzt[2] * m[12]) + (m_xyzt[6] * m[13]) + (m_xyzt[10] * m[14]) + (m_xyzt[14] * m[15]);

			return r;
		}

		m_mat_type operator * (const type* m) const {
			return ((*this) * m_mat_type(&m[0]));
		}


		// matrix + matrix = matrix
		m_mat_type operator + (const m_mat_type& m) const {
			m_mat_type m_r;

			for (unsigned int n = 0; n < MATRIX_SIZE; n += 4) {
				m_r[n + 0] = m_xyzt[n + 0] + m[n + 0];
				m_r[n + 1] = m_xyzt[n + 1] + m[n + 1];
				m_r[n + 2] = m_xyzt[n + 2] + m[n + 2];
				m_r[n + 3] = m_xyzt[n + 3] + m[n + 3];
			}

			return m_r;
		}
		// matrix - matrix = matrix
		m_mat_type operator - (const m_mat_type& m) const {
			return ((*this) + (m * type(-1)));
		}

		#if 0
		m_mat_type  operator +  (const m_mat_type& m) const { m_mat_type mr(false); mr.add_raw_matrix(m_xyzt); mr.add_raw_matrix(m.xyzt()); return    mr; }
		m_mat_type  operator *  (const m_mat_type& m) const { m_mat_type mr(false); mr.add_raw_matrix(m_xyzt); mr.mul_raw_matrix(m.xyzt()); return    mr; }
		m_mat_type& operator += (const m_mat_type& m) const {                                                     add_raw_matrix(m.xyzt()); return *this; }
		m_mat_type& operator *= (const m_mat_type& m) const {                                                     mul_raw_matrix(m.xyzt()); return *this; }
		#endif

		// matrix {*,+} matrix = matrix
		// slower than direct in-place updates but avoids duplication
		m_mat_type& operator *= (const m_mat_type& m) { return ((*this) = (*this) * m); }
		m_mat_type& operator += (const m_mat_type& m) { return ((*this) = (*this) + m); }


		// matrix * scalar = matrix
		m_mat_type operator * (const type s) const {
			m_mat_type m_r;

			for (unsigned int n = 0; n < MATRIX_SIZE; n += 4) {
				m_r[n + 0] = m_xyzt[n + 0] * s;
				m_r[n + 1] = m_xyzt[n + 1] * s;
				m_r[n + 2] = m_xyzt[n + 2] * s;
				m_r[n + 3] = m_xyzt[n + 3] * s;
			}

			return m_r;
		}

		m_mat_type& operator *= (const type s) { return ((*this) = (*this) * s); }


		type  operator [] (unsigned int idx) const { return m_xyzt[idx]; }
		type& operator [] (unsigned int idx)       { return m_xyzt[idx]; }

		const type* xyzt(unsigned int idx = 0) const { return &m_xyzt[idx]; }
		const type* data(unsigned int idx = 0) const { return &m_xyzt[idx]; }
		      type* xyzt(unsigned int idx = 0)       { return &m_xyzt[idx]; }
		      type* data(unsigned int idx = 0)       { return &m_xyzt[idx]; }


		type trace() const { return (m_xyzt[0] + m_xyzt[5] + m_xyzt[10] + m_xyzt[15]); }
		type det() const {
			type v = type(0);

			#define a(row, col) (*this)[(col - 1) * 4 + (row - 1)]
			v += (a(1,1) * a(2,2) * a(3,3) * a(4,4));  v += (a(1,1) * a(2,3) * a(3,4) * a(4,2));  v += (a(1,1) * a(2,4) * a(3,2) * a(4,3));
			v += (a(1,2) * a(2,1) * a(3,4) * a(4,3));  v += (a(1,2) * a(2,3) * a(3,1) * a(4,4));  v += (a(1,2) * a(2,4) * a(3,3) * a(4,1));
			v += (a(1,3) * a(2,1) * a(3,2) * a(4,4));  v += (a(1,3) * a(2,2) * a(3,4) * a(4,1));  v += (a(1,3) * a(2,4) * a(3,1) * a(4,2));
			v += (a(1,4) * a(2,1) * a(3,3) * a(4,2));  v += (a(1,4) * a(2,2) * a(3,1) * a(4,3));  v += (a(1,4) * a(2,3) * a(3,2) * a(4,1));
			v -= (a(1,1) * a(2,2) * a(3,4) * a(4,3));  v -= (a(1,1) * a(2,3) * a(3,2) * a(4,4));  v -= (a(1,1) * a(2,4) * a(3,3) * a(4,2));
			v -= (a(1,2) * a(2,1) * a(3,3) * a(4,4));  v -= (a(1,2) * a(2,3) * a(3,4) * a(4,1));  v -= (a(1,2) * a(2,4) * a(3,1) * a(4,3));
			v -= (a(1,3) * a(2,1) * a(3,4) * a(4,2));  v -= (a(1,3) * a(2,2) * a(3,1) * a(4,4));  v -= (a(1,3) * a(2,4) * a(3,2) * a(4,1));
			v -= (a(1,4) * a(2,1) * a(3,2) * a(4,3));  v -= (a(1,4) * a(2,2) * a(3,3) * a(4,1));  v -= (a(1,4) * a(2,3) * a(3,1) * a(4,2));
			#undef a

			return v;
		}


		// column accessors (NOTE: the t-vector is really a point!)
		m_vec_type get_x_vec() const { m_vec_type v; v.x() = m_xyzt[ 0]; v.y() = m_xyzt[ 1]; v.z() = m_xyzt[ 2]; v.w() = m_xyzt[ 3]; return v; }
		m_vec_type get_y_vec() const { m_vec_type v; v.x() = m_xyzt[ 4]; v.y() = m_xyzt[ 5]; v.z() = m_xyzt[ 6]; v.w() = m_xyzt[ 7]; return v; }
		m_vec_type get_z_vec() const { m_vec_type v; v.x() = m_xyzt[ 8]; v.y() = m_xyzt[ 9]; v.z() = m_xyzt[10]; v.w() = m_xyzt[11]; return v; }
		m_vec_type get_t_vec() const { m_vec_type v; v.x() = m_xyzt[12]; v.y() = m_xyzt[13]; v.z() = m_xyzt[14]; v.w() = m_xyzt[15]; return v; }
		m_pos_type get_t_pos() const { m_pos_type p; p.x() = m_xyzt[12]; p.y() = m_xyzt[13]; p.z() = m_xyzt[14]; p.w() = m_xyzt[15]; return p; }

		// random accessors (<idx> should be one of {0,1,2,3})
		m_vec_type get_row_vec(unsigned int idx) const { return (m_vec_type(m_xyzt[idx        ], m_xyzt[idx     + 4], m_xyzt[idx     + 8], m_xyzt[idx     + 12])); }
		m_vec_type get_col_vec(unsigned int idx) const { return (m_vec_type(m_xyzt[idx * 4 + 0], m_xyzt[idx * 4 + 1], m_xyzt[idx * 4 + 2], m_xyzt[idx * 4 +  3])); }


		m_mat_type orthonormalize() const {
			m_mat_type r = *this;
			m_vec_type xv = r.get_x_vec();
			m_vec_type yv = r.get_y_vec();
			m_vec_type zv;

			xv = xv.normalize_ref(); zv = xv.outer(yv);
			zv = zv.normalize_ref(); yv = zv.outer(xv);
			yv = yv.normalize_ref();

			r.set_x_vec(xv);
			r.set_y_vec(yv);
			r.set_z_vec(zv);
			return r;
		}

		m_mat_type translate(const m_vec_type& v) const {
			m_mat_type r = *this;
			r[12] += ((v.x() * r[0]) + (v.y() * r[4]) + (v.z() * r[ 8])); // same as tv.inner(rows[0])
			r[13] += ((v.x() * r[1]) + (v.y() * r[5]) + (v.z() * r[ 9])); // same as tv.inner(rows[1])
			r[14] += ((v.x() * r[2]) + (v.y() * r[6]) + (v.z() * r[10])); // same as tv.inner(rows[2])
			r[15] += ((v.x() * r[3]) + (v.y() * r[7]) + (v.z() * r[11])); // same as tv.inner(rows[3])
			return r;
		}
		m_mat_type scale(const m_vec_type& sv) const {
			m_mat_type r = *this;

			r[ 0] *= sv.x(); r[ 4] *= sv.y();
			r[ 1] *= sv.x(); r[ 5] *= sv.y();
			r[ 2] *= sv.x(); r[ 6] *= sv.y();

			r[ 8] *= sv.z(); // r[12] *= sv.w();
			r[ 9] *= sv.z(); // r[13] *= sv.w();
			r[10] *= sv.z(); // r[14] *= sv.w();

			return r;
		}

		m_mat_type transpose_rotation() const {
			m_mat_type r = *this;
			std::swap(r[1], r[4]);
			std::swap(r[2], r[8]);
			std::swap(r[6], r[9]);
			return r;
		}
		m_mat_type transpose_translation() const {
			m_mat_type r = *this;
			std::swap(r[ 3], r[12]);
			std::swap(r[ 7], r[13]);
			std::swap(r[11], r[14]);
			return r;
		}
		m_mat_type transpose() const {
			m_mat_type r = *this;
			r.transpose_rotation_ref();
			r.transpose_translation_ref();
			return r;
		}

		// assumes this matrix only performs translation and rotation
		m_mat_type invert_affine() const {
			m_mat_type r = *this;

			// transpose the rotation
			r.transpose_rotation_ref();
			// set the inverse translation; R*v ignores T
			r.set_t_vec(r * -r.get_t_vec());
			return r;
		}
		// generalized inverse for non-orthonormal 4x4 matrices
		// A^-1 = (1 / det(A)) (C^T)_{ij} = (1 / det(A)) C_{ji}
		// where C is the matrix of cofactors
		m_mat_type invert_projective(const type eps = M_FEPS) const;

		m_mat_type& orthonormalize_ref() { return ((*this) = orthonormalize()); }
		m_mat_type& translate_ref(const m_vec_type& tv) { return ((*this) = translate(tv)); }
		m_mat_type& scale_ref(const m_vec_type& sv) { return ((*this) = scale(sv)); }

		m_mat_type& transpose_ref() { return ((*this) = transpose()); }
		m_mat_type& transpose_rotation_ref() { return ((*this) = transpose_rotation()); }
		m_mat_type& transpose_translation_ref() { return ((*this) = transpose_translation()); }
		m_mat_type& invert_affine_ref() { return ((*this) = invert_affine()); }



		m_mat_type& rotate_x_ref(type angle) { return ((*this) = rotate_x(angle)); }
		m_mat_type& rotate_y_ref(type angle) { return ((*this) = rotate_y(angle)); }
		m_mat_type& rotate_z_ref(type angle) { return ((*this) = rotate_z(angle)); }

		// these perform rotations around the (fixed aka external aka extrinsic) global coordinate axes
		// for any external rotation, multiply by the matrix R(C)*R(B)*R(A) if rotation order is A,B,C
		m_mat_type& rotate_xyz_ext_ref(const m_vec_type& angles) { return ((*this) = rotate_xyz_ext(angles)); }
		m_mat_type& rotate_yxz_ext_ref(const m_vec_type& angles) { return ((*this) = rotate_yxz_ext(angles)); }
		m_mat_type& rotate_zxy_ext_ref(const m_vec_type& angles) { return ((*this) = rotate_zxy_ext(angles)); }
		m_mat_type& rotate_zyx_ext_ref(const m_vec_type& angles) { return ((*this) = rotate_zyx_ext(angles)); }
		// these perform rotations around the (moving aka internal aka intrinsic) local coordinate axes
		// for any internal rotation, multiply by the matrix R(A)*R(B)*R(C) if rotation order is A,B,C
		// (or equivalently by the external matrix corresponding to the rotation order C,B,A)
		m_mat_type& rotate_xyz_int_ref(const m_vec_type& angles) { return ((*this) = rotate_xyz_int(angles)); }
		m_mat_type& rotate_yxz_int_ref(const m_vec_type& angles) { return ((*this) = rotate_yxz_int(angles)); }
		m_mat_type& rotate_zxy_int_ref(const m_vec_type& angles) { return ((*this) = rotate_zxy_int(angles)); }
		m_mat_type& rotate_zyx_int_ref(const m_vec_type& angles) { return ((*this) = rotate_zyx_int(angles)); }


		// rotate in X,Y,Z order [Rext = R(Z)*R(Y)*R(X), Rint = R(X)*R(Y)*R(Z)]
		m_mat_type rotate_xyz_ext(const m_vec_type angles) const {
			const m_mat_type m  = *this;
			const m_mat_type rx = std::move(compose_x_rotation(angles[P_AXIS_IDX]));
			const m_mat_type ry = std::move(compose_y_rotation(angles[Y_AXIS_IDX]));
			const m_mat_type rz = std::move(compose_z_rotation(angles[R_AXIS_IDX]));
			return (m * (rz * ry * rx));
		}
		m_mat_type rotate_xyz_int(const m_vec_type angles) const {
			return (rotate_zyx_ext(angles));
		}


		// rotate in Y,X,Z order [Rext = R(Z)*R(X)*R(Y), Rint = R(Y)*R(X)*R(Z)]
		m_mat_type rotate_yxz_ext(const m_vec_type angles) const {
			const m_mat_type m  = *this;
			const m_mat_type rx = std::move(compose_x_rotation(angles[P_AXIS_IDX]));
			const m_mat_type ry = std::move(compose_y_rotation(angles[Y_AXIS_IDX]));
			const m_mat_type rz = std::move(compose_z_rotation(angles[R_AXIS_IDX]));
			return (m * (rz * rx * ry));
		}
		m_mat_type rotate_yxz_int(const m_vec_type angles) const {
			return (rotate_zxy_ext(angles));
		}


		// rotate in Z,X,Y order [Rext = R(Y)*R(X)*R(Z), Rint = R(Z)*R(X)*R(Y)]
		m_mat_type rotate_zxy_ext(const m_vec_type angles) const {
			const m_mat_type m  = *this;
			const m_mat_type rx = std::move(compose_x_rotation(angles[P_AXIS_IDX]));
			const m_mat_type ry = std::move(compose_y_rotation(angles[Y_AXIS_IDX]));
			const m_mat_type rz = std::move(compose_z_rotation(angles[R_AXIS_IDX]));
			return (m * (ry * rx * rz));
		}
		m_mat_type rotate_zxy_int(const m_vec_type angles) const {
			return (rotate_yxz_ext(angles));
		}


		// rotate in Z,Y,X order [Rext = R(X)*R(Y)*R(Z), Rint = R(Z)*R(Y)*R(X)]
		m_mat_type rotate_zyx_ext(const m_vec_type angles) const {
			const m_mat_type m  = *this;
			const m_mat_type rx = std::move(compose_x_rotation(angles[P_AXIS_IDX]));
			const m_mat_type ry = std::move(compose_y_rotation(angles[Y_AXIS_IDX]));
			const m_mat_type rz = std::move(compose_z_rotation(angles[R_AXIS_IDX]));
			return (m * (rx * ry * rz));
		}
		m_mat_type rotate_zyx_int(const m_vec_type angles) const {
			return (rotate_xyz_ext(angles));
		}


		m_mat_type rotate_x(type angle) const { return ((*this) * compose_x_rotation(angle)); }
		m_mat_type rotate_y(type angle) const { return ((*this) * compose_y_rotation(angle)); }
		m_mat_type rotate_z(type angle) const { return ((*this) * compose_z_rotation(angle)); }


		// rotate in local or global YZ-plane; phi or alpha radians
		// constructs but does *not* apply Rx, so it can be chained
		static m_mat_type compose_x_rotation(type angle) {
			angle = clamp_angle_rad(angle);

			// Rx = ([X = [1, 0, 0], Y = [0, ca, -sa], Z = [0, sa, ca]])
			const type ca = std::cos(angle);
			const type sa = std::sin(angle);

			m_mat_type rm;
			rm[ 5] = +ca; rm[ 9] = +sa;
			rm[ 6] = -sa; rm[10] = +ca;
			return rm;
		}

		// rotate in local or global XZ-plane; theta or beta radians
		// constructs but does *not* apply Ry, so it can be chained
		static m_mat_type compose_y_rotation(type angle) {
			angle = clamp_angle_rad(angle);

			// Ry = ([X = [ca, 0, sa], Y = [0, 1, 0], Z = [-sa, 0, ca]])
			const type ca = std::cos(angle);
			const type sa = std::sin(angle);

			m_mat_type rm;
			rm[ 0] = +ca; rm[ 8] = -sa;
			rm[ 2] = +sa; rm[10] = +ca;
			return rm;
		}

		// rotate in local or global XY-plane; psi or gamma radians
		// constructs but does *not* apply Rz, so it can be chained
		static m_mat_type compose_z_rotation(type angle) {
			angle = clamp_angle_rad(angle);

			// Rz = ([X = [ca, -sa, 0], Y = [sa, ca, 0], Z = [0, 0, 1]])
			const type ca = std::cos(angle);
			const type sa = std::sin(angle);

			m_mat_type rm;
			rm[0] = +ca; rm[4] = +sa;
			rm[1] = -sa; rm[5] = +ca;
			return rm;
		}



		// rotation about any arbitrary axis, ie. in the plane crossing
		// the origin whose normal is given by <axis.x, axis.y, axis.z>
		// any such rotation can be decomposed into rotations about the
		// three component axes (by aligning the arbitrary axis with one
		// of the coordinate axes [here z], then rotating, then undoing
		// step 1)
		//
		// NOTE: rot_axis equal to Y is not supported, use rotate_y_ref
		m_mat_type& rotate_axis_ref(const m_vec_type& axis_angle) {
			return ((*this) = rotate_axis(axis_angle));
		}

		m_mat_type rotate_axis(const m_vec_type& axis_angle) const {
			const m_vec_type& ay = m_vec_type::y_vector();
			const m_vec_type  vx = (axis_angle.outer(ay)).normalize(); // ignores .w
			const m_vec_type  vy = (axis_angle.outer(vx)).normalize();

			const m_mat_type m_fwd = std::move(m_mat_type(vx, vy, axis_angle)); // ignores z.w
			const m_mat_type m_inv = std::move(m_fwd.invert_affine());
			const m_mat_type m_rot = std::move(m_mat_type::compose_z_rotation(axis_angle.w()));

			return ((*this) * (m_fwd * m_rot * m_inv));
		}

		#if 0
		m_mat_type& rotate_axis_ref_tmp(const m_vec_type& axis_angle) {
			const type ca = std::cos(axis_angle.w());
			const type sa = std::sin(axis_angle.w());

			const m_vec_type& xyz_axis = m_vec_type::xyz_vector();
			const m_vec_type    x_axis = get_x_vec();
			const m_vec_type    y_axis = get_y_vec();
			const m_vec_type    z_axis = get_z_vec();

			// project rotation axis onto each principal axis; dot by default ignores .w (mul does not)
			const m_vec_type fwd_axis_x = (axis_angle * xyz_axis) * x_axis.inner(axis_angle);
			const m_vec_type fwd_axis_y = (axis_angle * xyz_axis) * y_axis.inner(axis_angle);
			const m_vec_type fwd_axis_z = (axis_angle * xyz_axis) * z_axis.inner(axis_angle);

			// NOTE: rgt and up define the rotational plane
			const m_vec_type rgt_axis_x = x_axis - fwd_axis_x;
			const m_vec_type rgt_axis_y = y_axis - fwd_axis_y;
			const m_vec_type rgt_axis_z = z_axis - fwd_axis_z;

			// does not preserve orthonormality for non-principal rotation axes
			const m_vec_type up_axis_x = (axis_angle.outer(rgt_axis_x)).normalize();
			const m_vec_type up_axis_y = (axis_angle.outer(rgt_axis_y)).normalize();
			const m_vec_type up_axis_z = (axis_angle.outer(rgt_axis_z)).normalize();

			// Rodriguez's (?) equation
			// r=point (to rotate), n=axis, n^t*r is dot-product in matrix-notation
			// then r' = n^t*r*n + cos(angle)*(r - n^t*r*n) + sin(angle)*(n cross r)
			//
			set_x_vec((fwd_axis_x + (rgt_axis_x * ca + up_axis_x * sa)));
			set_y_vec((fwd_axis_y + (rgt_axis_y * ca + up_axis_y * sa)));
			set_z_vec((fwd_axis_z + (rgt_axis_z * ca + up_axis_z * sa)));

			return *this;
		}
		#endif



		// note: v.w() is ignored
		m_mat_type& set_x_vec(const m_vec_type& v) { m_xyzt[ 0] = v.x(); m_xyzt[ 1] = v.y(); m_xyzt[ 2] = v.z(); m_xyzt[ 3] = type(0); return *this; }
		m_mat_type& set_y_vec(const m_vec_type& v) { m_xyzt[ 4] = v.x(); m_xyzt[ 5] = v.y(); m_xyzt[ 6] = v.z(); m_xyzt[ 7] = type(0); return *this; }
		m_mat_type& set_z_vec(const m_vec_type& v) { m_xyzt[ 8] = v.x(); m_xyzt[ 9] = v.y(); m_xyzt[10] = v.z(); m_xyzt[11] = type(0); return *this; }
		m_mat_type& set_t_vec(const m_vec_type& v) { m_xyzt[12] = v.x(); m_xyzt[13] = v.y(); m_xyzt[14] = v.z(); m_xyzt[15] = type(1); return *this; }

		m_mat_type& set_row_vec(unsigned int idx, const m_vec_type& v) { m_xyzt[idx     + 0] = v.x(); m_xyzt[idx     + 4] = v.y(); m_xyzt[idx     + 8] = v.z(); m_xyzt[idx     + 12] = v.w(); return *this; }
		m_mat_type& set_col_vec(unsigned int idx, const m_vec_type& v) { m_xyzt[idx * 4 + 0] = v.x(); m_xyzt[idx * 4 + 1] = v.y(); m_xyzt[idx * 4 + 2] = v.z(); m_xyzt[idx * 4 +  3] = v.w(); return *this; }

		m_mat_type& set_vectors(const m_vec_type& vp, unsigned int idx = 2) {
			m_vec_type va;

			// given a primary axis <vp> (interpreted as a rotated
			// identity-matrix x-, y-, or z-vector), construct the
			// orthogonal secondary and tertiary vectors through an
			// auxiliary
			//
			// if <vp> is +x, secondary is +y and tertiary is +z
			// if <vp> is +y, secondary is +z and tertiary is +x
			// if <vp> is +z, secondary is +x and tertiary is +y
			// (note the cyclicity)
			switch (idx) {
				case 0: { va = m_vec_type::z_vector(); } break;
				case 1: { va = m_vec_type::x_vector(); } break;
				case 2: { va = m_vec_type::y_vector(); } break;
				default: { assert(false); } break;
			}

			// auxiliary should not be colinear with primary
			assert(std::fabs(va.inner(vp)) < type(1));

			// v.outer(w) always points to the local "right"
			const m_vec_type vs = (vp.outer(va)).normalize();
			const m_vec_type vt = (vs.outer(vp)).normalize();

			switch (idx) {
				case 0: {
					set_x_vec( vp);
					set_y_vec(-vs);
					set_z_vec( vt);
				} break;
				case 1: {
					set_x_vec( vt);
					set_y_vec( vp);
					set_z_vec(-vs);
				} break;
				case 2: {
					set_x_vec(-vs);
					set_y_vec( vt);
					set_z_vec( vp);
				} break;
				default: {
					assert(false);
				} break;
			}

			return *this;
		}


		m_mat_type& add_raw_matrix(const type* m) {
			for (unsigned int i = 0; i < MATRIX_SIZE; i += 4) {
				m_xyzt[i + 0] += m[i + 0];
				m_xyzt[i + 1] += m[i + 1];
				m_xyzt[i + 2] += m[i + 2];
				m_xyzt[i + 3] += m[i + 3];
			}
			return *this;
		}

		m_mat_type& set_raw_matrix(const type* m) {
			for (unsigned int i = 0; i < MATRIX_SIZE; i += 4) {
				m_xyzt[i + 0] = m[i + 0];
				m_xyzt[i + 1] = m[i + 1];
				m_xyzt[i + 2] = m[i + 2];
				m_xyzt[i + 3] = m[i + 3];
			}
			return *this;
		}

		m_mat_type& set_diag_matrix(type diag_elem) {
			for (unsigned int n = 0; n < MATRIX_SIZE; n += 4) {
				m_xyzt[n + 0] = diag_elem * type(n ==  0);
				m_xyzt[n + 1] = diag_elem * type(n ==  4);
				m_xyzt[n + 2] = diag_elem * type(n ==  8);
				m_xyzt[n + 3] = diag_elem * type(n == 12);
			}

			return *this;
		}
		m_mat_type& set_null_matrix() { return (set_diag_matrix(0)); }
		m_mat_type& set_unit_matrix() { return (set_diag_matrix(1)); }

		void print(const char* tabs = "") const;


		unsigned int is_identity(const type eps = M_FEPS) const {
			const m_vec_type& xv = get_x_vec();
			const m_vec_type& yv = get_y_vec();
			const m_vec_type& zv = get_z_vec();
			const m_vec_type& tv = get_t_vec();

			unsigned int n = 0;

			n += ((!xv.equals(m_vec_type::x_vector(), m_vec_type::x_vector() * eps)) * (1 << 0));
			n += ((!yv.equals(m_vec_type::y_vector(), m_vec_type::y_vector() * eps)) * (1 << 1));
			n += ((!zv.equals(m_vec_type::z_vector(), m_vec_type::z_vector() * eps)) * (1 << 2));
			n += ((!tv.equals(m_vec_type::w_vector(), m_vec_type::w_vector() * eps)) * (1 << 3));

			return n;
		}

		unsigned int is_orthonormal(const type eps = M_FEPS) const {
			const m_vec_type& xv = get_x_vec();
			const m_vec_type& yv = get_y_vec();
			const m_vec_type& zv = get_z_vec();

			unsigned int n = 0;

			// test angles
			n += ((type(1) - fp_eq(xv.inner(yv), type(0), eps)) * (1 << 0));
			n += ((type(1) - fp_eq(yv.inner(zv), type(0), eps)) * (1 << 1));
			n += ((type(1) - fp_eq(xv.inner(zv), type(0), eps)) * (1 << 2));
			// test lengths
			n += ((type(1) - fp_eq(xv.sq_len(), type(1), eps)) * (1 << 3));
			n += ((type(1) - fp_eq(yv.sq_len(), type(1), eps)) * (1 << 4));
			n += ((type(1) - fp_eq(zv.sq_len(), type(1), eps)) * (1 << 5));

			return n;
		}


		m_vec_type get_angles_rh(const type eps = M_FEPS) const;
		m_vec_type get_angles_lh(const type eps = M_FEPS) const;


		static m_mat_type lerp(const m_mat_type& src_mat, const m_mat_type& dst_mat, const m_vec_type& alpha) {
			m_mat_type r;

			#if 0
			// interpolate translation independently
			const m_vec_type& int_tvec = lib_math::lerp(src_mat.get_t_vec(), dst_mat.get_t_vec(), alpha.y());
			const m_vec_type& int_zvec = lib_math::slerp(src_mat.get_z_vec(), dst_mat.get_z_vec(), alpha.x(), 0.00001f);
			const m_vec_type& int_yvec = lib_math::slerp(src_mat.get_y_vec(), dst_mat.get_y_vec(), alpha.x(), 0.00001f);
			const m_vec_type  int_xvec = int_zvec.outer(int_yvec);

			// construct lerp'ed matrix
			r.set_x_vec( int_xvec.normalize());
			r.set_y_vec((int_zvec.outer(int_xvec)).normalize());
			r.set_z_vec( int_zvec.normalize());
			r.set_t_vec( int_tvec);
			#endif

			#if 1
			const t_quaternion<type>& src_quat = t_quaternion<type>::rot_matrix_quat_cm(src_mat);
			const t_quaternion<type>& dst_quat = t_quaternion<type>::rot_matrix_quat_cm(dst_mat);
			const t_quaternion<type>& lrp_quat = t_quaternion<type>::slerp(src_quat, dst_quat, alpha.x(), 0.00005f);
			const         m_mat_type& lrp_matr = lrp_quat.to_rotation_matrix();

			r.set_x_vec(lrp_matr.get_x_vec());
			r.set_y_vec(lrp_matr.get_y_vec());
			r.set_z_vec(lrp_matr.get_z_vec());
			r.set_t_vec(lib_math::lerp(src_mat.get_t_vec(), dst_mat.get_t_vec(), alpha.y()));
			#endif

			return r;
		}

		static m_mat_type skew_sym(const m_vec_type& v) {
			m_mat_type m;
			m.set_x_vec(m_vec_type(type(0),   v.z(),  -v.y()));
			m.set_y_vec(m_vec_type( -v.z(), type(0),   v.x()));
			m.set_z_vec(m_vec_type(  v.y(),  -v.x(), type(0)));
			return m;
		}

		// same as vector::outer(v, w); useful for generalizing to higher dimensions (m is skew-symmetric)
		//   [   0 -v.z  v.y]   [w.x]   [   0 * w.x  -  v.z * w.y  +  v.y * w.z]   [v.y * w.z - v.z * w.y]
		//   [ v.z    0 -v.x] * [w.y] = [ v.z * w.x  +    0 * w.y  -  v.x * w.z] = [v.z * w.x - v.x * w.z]
		//   [-v.y  v.x    0]   [w.z]   [-v.y * w.x  +  v.x * w.y  +    0 * w.z]   [v.x * w.y - v.y * w.x]
		static m_vec_type outer_product(const m_vec_type& v, const m_vec_type& w) {
			return (skew_sym(v) * w);
		}

	private:
		// stores the elements in column-major order
		// m_xyzt[0 .. 3] represents the x-axis, etc
		type m_xyzt[MATRIX_SIZE];
	};


	typedef t_matrix44t< float> t_matrix44f;
	typedef t_matrix44t<double> t_matrix44d;

	// aliases for less typing
	typedef t_matrix44f t_mat44f;
	typedef t_matrix44d t_mat44d;
};

#endif

