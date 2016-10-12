#ifndef LIBGEOM_QUATERNION_HDR
#define LIBGEOM_QUATERNION_HDR

#include <cmath>
#include <cassert>
#include <cstddef>

#include "./math_defs.hpp"
#include "./vector.hpp"

namespace lib_math {
	template<typename type> struct t_quaternion {
	public:
		t_quaternion(type _a = type(0), type _b = type(0), type _c = type(0), type _d = type(0)) {
			a() = _a; b() = _b;
			c() = _c; d() = _d;
		}
		t_quaternion(type a, const t_vector4t<type>&  bcd) { *this = std::move(t_quaternion<type>(a       ,  bcd.x(),  bcd.y(),  bcd.z())); }
		t_quaternion(        const t_vector4t<type>& abcd) { *this = std::move(t_quaternion<type>(abcd.x(), abcd.y(), abcd.z(), abcd.w())); }


		#if 0
		t_quaternion<type> operator + (const t_quaternion& q) const { return (t_quaternion(abcd_vector() + q.abcd_vector())); }
		t_quaternion<type> operator - (const t_quaternion& q) const { return (t_quaternion(abcd_vector() - q.abcd_vector())); }
		#else
		t_quaternion<type> operator + (const t_quaternion& q) const { return (t_quaternion(a() + q.a(), b() + q.b(), c() + q.c(), d() + q.d())); }
		t_quaternion<type> operator - (const t_quaternion& q) const { return (t_quaternion(a() - q.a(), b() - q.b(), c() - q.c(), d() - q.d())); }
		#endif

		// Hamilton product = (a1 + b1*i + c1*j + d1*k) x (a2 + b2*i + c2*j + d2*k)
		#if 0
		t_quaternion<type> operator * (const t_quaternion& q) const {
			const t_vector4t<type>&  bcd =   bcd_vector();
			const t_vector4t<type>& qbcd = q.bcd_vector();
			const t_vector4t<type>  rbcd = (qbcd * a()) + (bcd * q.a()) + bcd.outer(qbcd);
			return (t_quaternion(a() * q.a() - bcd.inner(qbcd), rbcd));
		}
		#else
		t_quaternion<type> operator * (const t_quaternion& q) const {
			t_quaternion<type> r;
			r.a() = (a() * q.a())  -  (b() * q.b())  -  (c() * q.c())  -  (d() * q.d());
			r.b() = (a() * q.b())  +  (b() * q.a())  +  (c() * q.d())  -  (d() * q.c());
			r.c() = (a() * q.c())  -  (b() * q.d())  +  (c() * q.a())  +  (d() * q.b());
			r.d() = (a() * q.d())  +  (b() * q.c())  -  (c() * q.b())  +  (d() * q.a());
			return r;
		}
		#endif

		t_quaternion<type> operator * (const type s) const { return (t_quaternion(a() * s, b() * s, c() * s, d() * s)); }
		t_quaternion<type> operator / (const type s) const { return ((*this) * (type(1) / s)); }

		// conjugate; -q = q* = -0.5 * (q + i*q*i + j*q*j + k*q*k)
		// scalar part in terms of the conjugate is (q + (-q)) * 0.5
		// vector part in terms of the conjugate is (q - (-q)) * 0.5
		// additive inverse of a vector bcd in R3 is the same as its conjugate in
		// quaternion form; geometry of R3 is reflected in the algebraic structure
		// of quaternions
		t_quaternion<type> operator - () const {
			return (t_quaternion(a(), -b(), -c(), -d()));
		}


		// rotate an R3 vector <v> by interpreting it as a pure quaternion
		// any unit-length quaternion q = (a = w, r = {b, c, d}) represents a rotation;
		// the product q*v*-q can also be written as v + cross(2*r, (cross(r,v) + v*w))
		// other methods are "matrix * vector" or "(axis-angle to matrix) * vector"
		#if 0
		t_quaternion<type> rotate_vector(const t_vector4t<type>& v) const {
			const t_vector4t<type>& bcd = bcd_vector();
			const t_vector4t<type>  lhs = bcd * type(2);
			const t_vector4t<type>  rhs = bcd.outer(v) + v * a();
			return (v + lhs.outer(rhs));
		}
		#else
		t_quaternion<type> rotate_vector(const t_vector4t<type>& v) const {
			return (((*this) * t_quaternion(type(0), v)) * -(*this));
		}
		#endif


		// only valid for a quaternion whose norm is close to 1
		// assume  (1 + eps) = |q|, then eps = |q|-1 and we get
		//     sqrt(1 + eps) ~=   (1 + (  eps)*0.5)
		//   1/sqrt(1 + eps) ~=   (1 - (  eps)*0.5) = 1 - (|q|-1)*0.5
		//   1/sqrt(1 + eps) ~= 1/(1 - (|q|-1)*0.5) = 1.5 - 0.5 * |q|*|q|
		t_quaternion<type>  anormalize() const { return ( (*this) * (type(1.5) - type(0.5) * sq_norm())); }
		t_quaternion<type>   normalize() const { return ( (*this) * (type(1) / norm())); }
		t_quaternion<type>  reciprocal() const { return (-(*this) * (type(1) / sq_norm())); }
		t_quaternion<type>& anormalize_ref() { *this = anormalize(); }
		t_quaternion<type>&  normalize_ref() { *this =  normalize(); }
		t_quaternion<type>& reciprocal_ref() { *this = reciprocal(); }
		t_quaternion<type>& axis_angle_ref(const t_vector4t<type>& axis_angle) { *this = axis_angle_quat(axis_angle); }

		t_quaternion<type> bcd_multiply(const t_quaternion<type>& q) const {
			const t_vector4t<type>&  bcd =   bcd_vector();
			const t_vector4t<type>& qbcd = q.bcd_vector();
			return (t_quaternion(type(0), bcd.outer(qbcd) - bcd.inner(qbcd)));
		}


		static t_quaternion<type> rot_matrix_quat_cm(const t_matrix44t<type>& rot_matrix) {
			return (rot_matrix_quat_rm(rot_matrix.transpose()));
		}
		// convert a (assumed row-major) rotation-matrix to a unit-quaternion
		static t_quaternion<type> rot_matrix_quat_rm(const t_matrix44t<type>& rot_matrix) {
			t_quaternion<type> q;

			const t_vector4t<type>& xv = rot_matrix.get_x_vec();
			const t_vector4t<type>& yv = rot_matrix.get_y_vec();
			const t_vector4t<type>& zv = rot_matrix.get_z_vec();

			// we can write the conversion as
			//   ca = sqrt(tr)
			//   sa = 0.5/ca
			//    w = ca*0.5
			// or as
			//   ca' = 1.0/sqrt(tr)
			//   sa' = ca'*0.5
			//    w' = sa' * tr
			// since
			//    0.5/ca = 1.0/(2.0*ca) = (1.0/ca) * 0.5 = (1.0/sqrt(tr)) * 0.5 = ca' * 0.5 = sa'
			//    ca*0.5 = sqrt(tr)*0.5 = (tr/sqrt(tr))*0.5 = ((1.0/sqrt(tr))*0.5)*tr = (ca'*0.5)*tr = sa'*tr
			//
			// const type ca = std::sqrt(type(1) + (xv.x() + yv.y() + zv.z()));
			// const type sa = type(0.5) / ca;
			//
			// note: a 180-degree rotation around (e.g.) the y-axis can cause tr=0
			// if this happens the quaternion also has to be explicitly normalized
			const type tr = type(1) + (xv.x() + yv.y() + zv.z());
			const type ca = lib_math::inv_sqrt_sse(std::max(type(1e-4), tr));
			const type sa = ca * type(0.5);

			q.a() = sa * tr;
			q.b() = sa * (zv.y() - yv.z());
			q.c() = sa * (xv.z() - zv.x());
			q.d() = sa * (yv.x() - xv.y());
			return (q.normalize());
		}

		// encode a rotation of <angle> radians around the unit-vector <axis>
		// when rotating a vector we multiply it both by this quaternion and
		// its conjugate, rotating into and then out of the fourth dimension
		// (with each multiplication covering half the total angle)
		static t_quaternion<type> axis_angle_quat(const t_vector4t<type>& axis_angle) {
			t_quaternion<type> q;

			// exp(theta*0.5 * (dot(axis, {i,j,k}))) = cos(theta*0.5) + dot(axis, {i,j,k}) * sin(theta*0.5)
			const type ca = std::cos(axis_angle.w() * type(0.5));
			const type sa = std::sin(axis_angle.w() * type(0.5));

			q.a() = ca;
			q.b() = sa * axis_angle.x();
			q.c() = sa * axis_angle.y();
			q.d() = sa * axis_angle.z();
			return q;
		}


		static t_quaternion<type> nlerp(const t_quaternion& q0, const t_quaternion& q1, type wgt, type /*eps*/) {
			// quaternions q=(a,v) and -q=(-a,-v) map to the same rotation; pick the short path
			const type w0 = (type(1) - wgt);
			const type w1 = (type(0) + wgt) * lib_math::sign(q0.calc_dot(q1));

			return ((q0 * w0 + q1 * w1).normalize());
		}

		#if 0
		static t_quaternion<type> slerp(const t_quaternion& q0, const t_quaternion& q1, type wgt, type eps) {
			const type  cosa = lib_math::clamp(q0.calc_dot(q1), type(-1), type(1));
			const type  sign = lib_math::sign(cosa);
			const type acosa = std::acos(cosa * sign);

			// linearly interpolate (with normalization) if q0 and q1 are ~colinear
			// otherwise spherically interpolate; for any unit quaternion in versor
			// form q = cos(a) + v * sin(a) we get q^t = cos(a * t) + v * sin(a * t)
			// where a * t = angle * wgt
			//
			// [ method | commutative | const-vel | min-torque]
			// [ qslerp |          no |       yes |        yes]
			// [ qnlerp |         yes |        no |        yes]
			// [ qllerp |         yes |       yes |         no]
			//
			if ((cosa * sign) > (type(1) - eps))
				return (nlerp(q0, q1, wgt, eps));

			// note: can also express slerp(q0, q1, t) as (q1 * inverse(q0))^t * q0
			const t_quaternion<type> q2 = (q1 - (q0 * cosa)) / std::sin(acosa); // v1 - v0*dot(v0,v1)
			const t_quaternion<type> qi = q0 * std::cos(acosa * wgt) + q2 * std::sin(acosa * wgt * sign);
			return (qi.normalize());
		}
		#else
		static t_quaternion<type> slerp(const t_quaternion& q0, const t_quaternion& q1, type wgt, type eps) {
			const type cosa = lib_math::clamp(q0.calc_dot(q1), type(-1), type(1));
			const type sign = lib_math::sign(cosa);

			if ((cosa * sign) < (type(1) - eps)) {
				const type acosa = std::acos(cosa * sign);
				const type isina = type(1) / std::sin(acosa);

				const type w0 = std::sin((type(1) - wgt) * acosa) * isina;
				const type w1 = std::sin((type(0) + wgt) * acosa) * isina;

				return ((q0 * w0) + (q1 * w1 * sign));
			}

			return (nlerp(q0, q1, wgt, eps));
		}
		#endif


		// convert a unit-length quaternion to an (xyz=axis, w=angle) vector
		t_vector4t<type> to_axis_angle_vector() const {
			t_vector4t<type> v;

			const type angle = calc_angle();
			const type scale = type(1) / std::sin(angle * type(0.5));

			v.x() = scale * b();
			v.y() = scale * c();
			v.z() = scale * d();
			v.w() = angle;
			return v;
		}

		// convert a unit-quaternion to a (column-major) rotation-matrix
		t_matrix44t<type> to_rotation_matrix() const {
			t_matrix44t<type> r;

			constexpr type t1 = 1;
			constexpr type t2 = 2;

			const type bb = b() * b();
			const type cc = c() * c();
			const type dd = d() * d();

			#if 1
			const type aa = a() * a() * t1; // use t1 to kill warning
			const type ab = a() * b() * t2;
			const type ac = a() * c() * t2;
			const type ad = a() * d() * t2;
			const type bc = b() * c() * t2;
			const type bd = b() * d() * t2;
			const type cd = c() * d() * t2;
			#else
			const type ab = a() * b() * t1;
			const type ac = a() * c() * t1;
			const type ad = a() * d() * t1;
			const type bc = b() * c() * t1;
			const type bd = b() * d() * t1;
			const type cd = c() * d() * t1;
			#endif

			#if 1
			r.set_x_vec({aa + bb - cc - dd, bc + ad, bd - ac});
			r.set_y_vec({bc - ad, aa - bb + cc - dd, cd + ab});
			r.set_z_vec({bd + ac, cd - ab, aa - bb - cc + dd});
			#else
			// alternative (equivalent) formulation; non-zero quaternions
			// behave like homogeneous coordinates for rotation matrices
			r.set_x_vec({t1 - t2 * (cc + dd),      t2 * (bc + ad),      t2 * (bd - ac)});
			r.set_y_vec({     t2 * (bc - ad), t1 - t2 * (bb + dd),      t2 * (cd + ab)});
			r.set_z_vec({     t2 * (bd + ac),      t2 * (cd - ab), t1 - t2 * (bb + cc)});
			#endif
			return r;
		}

		t_vector4t<type> abcd_vector() const { return (t_vector4t<type>(a(), b(), c(), d())); }
		t_vector4t<type>  bcd_vector() const { return (t_vector4t<type>(b(), c(), d())); }

		// same as the Euclidean norm if considering quaternions as
		// a vector space; note that norm(q * s) = norm(q) * abs(s)
		type    norm() const { return (std::sqrt(sq_norm())); }
		type sq_norm() const { return (calc_dot(*this)); }
		// q multiplied by its own conjugate is not a scalar here
		// type sq_norm() const { return ((*this) * -(*this)); }
		// type sq_norm() const { return (a() * a()) + sq_norm_bcd()); }
		type    norm_bcd() const { return (std::sqrt(sq_norm_bcd())); }
		type sq_norm_bcd() const { return (calc_dot(*this) - (a() * a())); }

		// both of these angle calculations require *this to be normalized
		// type calc_angle() const { return (type(2) * std::asin(norm_bcd())); }
		type calc_angle() const { return (type(2) * std::acos(a())); }
		type calc_dist(const t_quaternion& q) const { return (((*this) - q).norm()); }
		type calc_dot(const t_quaternion& q) const { return (a() * q.a() + b() * q.b() + c() * q.c() + d() * q.d()); }

		type  a() const { return m_abcd[0]; } // w
		type  b() const { return m_abcd[1]; } // x
		type  c() const { return m_abcd[2]; } // y
		type  d() const { return m_abcd[3]; } // z
		type& a()       { return m_abcd[0]; }
		type& b()       { return m_abcd[1]; }
		type& c()       { return m_abcd[2]; }
		type& d()       { return m_abcd[3]; }

		// pure quaternions have a==0 and at least one non-zero imaginary component
		// real quaternions have a!=0 and         only     zero imaginary components
		// bcd-part is imaginary but can be treated as an ordinary vector in 3D space
		#if 0
		bool is_pure() const { return (a() == type(0) && bcd_vector() != t_vector4t<type>::zero_vector()); }
		bool is_real() const { return (a() != type(0) && bcd_vector() == t_vector4t<type>::zero_vector()); }
		bool is_unit() const { return (a() == type(1) && bcd_vector() == t_vector4t<type>::zero_vector()); }
		#else
		bool is_pure() const { return (a() == type(0) && (b() != type(0) || c() != type(0) || c() != type(0))); }
		bool is_real() const { return (a() != type(0) && (b() == type(0) && c() == type(0) && c() == type(0))); }
		bool is_unit() const { return (a() == type(1) && (b() == type(0) && c() == type(0) && c() == type(0))); }
		#endif

	private:
		// non-versor form; basis elements {1,i,j,k} are 4-vectors
		//   [0] is the coefficient for the r-axis (1, 0, 0, 0)
		//   [1] is the coefficient for the i-axis (0, 1, 0, 0)*i; i*i = -1
		//   [2] is the coefficient for the j-axis (0, 0, 1, 0)*j; j*j = -1
		//   [3] is the coefficient for the k-axis (0, 0, 0, 1)*k; k*k = -1
		//   q = {a,b,c,d} x {1,i,j,k} = {a*1, b*i, c*j, d*k}
		// from the identities and i*j*k = -1 it also follows that
		//   i*j =  k, j*k =  i, k*i =  j
		//   j*i = -k, k*j = -i, i*k = -j
		type m_abcd[4];
	};


	typedef t_quaternion< float> t_quaternion4f;
	typedef t_quaternion<double> t_quaternion4d;

	typedef t_quaternion4f t_quat4f;
	typedef t_quaternion4d t_quat4d;
};

#endif

