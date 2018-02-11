#ifndef LIBGEOM_VECTOR_HDR
#define LIBGEOM_VECTOR_HDR

#include <cassert>
#include <cmath>
#include <cstdint>

#include <algorithm>

#include "./math_defs.hpp"
#include "./template_funcs.hpp"

namespace lib_math {
	template<typename type> struct t_tuple4t;
	template<typename type> struct t_point4t;
	template<typename type> struct t_matrix44t;

	// represents a free vector (magnitude & direction)
	template<typename type> struct t_vector4t {
	public:
		t_vector4t<type>(const type _x = type(0), const type _y = type(0), const type _z = type(0), const type _w = type(0)) {
			x() = _x; y() = _y;
			z() = _z; w() = _w;
		}
		t_vector4t<type>(const type* v) { *this = std::move(t_vector4t<type>(v[0], v[1], v[2], v[3])); }
		t_vector4t<type>(const t_tuple4t<type>& t) { *this = std::move(t_vector4t<type>(t.x(), t.y(), t.z(), t.w())); }
		t_vector4t<type>(const t_vector4t<type>& v) { *this = v; }

		bool operator == (const t_vector4t<type>& v) const { return ( equals(v)); }
		bool operator != (const t_vector4t<type>& v) const { return (!equals(v)); }

		t_vector4t<type> operator + (const t_vector4t<type>& v) const { return (t_vector4t<type>( x() + v.x(),  y() + v.y(),  z() + v.z(),  w() + v.w())); }
		t_vector4t<type> operator - (const t_vector4t<type>& v) const { return (t_vector4t<type>( x() - v.x(),  y() - v.y(),  z() - v.z(),  w() - v.w())); }
		t_vector4t<type> operator * (const t_vector4t<type>& v) const { return (t_vector4t<type>( x() * v.x(),  y() * v.y(),  z() * v.z(),  w() * v.w())); } // UNOFFICIAL
		t_vector4t<type> operator / (const t_vector4t<type>& v) const { return (t_vector4t<type>( x() / v.x(),  y() / v.y(),  z() / v.z(),  w() / v.w())); } // UNOFFICIAL
		t_vector4t<type> operator - (                         ) const { return (t_vector4t<type>(-x(),         -y(),         -z()        , -w()        )); }

		t_vector4t<type>& operator  = (const t_vector4t<type>& v) { x()  = v.x(); y()  = v.y(); z()  = v.z(); w()  = v.w(); return *this; }
		t_vector4t<type>& operator += (const t_vector4t<type>& v) { x() += v.x(); y() += v.y(); z() += v.z(); w() += v.w(); return *this; }
		t_vector4t<type>& operator -= (const t_vector4t<type>& v) { x() -= v.x(); y() -= v.y(); z() -= v.z(); w() -= v.w(); return *this; }
		t_vector4t<type>& operator *= (const t_vector4t<type>& v) { x() *= v.x(); y() *= v.y(); z() *= v.z(); w() *= v.w(); return *this; } // UNOFFICIAL
		t_vector4t<type>& operator /= (const t_vector4t<type>& v) { x() /= v.x(); y() /= v.y(); z() /= v.z(); w() /= v.w(); return *this; } // UNOFFICIAL

		// t_vector4t<type> operator + (const type s) const { return (t_vector4t<type>(x() + s, y() + s, z() + s, w() + s)); } // UNOFFICIAL
		// t_vector4t<type> operator - (const type s) const { return (t_vector4t<type>(x() - s, y() - s, z() - s, w() - s)); } // UNOFFICIAL
		t_vector4t<type> operator * (const type s) const { return (t_vector4t<type>(x() * s, y() * s, z() * s, w() * s)); }
		t_vector4t<type> operator / (const type s) const { return (t_vector4t<type>(x() / s, y() / s, z() / s, w() / s)); }

		// t_vector4t<type>& operator += (const type s) { x() += s; y() += s; z() += s; w() += s; return *this; } // UNOFFICIAL
		// t_vector4t<type>& operator -= (const type s) { x() -= s; y() -= s; z() -= s; w() -= s; return *this; } // UNOFFICIAL
		t_vector4t<type>& operator *= (const type s) { x() *= s; y() *= s; z() *= s; w() *= s; return *this; }
		t_vector4t<type>& operator /= (const type s) { x() /= s; y() /= s; z() /= s; w() /= s; return *this; }


		// returns the vector (v * M^T) where v=this is interpreted as a    row-vector
		// same as the vector (M * v  ) where v=this is interpreted as a column-vector
		// M must be manually transposed here (in general, M * v == v * M^T != v^T * M
		// since column-major post-multiplication equals row-major pre-multiplication)
		// UNOFFICIAL! (matrices operate on points)
		t_vector4t<type> operator * (const t_matrix44t<type>& m) const {
			const type _x = x() * m[0] + y() * m[1] + z() * m[ 2];
			const type _y = x() * m[4] + y() * m[5] + z() * m[ 6];
			const type _z = x() * m[8] + y() * m[9] + z() * m[10];
			return (t_vector4t<type>(_x, _y, _z, 0));
		}

		// returns the matrix M = v * w^T where v=this is interpreted as a column-vector
		// and w^T as a row-vector; as such (v^T * w) is a scalar and (v * w^T) a matrix
		// UNOFFICIAL! (matrices operate on points)
		t_matrix44t<type> operator ^ (const t_vector4t<type>& w) const {
			const t_vector4t<type> xv = t_vector4t<type>(x() * w.x(), y() * w.x(), z() * w.x(), type(0));
			const t_vector4t<type> yv = t_vector4t<type>(x() * w.y(), y() * w.y(), z() * w.y(), type(0));
			const t_vector4t<type> zv = t_vector4t<type>(x() * w.z(), y() * w.z(), z() * w.z(), type(0));
			return (t_matrix44t<type>(xv, yv, zv, t_vector4t<type>::w_vector()));
		}


		// cf. point::to_vector
		t_tuple4t<type> to_tuple() const { return (t_tuple4t<type>(*this)); }
		t_point4t<type> to_point() const { assert(w() == type(0)); return (t_point4t<type>::zero_point() + (*this)); }

		// inner product is only defined between vectors and vectors;
		// w-components are normally zero and do not contribute to it
		// (so dedicated "*_xyz" convenience wrappers can be omitted)
		type inner_product(const t_vector4t<type>& v, const t_vector4t<type>& mask = xyz_vector()) const { return (inner(v, mask)); }
		type inner        (const t_vector4t<type>& v, const t_vector4t<type>& mask = xyz_vector()) const {
			return ((x() * v.x() * mask.x()) + (y() * v.y() * mask.y()) + (z() * v.z() * mask.z()) + (w() * v.w() * mask.w()));
		}

		// scalar-triple (mixed) product
		type triple_product(const t_vector4t<type>& v, t_vector4t<type>& w) const { return (triple(v, w)); }
		type triple        (const t_vector4t<type>& v, t_vector4t<type>& w) const { return ((v.outer(w)).inner(*this)); }

		type sq_magnit(const t_vector4t<type>& mask = xyz_vector()) const { return (sq_len(mask)); }
		type    magnit(const t_vector4t<type>& mask = xyz_vector()) const { return (   len(mask)); }
		type sq_length(const t_vector4t<type>& mask = xyz_vector()) const { return (sq_len(mask)); }
		type    length(const t_vector4t<type>& mask = xyz_vector()) const { return (   len(mask)); }
		type sq_len   (const t_vector4t<type>& mask = xyz_vector()) const { return (inner(*this, mask)); }
		type    len   (const t_vector4t<type>& mask = xyz_vector()) const { return (lib_math::sqrt(sq_len(mask))); }

		#if 0
		type dot_xy(const t_vector4t<type>& v) const { return (inner(v, xy_vector())); }
		type dot_yz(const t_vector4t<type>& v) const { return (inner(v, yz_vector())); }
		type dot_xz(const t_vector4t<type>& v) const { return (inner(v, xz_vector())); }

		type sq_magnit_xy() const { return (inner(*this, xy_vector())); }
		type sq_magnit_yz() const { return (inner(*this, yz_vector())); }
		type sq_magnit_xz() const { return (inner(*this, xz_vector())); }

		type magnit_xy() const { return (lib_math::sqrt(sq_magnit_xy())); }
		type magnit_yz() const { return (lib_math::sqrt(sq_magnit_yz())); }
		type magnit_xz() const { return (lib_math::sqrt(sq_magnit_xz())); }
		#endif


		type delta_angle_xy(const t_vector4t<type>& w) const {
			const t_vector4t<type>& v = *this;
			const t_vector4t<type>  c = v.outer_product(w);

			// angle to rotate <this> (anti-clockwise) to <w>, in [-PI, PI>
			// both vectors are assumed not to have a z-component (i.e. lie
			// in the xy-plane) but can be non-normalized
			const type angle = std::acos(inner_product(w) / (v.length() * w.length()));
			const type zsign = lib_math::sign(c.z());

			return (angle * zsign);
		}


		// returns length, normalizes at same time
		// no point adding a const-version of this
		// ("length" == "magnitude")
		type magnit_normalize_ref(const t_vector4t<type>& mask = xyz_vector(), const type eps = M_FEPS) {
			type sq_mag = sq_magnit(mask);
			type    mag = sq_mag;

			if (normalizable(mask, eps * eps)) {
				(*this) *= mask;
				(*this) /= (mag = lib_math::sqrt(sq_mag));
			}

			return mag;
		}


		// note: these are EXTRINSIC rotation operators
		// they rotate <this> around the WORLD'S (fixed)
		// coordinate axes (ie. rot_y_ in xz-plane etc.)
		t_vector4t<type>& rotate_x_ext_ref(type angle) {
			const type ca = std::cos(angle);
			const type sa = std::sin(angle);
			const type ry = ca * y() - sa * z();
			const type rz = sa * y() + ca * z();

			// counter-clockwise rotation (+y becomes +z) for positive angles
			return ((*this) = std::move(t_vector4t<type>(x(), ry, rz)));
		}
		t_vector4t<type>& rotate_y_ext_ref(type angle) {
			const type ca = std::cos(angle);
			const type sa = std::sin(angle);
			const type rx =  ca * x() + sa * z();
			const type rz = -sa * x() + ca * z();

			// counter-clockwise rotation (+z becomes +x) for positive angles
			// +x is (1, 0, 0) and points RIGHT (seen from above)
			// +z is (0, 0, 1) and points DOWN (seen from above)
			return ((*this) = std::move(t_vector4t<type>(rx, y(), rz)));
		}
		t_vector4t<type>& rotate_z_ext_ref(type angle) {
			const type ca = std::cos(angle);
			const type sa = std::sin(angle);
			const type rx = ca * x() - sa * y();
			const type ry = sa * x() + ca * y();

			// counter-clockwise rotation (+x becomes +y) for positive angles
			return ((*this) = std::move(t_vector4t<type>(rx, ry, z())));
		}

		t_vector4t<type>& rotate_ypr_ref(const t_vector4t<type>& ypr) {
			x() = std::cos(ypr.x()) * std::cos(ypr.y());
			y() = std::sin(ypr.y());
			z() = std::sin(ypr.x()) * std::cos(ypr.y());
			return *this;
		}


		// project the vector b onto a=this; a should be normalized
		t_vector4t<type> project(const t_vector4t<type>& b) const {
			const t_vector4t<type>& a = *this;
			const t_vector4t<type>  p = a * a.inner_product(b);
			return p;
		}

		// reflect this vector about the vector <N>
		t_vector4t<type> reflect(const t_vector4t<type>& n) const {
			const t_vector4t<type>& v = *this;
			const t_vector4t<type>  r = n * (n.inner(v) * type(2));
			return (v - r);
		}

		// refract this vector about the vector <N>
		//   n1 is the refr. index of the material being exited
		//   n2 is the refr. index of the material being entered
		t_vector4t<type> refract(const t_vector4t<type>& n, type r1, type r2) const {
			const t_vector4t<type>& v = *this;

			const type r  = r1 / r2;
			const type ct = -(n.inner(v));
			const type st = type(1) - (r * r * (type(1) - ct * ct));

			if (st > type(0))
				return ((v * r) + (n * (r * ct - lib_math::sqrt(st))));

			// total internal reflection
			return n;
		}


		// NOTE: incorrect use of mask here
		t_vector4t<type> outer_product(const t_vector4t<type>& vec, const t_vector4t<type>& mask = xyz_vector()) const { return (outer(vec, mask)); }
		t_vector4t<type> outer        (const t_vector4t<type>& vec, const t_vector4t<type>& mask = xyz_vector()) const {
			const type _x = ((y() * vec.z()) - (z() * vec.y())) * mask.x();
			const type _y = ((z() * vec.x()) - (x() * vec.z())) * mask.y();
			const type _z = ((x() * vec.y()) - (y() * vec.x())) * mask.z();
			return (t_vector4t<type>(_x, _y, _z, type(0)));
		}

		// note: these divide instead of multiplying by
		// the reciprocal of length (<type> can be int)
		t_vector4t<type> normalize(const t_vector4t<type>& mask = xyz_vector(), const type eps = M_FEPS) const {
			t_vector4t<type> r = *this;

			if (normalizable(mask, eps * eps))
				return ((r * mask) / lib_math::sqrt(sq_magnit(mask)));

			return r;
		}

		t_vector4t<type>& normalize_ref() { return ((*this) = normalize()); }
		t_vector4t<type>& outer_ref(const t_vector4t<type>& v) { return ((*this) = outer(v)); }
		t_vector4t<type>& outer_product_ref(const t_vector4t<type>& v) { return (outer_ref(v)); }



		#if 0
		// these are redundant with multiplicative-masks
		t_vector4t<type> xz() const { return (t_vector4t<type>(    x(), type(0),     z(), w())); }
		t_vector4t<type> xy() const { return (t_vector4t<type>(    x(),     y(), type(0), w())); }
		t_vector4t<type> yz() const { return (t_vector4t<type>(type(0),     y(),     z(), w())); }
		#endif


		#if 0
		t_vector4t<type> project_xz() const { return (project(xz_vector())); }
		t_vector4t<type> project_xy() const { return (project(xy_vector())); }
		t_vector4t<type> project_yz() const { return (project(yz_vector())); }

		t_vector4t<type>& project_xz_ref() { return ((*this) = project_xz()); }
		t_vector4t<type>& project_xy_ref() { return ((*this) = project_xy()); }
		t_vector4t<type>& project_yz_ref() { return ((*this) = project_yz()); }
		#endif

		t_vector4t<type> abs(                         ) const { return (t_vector4t<type>(std::abs(      x()),        std::abs(     y()),         std::abs(      z()),        w())); }
		t_vector4t<type> max(const t_vector4t<type>& v) const { return (t_vector4t<type>(std::max<type>(x(), v.x()), std::max<type>(y(), v.y()), std::max<type>(z(), v.z()), w())); }
		t_vector4t<type> min(const t_vector4t<type>& v) const { return (t_vector4t<type>(std::min<type>(x(), v.x()), std::min<type>(y(), v.y()), std::min<type>(z(), v.z()), w())); }
		t_vector4t<type> pow(const t_vector4t<type>& v) const { return (t_vector4t<type>(std::pow<type>(x(), v.x()), std::pow<type>(y(), v.y()), std::pow<type>(z(), v.z()), w())); }
		t_vector4t<type> mod(const t_vector4t<type>& v) const { return (t_vector4t<type>(std::fmod     (x(), v.x()), std::fmod     (y(), v.y()), std::fmod     (z(), v.z()), w())); }
		// t_vector4t<type> mod(const t_vector4t<type>& v) const { return (t_vector4t<type>(std::mod<type>(x(), v.x()), std::mod<type>(y(), v.y()), std::mod<type>(z(), v.z()), w())); }

		t_vector4t<type> signs() const { return (t_vector4t<type>(signum(x()), signum(y()), signum(z()), signum(w()))); }
		t_vector4t<type> clamp(const t_vector4t<type>& vmin, const t_vector4t<type>& vmax) const {
		t_vector4t<type> r;
			r.x() = lib_math::clamp(x(), vmin.x(), vmax.x());
			r.y() = lib_math::clamp(y(), vmin.y(), vmax.y());
			r.z() = lib_math::clamp(z(), vmin.z(), vmax.z());
			r.w() = lib_math::clamp(w(), vmin.w(), vmax.w());
			return r;
		}


		template<typename t_rng> t_vector4t<type>  randomize(t_rng& rng) const {
			return (t_vector4t<type>((rng() * type(2)) - type(1), (rng() * type(2)) - type(1), (rng() * type(2)) - type(1)));
		}
		template<typename t_rng> t_vector4t<type> urandomize(t_rng& rng, const type eps = M_FEPS) const {
			t_vector4t<type> r = xyz_vector();

			while (std::abs(r.sq_len() - type(1)) > eps) {
				r.x() = (rng() * type(2)) - type(1);
				r.y() = (rng() * type(2)) - type(1);
				r.z() = (rng() * type(2)) - type(1);
			}

			return r;
		}

		template<typename t_rng> t_vector4t<type>&  randomize_ref(t_rng& rng                         ) { return ((*this) =  randomize(   )); }
		template<typename t_rng> t_vector4t<type>& urandomize_ref(t_rng& rng, const type eps = M_FEPS) { return ((*this) = urandomize(eps)); }


		uint32_t hash(const t_vector4t<type>& mask = ones_vector()) const {
			const t_tuple4t<type> t(*this);
			const t_tuple4t<type> m( mask);
			return (t.hash(m));
		}

		bool equals(const t_vector4t<type>& vec, const t_vector4t<type>& eps = eps_vector()) const {
			const t_tuple4t<type> t(*this);
			const t_tuple4t<type> v(  vec);
			const t_tuple4t<type> e(  eps);
			return (t.equals(v, e));
		}
		// vector is normalizable if its (masked) magnitude != {0 or 1}
		bool normalizable(const t_vector4t<type>& mask, const type eps = M_FEPS) const {
			return (sq_magnit(mask) > eps && std::abs(type(1) - sq_magnit(mask)) > eps);
		}

		#if 0
		bool normalizable_xz(const type eps = M_FEPS) const { return (normalizable(xz_vector(), eps * eps)); }
		bool normalizable_xy(const type eps = M_FEPS) const { return (normalizable(xy_vector(), eps * eps)); }
		bool normalizable_yz(const type eps = M_FEPS) const { return (normalizable(yz_vector(), eps * eps)); }
		#endif


		void print() const { (t_tuple4t<type>(*this)).print(); }
		void sassert() const { (t_tuple4t<type>(*this)).sassert(); }


		const type* xyzw(unsigned int idx = 0) const { return &m_xyzw[idx]; }
		const type* data(unsigned int idx = 0) const { return &m_xyzw[idx]; }
		      type* xyzw(unsigned int idx = 0)       { return &m_xyzw[idx]; }
		      type* data(unsigned int idx = 0)       { return &m_xyzw[idx]; }

		type  operator [] (uint32_t n) const { return m_xyzw[n]; }
		type& operator [] (uint32_t n)       { return m_xyzw[n]; }

		type  x() const { return m_xyzw[0]; }
		type  y() const { return m_xyzw[1]; }
		type  z() const { return m_xyzw[2]; }
		type  w() const { return m_xyzw[3]; }
		type& x()       { return m_xyzw[0]; }
		type& y()       { return m_xyzw[1]; }
		type& z()       { return m_xyzw[2]; }
		type& w()       { return m_xyzw[3]; }


		static const t_vector4t<type>&    pi_vector();
		static const t_vector4t<type>&   eps_vector();
		static const t_vector4t<type>&   inf_vector();
		static const t_vector4t<type>&  zero_vector();
		static const t_vector4t<type>&  ones_vector();
		static const t_vector4t<type>& error_vector();
		static const t_vector4t<type>&     x_vector();
		static const t_vector4t<type>&     y_vector();
		static const t_vector4t<type>&     z_vector();
		static const t_vector4t<type>&     w_vector();
		static const t_vector4t<type>&    xz_vector();
		static const t_vector4t<type>&    xy_vector();
		static const t_vector4t<type>&    yz_vector();
		static const t_vector4t<type>&   xyz_vector();
		static const t_vector4t<type>&  xyzw_vector();

	private:
		type m_xyzw[4];
	};


	typedef t_vector4t< int32_t> t_vector4i;
	typedef t_vector4t<uint32_t> t_vector4ui;

	typedef t_vector4t< int64_t> t_vector4li;
	typedef t_vector4t<uint64_t> t_vector4lui;

	typedef t_vector4t< float> t_vector4f;
	typedef t_vector4t<double> t_vector4d;

	// aliases for less typing
	typedef t_vector4i  t_vec4i;
	typedef t_vector4ui t_vec4ui;

	typedef t_vector4li  t_vec4li;
	typedef t_vector4lui t_vec4lui;

	typedef t_vector4f t_vec4f;
	typedef t_vector4d t_vec4d;
};

#endif

