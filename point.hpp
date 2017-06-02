#ifndef LIBGEOM_POINT_HDR
#define LIBGEOM_POINT_HDR

#include <cassert>
#include <cstdint>
#include <utility>

#include "./math_defs.hpp"

namespace lib_math {
	// we do *not* want these to inherit from t_tuple
	// for type-safety reasons: (point - point) should
	// always be a vector, etc.
	template<typename type> struct t_tuple4t;
	template<typename type> struct t_vector4t;

	template<typename type> struct t_point4t {
	public:
		t_point4t<type>(const type _x = type(0), const type _y = type(0), const type _z = type(0), const type _w = type(1)) {
			x() = _x; y() = _y;
			z() = _z; w() = _w;
		}
		t_point4t<type>(const type* v) { *this = std::move(t_point4t<type>(v[0], v[1], v[2], v[3])); }
		t_point4t<type>(const t_tuple4t<type>& t) { *this = std::move(t_point4t<type>(t.x(), t.y(), t.z(), t.w())); }
		t_point4t<type>(const t_point4t<type>& p) { *this = p; }

		bool operator == (const t_point4t<type>& p) const { return ( equals(p)); }
		bool operator != (const t_point4t<type>& p) const { return (!equals(p)); }

		// rules:
		//   type(point + vector) := type(point)
		//   type(point - vector) := type(point + -vector)
		//   type(point - point) := type(vector)
		//   type(point + point) := UNDEFINED
		//   points: w=1, vectors: w=0
		t_point4t<type> operator + (const t_vector4t<type>& v) const { return (t_point4t<type>(x() + v.x(), y() + v.y(), z() + v.z(), w() + v.w())); }
		t_point4t<type> operator - (const t_vector4t<type>& v) const { return (t_point4t<type>(x() - v.x(), y() - v.y(), z() - v.z(), w() - v.w())); }
		t_vector4t<type> operator - (const t_point4t<type>& p) const { return (t_vector4t<type>(x() - p.x(), y() - p.y(), z() - p.z(), w() - p.w())); }

		t_point4t<type>& operator  = (const t_point4t<type>& p)  { x()  = p.x(); y()  = p.y(); z()  = p.z(); w()  = p.w(); return *this; }
		t_point4t<type>& operator += (const t_vector4t<type>& v) { x() += v.x(); y() += v.y(); z() += v.z(); w() += v.w(); return *this; }
		t_point4t<type>& operator -= (const t_vector4t<type>& v) { x() -= v.x(); y() -= v.y(); z() -= v.z(); w() -= v.w(); return *this; }

		// t_point4t<type> operator - () const { return (t_point4t<type>(-x(), -y(), -z(), -w())); }

		// cf. vector::to_point
		t_tuple4t<type> to_tuple() const { return (t_tuple4t<type>(*this)); }
		t_vector4t<type> to_vector() const { assert(w() == type(1)); return ((*this) - zero_point()); }
		t_point4t<type> to_ndc() const { assert(w() != type(0)); return {x() / w(), y() / w(), z() / w(), type(1)}; }


		uint32_t hash(const t_point4t<type>& mask = ones_point()) const {
			const t_tuple4t<type> t(*this);
			const t_tuple4t<type> m( mask);
			return (t.hash(m));
		}

		bool equals(const t_point4t<type>& pnt, const t_point4t<type>& eps = eps_point()) const {
			const t_tuple4t<type> t(*this);
			const t_tuple4t<type> p(  pnt);
			const t_tuple4t<type> e(  eps);
			return (t.equals(p, e));
		}


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


		static const t_point4t<type>&  eps_point();
		static const t_point4t<type>& null_point();
		static const t_point4t<type>& zero_point();
		static const t_point4t<type>& ones_point();

		static type eps_scalar() { return (eps_point.x()); }

	private:
		type m_xyzw[4];
	};


	typedef t_point4t< int32_t> t_point4i;
	typedef t_point4t<uint32_t> t_point4ui;

	typedef t_point4t< int64_t> t_point4li;
	typedef t_point4t<uint64_t> t_point4lui;

	typedef t_point4t< float> t_point4f;
	typedef t_point4t<double> t_point4d;

	// aliases for less typing
	typedef t_point4i  t_pos4i;
	typedef t_point4ui t_pos4ui;

	typedef t_point4li  t_pos4li;
	typedef t_point4lui t_pos4lui;

	typedef t_point4f t_pos4f;
	typedef t_point4d t_pos4d;
};

#endif

