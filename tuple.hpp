#ifndef LIBGEOM_TUPLE_HDR
#define LIBGEOM_TUPLE_HDR

#include <cassert>
#include <cmath>
#include <cstdint>

#include "./math_defs.hpp"
#include "./template_funcs.hpp"

namespace lib_math {
	template<typename type> struct t_point4t;
	template<typename type> struct t_vector4t;

	template<typename type> struct t_tuple4t {
	public:
		t_tuple4t<type>(const type _x = type(0), const type _y = type(0), const type _z = type(0), const type _w = type(0)) {
			x() = _x; y() = _y;
			z() = _z; w() = _w;
		}
		// safer than const type*
		t_tuple4t<type>(const type* v) { *this = std::move(t_tuple4t<type>(v[0], v[1], v[2], v[3])); }
		t_tuple4t<type>(const t_point4t<type>& p) { *this = std::move(t_tuple4t<type>(p.x(), p.y(), p.z(), p.w())); }
		t_tuple4t<type>(const t_vector4t<type>& v) { *this = std::move(t_tuple4t<type>(v.x(), v.y(), v.z(), v.w())); }
		t_tuple4t<type>(const t_tuple4t<type>& t) { *this = t; }

		bool operator == (const t_tuple4t<type>& t) const { return ( equals(t)); }
		bool operator != (const t_tuple4t<type>& t) const { return (!equals(t)); }

		t_tuple4t<type>& operator = (const t_tuple4t<type>& t) {
			x() = t.x(); y() = t.y();
			z() = t.z(); w() = t.w();
			return *this;
		}

		t_tuple4t<type> operator + (const t_tuple4t<type>& t) const { return (t_tuple4t<type>(x() + t.x(), y() + t.y(), z() + t.z(), w() + t.w())); }
		t_tuple4t<type> operator - (const t_tuple4t<type>& t) const { return (t_tuple4t<type>(x() - t.x(), y() - t.y(), z() - t.z(), w() - t.w())); }

		t_tuple4t<type> operator * (const type s) const { return (t_tuple4t<type>(x() * s, y() * s, z() * s, w() * s)); }
		t_tuple4t<type> operator / (const type s) const { return (t_tuple4t<type>(x() / s, y() / s, z() / s, w() / s)); }


		uint32_t hash(const t_tuple4t<type>& mask = ones_tuple()) const {
			const uint32_t hx = (*reinterpret_cast<uint32_t*>(const_cast<type*>(&m_xyzw[0]))) * mask.x();
			const uint32_t hy = (*reinterpret_cast<uint32_t*>(const_cast<type*>(&m_xyzw[1]))) * mask.y();
			const uint32_t hz = (*reinterpret_cast<uint32_t*>(const_cast<type*>(&m_xyzw[2]))) * mask.z();
			const uint32_t hw = (*reinterpret_cast<uint32_t*>(const_cast<type*>(&m_xyzw[3]))) * mask.w();
			return (hx ^ hy ^ hz ^ hw);
		}

		static uint32_t hash(const t_point4t<type>& pnt, const t_point4t<type>& mask) {
			const t_tuple4t<type> p( pnt);
			const t_tuple4t<type> m(mask);
			return (p.hash(m));
		}
		static uint32_t hash(const t_vector4t<type>& vec, const t_vector4t<type>& mask) {
			const t_tuple4t<type> v( vec);
			const t_tuple4t<type> m(mask);
			return (v.hash(m));
		}


		bool equals(const t_tuple4t<type>& t, const t_tuple4t<type>& eps = eps_tuple()) const {
			unsigned int mask = 0;
			mask += lib_math::fp_eq<type>(x(), t.x(), eps.x());
			mask += lib_math::fp_eq<type>(y(), t.y(), eps.y());
			mask += lib_math::fp_eq<type>(z(), t.z(), eps.z());
			mask += lib_math::fp_eq<type>(w(), t.w(), eps.w());
			return (mask == 4);
		}

		bool is_point() const { return (w() == type(1)); }
		bool is_vector() const { return (w() == type(0)); }

		static bool equals(const t_point4t<type>& pnt_a, const t_point4t<type>& pnt_b, const t_point4t<type>& eps) {
			return (equals(t_tuple4t<type>(pnt_a), t_tuple4t<type>(pnt_b), t_tuple4t<type>(eps)));
		}
		static bool equals(const t_vector4t<type>& vec_a, const t_vector4t<type>& vec_b, const t_vector4t<type>& eps) {
			return (equals(t_tuple4t<type>(vec_a), t_tuple4t<type>(vec_b), t_tuple4t<type>(eps)));
		}


		void print() const;
		void sassert() const {
			assert(!std::isnan(x()) && !std::isinf(x()));
			assert(!std::isnan(y()) && !std::isinf(y()));
			assert(!std::isnan(z()) && !std::isinf(z()));
			assert(!std::isnan(w()) && !std::isinf(w()));
		}

		const type* xyzw() const { return &m_xyzw[0]; }
		      type* xyzw()       { return &m_xyzw[0]; }

		// no need to make functors from tuples
		// operator const type* () const { return xyzw(); }
		// operator       type* ()       { return xyzw(); }

		type  operator [] (unsigned int n) const { return m_xyzw[n]; }
		type& operator [] (unsigned int n)       { return m_xyzw[n]; }

		// same as operator []
		// type  i(std::size_t n) const { return m_xyzw[n]; }
		// type& i(std::size_t n)       { return m_xyzw[n]; }
		type  x() const { return m_xyzw[0]; }
		type  y() const { return m_xyzw[1]; }
		type  z() const { return m_xyzw[2]; }
		type  w() const { return m_xyzw[3]; }
		type& x()       { return m_xyzw[0]; }
		type& y()       { return m_xyzw[1]; }
		type& z()       { return m_xyzw[2]; }
		type& w()       { return m_xyzw[3]; }

		static const t_tuple4t<type>& zero_tuple();
		static const t_tuple4t<type>& ones_tuple();
		static const t_tuple4t<type>&  eps_tuple();

	private:
		type m_xyzw[4];
	};


	typedef t_tuple4t< int32_t> t_tuple4i;
	typedef t_tuple4t<uint32_t> t_tuple4ui;

	typedef t_tuple4t< int64_t> t_tuple4li;
	typedef t_tuple4t<uint64_t> t_tuple4lui;

	typedef t_tuple4t< float> t_tuple4f;
	typedef t_tuple4t<double> t_tuple4d;

	// aliases for less typing
	typedef t_tuple4i  t_tup4i;
	typedef t_tuple4ui t_tup4ui;

	typedef t_tuple4li  t_tup4li;
	typedef t_tuple4lui t_tup4lui;

	typedef t_tuple4f t_tup4f;
	typedef t_tuple4d t_tup4d;
};

#endif

