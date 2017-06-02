#ifndef LIBGEOM_TEMPLATE_FUNCS_HDR
#define LIBGEOM_TEMPLATE_FUNCS_HDR

#include <xmmintrin.h>

#include <algorithm>
#include <cmath>

namespace lib_sys {
	static uint8_t is_little_endian_arch() {
		const int16_t i = 1;
		const int8_t* b = reinterpret_cast<const int8_t*>(&i);
		return (b[0] == 1);
	}
};

namespace lib_math {
	static const uint8_t IS_LITTLE_ENDIAN = lib_sys::is_little_endian_arch();

	static constexpr size_t POWERS_OF_TWO[] = {1,  2,   4,    8,    16,     32,      64,      128};
	static constexpr size_t POWERS_OF_TEN[] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000};

	template<typename type> struct t_matrix44t;
	template<typename type> struct t_vector4t;

	// note: these three are already in gcc's cmath
	template<typename type> type copysign_builtin(type v, type w) { return (__builtin_copysignf(v, w)); }
	template<typename type> type     sqrt_builtin(type v) { return (          __builtin_sqrtf(v)); }
	template<typename type> type inv_sqrt_builtin(type v) { return (type(1) / __builtin_sqrtf(v)); }

	#if 0
	// not needed on x64 architectures
	#define attr_align_stack __attribute__ ((force_align_arg_pointer))
	#else
	#define attr_align_stack
	#endif
	template<typename type> type  attr_align_stack      sqrt_sse(type v) {
		__m128 vec = _mm_set_ss(v); vec = _mm_sqrt_ss(vec); return (_mm_cvtss_f32(vec));
	}
	template<typename type> type  attr_align_stack  inv_sqrt_sse(type v) {
		__m128 vec = _mm_set_ss(v); vec = _mm_rsqrt_ss(vec); return (_mm_cvtss_f32(vec));
	}

	template<typename float4>  attr_align_stack  __m128 float4_to_m128(const float4& f) { return (_mm_loadu_ps(&f[0])); }
	template<typename float4>  attr_align_stack  float4 m128_to_float4(const __m128& m) {
		union { __m128 m; float f[4]; } u;
		u.m = m;
		return (float4(u.f[0], u.f[1], u.f[2], u.f[3]));
	}


	template<typename type> type     sqrt(type v) { return (               std::sqrt(v)); }
	template<typename type> type inv_sqrt(type v) { return (type(1) / lib_math::sqrt(v)); }

	template<typename type> type min3(type a, type b, type c) { return (std::min<type>(a, std::min<type>(b, c))); }
	template<typename type> type max3(type a, type b, type c) { return (std::max<type>(a, std::max<type>(b, c))); }
	template<typename type> type sign(type v) { return ((v >= type(0)) * type(2) - type(1)); }

	template<typename type> int8_t ssign(const type v) { return ((v > type(0)) - (v < type(0))); }
	template<typename type> int8_t bsign(const type v) {
		static_assert(std::is_signed<type>::value, "");

		#if 0
		switch (sizeof(type)) {
			// breaks strict-aliasing rules
			case 1: { return (((*reinterpret_cast<const int8_t *>(&v) >>  7) & 1) * -2 + 1); } break;
			case 2: { return (((*reinterpret_cast<const int16_t*>(&v) >> 15) & 1) * -2 + 1); } break;
			case 4: { return (((*reinterpret_cast<const int32_t*>(&v) >> 31) & 1) * -2 + 1); } break;
			case 8: { return (((*reinterpret_cast<const int64_t*>(&v) >> 63) & 1) * -2 + 1); } break;
			default: {} break;
		}
		return 0;
		#else
		// sensitive to endian-ness
		constexpr size_t byte_idx = sizeof(type) - 1;

		const  int8_t sign_bit = (reinterpret_cast<const int8_t*>(&v)[byte_idx * IS_LITTLE_ENDIAN] >> 7) & 1;
		const  int8_t sign_val = sign_bit * -2 + 1;
		return sign_val;
		#endif
	}

	// (type(1) - wgt) works only if <type> is a scalar, so we rewrite it to a form accepting
	// both scalars and vectors (X * (1 - a) + Y * a == X - X * a + Y * a == X - (X + Y) * a)
	// the latter however still requires wgt to be specified as a vector which is painful
	// template<typename type> type lerp(type v0, type v1, type wgt) { return (v0 * (type(1) - wgt) + v1 * wgt); }
	// template<typename type> type lerp(type v0, type v1, type wgt) { return (v0 - (v0 + v1) * wgt); }
	template<typename type> type lerp(type v0, type v1, float wgt) { return (v0 * (1.0f - wgt) + v1 * wgt); }
	template<typename type> type scale(type v, type v0, type v1) { return ((v - v0) / (v1 - v0)); }
	// convert a value x on scale [a,b] to a value y on scale [c,d]
	template<typename type> type rescale(type x, type a, type b, type c, type d) { return (c + (d - c) * scale(x, a, b)); }
	template<typename type> type clamp(type v, type v0, type v1) { return (std::max<type>(v0, std::min<type>(v1, v))); }
	template<typename type> type square(type v) { return (v * v); }

	template<typename type> type abss(type v) { return (lib_math::sign(v) * v); }
	template<typename type> type absl(type v) { return (lib_math::lerp(v, -v, v < type(0))); }
	template<typename type> type absm(type v) { return (std::max(v, -v)); }


	template<typename type> t_vector4t<type> nlerp(const t_vector4t<type>& v0, const t_vector4t<type>& v1, type wgt) {
		assert(v0 != -v1);
		const t_vector4t<type>& vi = lib_math::lerp(v0, v1, wgt);
		const t_vector4t<type>  vn = vi.normalize();
		return vn;
	}
	template<typename type> t_vector4t<type> slerp(const t_vector4t<type>& v0, const t_vector4t<type>& v1, type wgt, type eps) {
		// calculate the angle subtended by the arc spanned by unit-vectors
		// v0 and v1 (angle=acos(dot(v0, v1))); the interpolated vector (vi)
		// is a linear combination of the orthonormal basis (v0, v2) where
		//   v2 = normalize(v1 - v0 * dot(v0, v1))
		//   vi = v0 * cos(angle * wgt) + v2 * sin(angle * wgt)
		const type  cosa = lib_math::clamp(v0.inner(v1), type(-1), type(1));
		const type acosa = std::acos(cosa);
		const type isina = type(1) / std::sin(acosa);

		// with an added normalize, this reduces to the general nlerp
		if (std::fabs(cosa) > (type(1) - eps))
			return (lib_math::lerp(v0, v1, wgt));

		// note: explicit normalization is more accurate
		// const t_vector4t<type> v2 = (v1 - v0 * cosa) * isina;
		// const t_vector4t<type> vi = v0 * std::cos(acosa * wgt) + v2 * std::sin(acosa * wgt);
		// return vi;
		// standard formulation
		const t_vector4t<type> w0 = v0 * (std::sin((type(1) - wgt) * acosa) * isina);
		const t_vector4t<type> w1 = v1 * (std::sin((type(0) + wgt) * acosa) * isina);
		return (w0 + w1);
	}

	#if 0
	template<typename type> t_vector4t<type> slerp_aux(
		const t_vector4t<type>& v,
		const t_vector4t<type>& w,
		const type wgt,
		const type eps
	) {
		const type cos_angle = lib_math::clamp(v.inner(w), type(-1), type(1));
		const type ccw_slerp = lib_math::sign(v.inner(w.outer(t_vector4t<type>::y_vector())));
		const type angle_max = std::acos(cos_angle); // radians
		const type angle_int = angle_max * wgt * ccw_slerp;

		// if dot(v, w) is  1, no interpolation is required
		// if dot(v, w) is -1, no interpolation is possible
		if (cos_angle <= (type(-1) + eps))
			return t_vector4t<type>::zero_vector();
		if (cos_angle >= (type(1) - eps))
			return v;

		// rotate in world xz-plane, let caller transform back
		return (v.rotate_y_ext(angle_int));
	}

	template<typename type> t_vector4t<type> slerp(const t_vector4t<type>& vz, const t_vector4t<type>& vw, type wgt, type eps) {
		// Nth-order polynomial vector interpolation over angles
		//
		// two (non-colinear) vectors v and w uniquely span a plane with
		// normal (v cross w); we could rotate directly in this plane to
		// do interpolation but take the simple approach of constructing
		// an axis-system around N and making this align with the world
		// axis-system by transposing it, s.t. rotating v around N to w
		// equals rotating v' around world-y to w'
		const t_vector4t<type> vy = (vz.outer(vw)).normalize();
		const t_vector4t<type> vx = (vz.outer(vy)).normalize();

		// compose axis-system; vy is plane normal
		const t_matrix44t<type> m1 = t_matrix44t<type>(vx, vy, vz);
		const t_matrix44t<type> m2 = m1.transpose_rotation();

		// inv-transform source and target to axis-aligned vectors
		// interpolate in the local space, transform back to world
		return (m1 * slerp_aux(m2 * vz, m2 * vw, wgt, eps));
	}
	#endif


	template<typename type> t_vector4t<type> angle_to_vector_xy(type angle) { return (t_vector4t<type>(std::cos(angle), std::sin(angle),         type(0))); }
	template<typename type> t_vector4t<type> angle_to_vector_xz(type angle) { return (t_vector4t<type>(std::cos(angle),         type(0), std::sin(angle))); }
	template<typename type> t_vector4t<type> angle_to_vector_yz(type angle) { return (t_vector4t<type>(        type(0), std::cos(angle), std::sin(angle))); }
	template<typename type> type vector_to_angle_xy(const t_vector4t<type>& vector) { return (std::atan2(vector.y(), vector.x())); }
	template<typename type> type vector_to_angle_xz(const t_vector4t<type>& vector) { return (std::atan2(vector.z(), vector.x())); }
	template<typename type> type vector_to_angle_yz(const t_vector4t<type>& vector) { return (std::atan2(vector.z(), vector.y())); }

	template<typename type> constexpr type deg_to_rad() { return (M_PI / type(180)); }
	template<typename type> constexpr type rad_to_deg() { return (type(180) / M_PI); }

	template<typename type> type angle_lerp(type v0, type v1, type wgt, type eps) {
		// NOTE:
		//   angles can not be continuously defined, so this approach is semi-valid
		//   (e.g. atan2 has a transition at PI/-PI which breaks lerp, transforming
		//   by (a + 2PI) % 2PI only shifts the discontinuity to 2PI/0)
		//   comparing signs is only valid for angles in [-PI, PI], not [0, 2PI]
		if (lib_math::sign(v0) == lib_math::sign(v1))
			return (lib_math::lerp(v0, v1, wgt));
		const t_vector4t<type>& vv0 = lib_math::angle_to_vector_xz(v0);
		const t_vector4t<type>& vv1 = lib_math::angle_to_vector_xz(v1);
		const t_vector4t<type>& vvi = lib_math::slerp(vv0, vv1, wgt, eps);
		return (lib_math::vector_to_angle_xz(vvi));
	}

	// clamp an angle in radians to the range [0, 2PI]
	template<typename type> type clamp_angle_rad(type raw_angle) {
		constexpr float max_angle = M_PI + M_PI;
		raw_angle = std::fmod(raw_angle, max_angle);
		raw_angle += (max_angle * (raw_angle < type(0)));
		return raw_angle;
	}
	// clamp an angle in degrees to the range [0, 360]
	template<typename type> type clamp_angle_deg(type raw_angle) {
		return (lib_math::clamp_angle_rad(raw_angle * lib_math::deg_to_rad<type>()) * lib_math::rad_to_deg<type>());
	}

	// behaves as relative-tolerance comparison when abs(a) and abs(b)
	// are both > 1, otherwise behaves as absolute-tolerance comparison
	// (alternatively use "abs(a - b) <= (eps * (1 + abs(a) + abs(b)))")
	//
	template<typename type> bool fp_eq(type a, type b, type eps) {
		return (a == b || std::abs(a - b) <= (eps * max3<type>(std::abs(a), std::abs(b), type(1))));
	}

	// unlike sign(), treats zero specially
	#if 0
	template<typename type> type signum(type v) {
		if (v > type(0)) return (type( 1));
		if (v < type(0)) return (type(-1));
		return (type(0));
	}
	#else
	template<typename type> type signum(type v) {
		const type gtz = (v >  type(0)) * type(2) - type(1);
		const type eqz = (v == type(0)) * type(1);
		return (gtz + eqz);
	}
	#endif

	template<typename type> type round(type v, size_t n) {
		if (n > 0) {
			// round number to <n> decimals
			const int i = std::min(7, int(sizeof(POWERS_OF_TEN) / sizeof(POWERS_OF_TEN[0])) - 1);
			const int n = lib_math::clamp(n, 0, i);

			const type vinteg = std::floor(v);
			const type vfract = v - vinteg;

			return (vinteg + std::floor((vfract * POWERS_OF_TEN[n]) + type(0.5f)) / POWERS_OF_TEN[n]);
		}

		return (std::floor(v + type(0.5f)));
	}
};

#endif

