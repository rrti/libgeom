#ifndef LIBGEOM_TEMPLATE_FUNCS_HDR
#define LIBGEOM_TEMPLATE_FUNCS_HDR

#include <xmmintrin.h>
#include <cmath>
#include <vector>

namespace lib_math {
	static const size_t POWERS_OF_TWO[] = {1,  2,   4,    8,    16,     32,      64,      128};
	static const size_t POWERS_OF_TEN[] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000};

	template<typename type> type  __attribute__ ((force_align_arg_pointer))  sqrt_sse(type v) {
		__m128 vec = _mm_set_ss(v); vec = _mm_sqrt_ss(vec); return (_mm_cvtss_f32(vec));
	}
	template<typename type> type  __attribute__ ((force_align_arg_pointer))  inv_sqrt_sse(type v) {
		__m128 vec = _mm_set_ss(v); vec = _mm_rsqrt_ss(vec); return (_mm_cvtss_f32(vec));
	}

	template<typename type> type     sqrt(type v) { return (               std::sqrt(v)); }
	template<typename type> type inv_sqrt(type v) { return (type(1) / lib_math::sqrt(v)); }

	template<typename type> type abs(type v) { return (std::max(v, -v)); }
	template<typename type> type min3(type a, type b, type c) { return (std::min<type>(a, std::min<type>(b, c))); }
	template<typename type> type max3(type a, type b, type c) { return (std::max<type>(a, std::max<type>(b, c))); }
	template<typename type> type sign(type v) { return ((v >= type(0)) * type(2) - type(1)); }
	template<typename type> type lerp(type vmin, type vmax, type alpha) { return (vmin * (type(1) - alpha) + vmax * alpha); }
	template<typename type> type norm(type v, type vmin, type vmax) { return ((v - vmin) / (vmax - vmin)); }
	template<typename type> type clamp(type v, type vmin, type vmax) { return (std::max<type>(vmin, std::min<type>(vmax, v))); }
	template<typename type> type square(type v) { return (v * v); }

	// clamp an angle in radians to the range [0, 2PI]
	template<typename type> type clamp_angle_rad(type raw_angle) {
		constexpr float max_angle = M_PI + M_PI;
		raw_angle = std::fmod(raw_angle, max_angle);
		raw_angle += (max_angle * (raw_angle < type(0)));
		return raw_angle;
	}
	// clamp an angle in degrees to the range [0, 360]
	template<typename type> type clamp_angle_deg(type raw_angle) {
		constexpr float deg_to_rad = M_PI / 180.0f;
		constexpr float rad_to_deg = 180.0f / M_PI;
		return (clamp_angle_rad(raw_angle * deg_to_rad) * rad_to_deg);
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
		const int gtz = (v >  type(0)) * 2 - 1;
		const int eqz = (v == type(0)) * 1;
		return (gtz + eqz);
	}
	#endif

	template<typename type> type round(type v, size_t n) {
		if (n > 0) {
			// round number to <n> decimals
			const int i = std::min(7, int(sizeof(POWERS_OF_TEN) / sizeof(POWERS_OF_TEN[0])) - 1);
			const int n = clamp(n, 0, i);

			const type vinteg = std::floor(v);
			const type vfract = v - vinteg;

			return (vinteg + std::floor((vfract * POWERS_OF_TEN[n]) + type(0.5f)) / POWERS_OF_TEN[n]);
		}

		return (std::floor(v + type(0.5f)));
	}


	template<typename type> type vec_sum(const std::vector<type>& v) {
		type r = type();
		for (size_t n = 0; n < v.size(); n++) {
			r += v[n];
		}
		return r;
	}

	template<typename type> type vec_avg(const std::vector<type>& v) {
		return (vec_sum<type>(v) / std::max(1lu, v.size()));
	}
};

#endif

