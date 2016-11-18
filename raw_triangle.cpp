#include <algorithm>
#include <cassert>

#include "./raw_triangle.hpp"
#include "./idx_triangle.hpp"
#include "./matrix.hpp"
#include "./plane.hpp"
#include "./ray.hpp"
#include "./template_funcs.hpp"

namespace lib_math {
	t_raw_triangle::t_raw_triangle() {
		m_verts[0] = m_point_type::zero_point();
		m_verts[1] = m_point_type::zero_point();
		m_verts[2] = m_point_type::zero_point();
		m_normal   = m_vector_type::zero_vector();
		m_radius   = m_coor_type(0);
	}
	t_raw_triangle::t_raw_triangle(const m_point_type& vert_a, const m_point_type& vert_b, const m_point_type& vert_c) {
		// if looking at triangle from the side in which
		// vertices are ordered CCW, normal points along
		// viewing direction (away from eye)
		// *this = t_raw_triangle(vert_a, vert_b, vert_c, ((vert_c - vert_a).normalize()).outer((vert_b - vert_a).normalize()).normalize());
		*this = t_raw_triangle(vert_a, vert_b, vert_c, ((vert_c - vert_a).outer(vert_b - vert_a)).normalize());
	}
	t_raw_triangle::t_raw_triangle(const m_point_type& vert_a, const m_point_type& vert_b, const m_point_type& vert_c, const m_vector_type& n) {
		m_verts[0] = vert_a;
		m_verts[1] = vert_b;
		m_verts[2] = vert_c;
		m_normal   = n * m_vector_type::xyz_vector();

		assert(m_normal.w() == m_coor_type(0));

		if (n.w() <= m_coor_type(0)) {
			// regular normal
			m_radius = calc_radius();
		} else {
			// extended normal
			m_radius = n.w();
		}
	}

	t_raw_triangle& t_raw_triangle::operator = (const t_raw_triangle& triangle) {
		m_verts[0] = triangle.get_vertex(0);
		m_verts[1] = triangle.get_vertex(1);
		m_verts[2] = triangle.get_vertex(2);
		m_normal   = triangle.get_normal();
		m_radius   = triangle.get_radius();
		return *this;
	}

	t_raw_triangle& t_raw_triangle::invert_winding() {
		const m_point_type p = m_verts[0];

		m_verts[0] = m_verts[2];
		m_verts[2] = p;
		m_normal   = -m_normal;
		return *this;
	}

	t_raw_triangle t_raw_triangle::transform(const m_matrix_type& matrix) const {
		t_raw_triangle r;
		r.get_vertex(0) = matrix * get_vertex(0);
		r.get_vertex(1) = matrix * get_vertex(1);
		r.get_vertex(2) = matrix * get_vertex(2);
		r.get_normal()  = matrix * m_normal;
		r.get_radius()  = m_radius;
		return r;
	}


	t_raw_triangle::m_coor_type t_raw_triangle::calc_radius() const {
		// use index-triangle ctor to calculate radius
		const t_tuple4ui indices = t_tuple4ui(0, 1, 2);
		const t_idx_triangle triangle = t_idx_triangle(&m_verts[0], indices, m_normal);
		return (triangle.get_radius());
	}

	t_raw_triangle::m_coor_type t_raw_triangle::origin_distance() const {
		const m_coor_type p_a = ((m_verts[0]).to_vector()).inner(m_normal);
		const m_coor_type p_b = ((m_verts[1]).to_vector()).inner(m_normal);
		const m_coor_type p_c = ((m_verts[2]).to_vector()).inner(m_normal);
		// can return any of p_*, averaging just reduces error
		return ((p_a + p_b + p_c) / m_coor_type(3));
	}


	bool t_raw_triangle::point_in_triangle(const m_point_type& point, const m_coor_type eps) const {
		// construct planes (parallel to m_normal) for each triangle
		// edge (note that we can not use <point> itself to do this)
		const m_vector_type ab_plane_n = (m_verts[1] - m_verts[0]).outer(m_normal);
		const m_vector_type bc_plane_n = (m_verts[2] - m_verts[1]).outer(m_normal);
		const m_vector_type ca_plane_n = (m_verts[0] - m_verts[2]).outer(m_normal);

		const m_vector_type vp =      point.to_vector();
		const m_vector_type vb = m_verts[1].to_vector();
		const m_vector_type vc = m_verts[2].to_vector();
		const m_vector_type va = m_verts[0].to_vector();

		// test if <point> is above each plane, use GEQ so edges count
		// NOTE:
		//   the edge-plane normals (*_n) are not normalized, therefore
		//   the projected point-plane distances are not scaled correctly
		//   (however signs are preserved)
		//   CALLER must ensure this is only called when a point is within
		//   epsilon-distance from triangle-plane, otherwise acts for prism
		const bool above_ab_plane = ((vp.inner(ab_plane_n) - vb.inner(ab_plane_n)) >= -eps);
		const bool above_bc_plane = ((vp.inner(bc_plane_n) - vc.inner(bc_plane_n)) >= -eps);
		const bool above_ca_plane = ((vp.inner(ca_plane_n) - va.inner(ca_plane_n)) >= -eps);

		return (above_ab_plane == above_bc_plane && above_bc_plane == above_ca_plane);
	}

	bool t_raw_triangle::point_in_triangle_bc(const m_point_type& point, const m_coor_type eps) const {
		// barycentric containment test
		//
		// any point <p> inside a triangle can be expressed as a linear
		// combination of its vertices <a>, <b> and <c> by three weights
		// alpha, beta, and gamma (p = alpha * a + beta * b + gamma * c)
		// which must all fall in the range [0,1] and sum to 1
		//
		// to find alpha and gamma analytically we can solve either the
		// linear system a + vba*alpha = p - vca*gamma (three equations
		// with two unknowns) or a + vca*gamma = p - vba*alpha (equally
		// constrained); this amounts to constructing the parallelogram
		// defined by <a> and <p> using vba and vca as basis-vectors

		#if 0
		const m_vector_type vba = m_verts[1] - m_verts[0]; // b - a
		const m_vector_type vca = m_verts[2] - m_verts[0]; // c - a

		// note: all ortho-vectors must point inward
		const m_vector_type vbao = vba.outer(m_normal) * m_coor_type( 1);
		const m_vector_type vcao = vca.outer(m_normal) * m_coor_type(-1);

		// .x=min(alpha), .y=max(alpha), .z=min(gamma), .w=max(gamma)
		m_vector_type ranges = {m_coor_type(-1), m_coor_type(1), m_coor_type(-1), m_coor_type(1)};
		m_vector_type params;

		// find alpha
		while ((ranges.y() - ranges.x()) > eps) {
			params.x() = (ranges.x() + ranges.y()) * 0.5f;

			const m_point_type tp = point - (vba * params.x());
			const m_vector_type tv = tp - m_verts[0];

			if (tv.inner(vcao) >= m_coor_type(0)) {
				ranges.x() = params.x();
			} else {
				ranges.y() = params.x();
			}
		}

		// find gamma
		while ((ranges.w() - ranges.z()) > eps) {
			params.y() = (ranges.z() + ranges.w()) * 0.5f;

			const m_point_type tp = point - (vca * params.y());
			const m_vector_type tv = tp - m_verts[0];

			if (tv.inner(vbao) >= m_coor_type(0)) {
				ranges.z() = params.y();
			} else {
				ranges.w() = params.y();
			}
		}

		// note: should be LEQ 1 (edge-inclusive), but that allows external points colinear with edges
		return (params.x() >= m_coor_type(0) && params.y() >= m_coor_type(0) && (params.x() + params.y()) < m_coor_type(1));

		#else

		const m_point_type& p = point;
		const m_point_type  q = p - m_normal * (p - m_verts[0]).inner(m_normal);

		for (unsigned int n = 0; n < 3; n++) {
			const m_vector_type vp  =                    q - m_verts[n]; // p - a
			const m_vector_type vba = m_verts[(n + 1) % 3] - m_verts[n]; // b - a
			const m_vector_type vca = m_verts[(n + 2) % 3] - m_verts[n]; // c - a

			// avoid the single possible degenerate case
			if (lib_math::absm(vba.x()) < eps)
				continue;

			// TODO: t_ray::calc_closest_point(ray)
			//   a.x + vca.x*gamma == p.x - vba.x*alpha -> p.x - a.x == vca.x*gamma + vba.x*alpha
			//   a.y + vca.y*gamma == p.y - vba.y*alpha -> p.y - a.y == vca.y*gamma + vba.y*alpha
			//   a.z + vca.z*gamma == p.z - vba.z*alpha -> p.z - a.z == vca.z*gamma + vba.z*alpha
			//
			//   vpx == vca.x*gamma + vba.x*alpha -> vpx - vca.x*gamma == vba.x*alpha -> (vpx - vca.x*gamma)/vba.x == alpha
			//   vpy == vca.y*gamma + vba.y*alpha -> vpy - vca.y*gamma == vba.y*alpha
			//
			//   (vpy - vca.y*gamma)       == vba.y * ((vpx - vca.x*gamma)/vba.x)
			//   (vpy - vca.y*gamma)*vba.x == vba.y *  (vpx - vca.x*gamma)
			//
			//   vpy*vba.x - vca.y*vba.x*gamma == vpx*vba.y - vca.x*vba.y*gamma
			//    ky       - qyx        *gamma ==  kx       - qxy        *gamma
			//
			//    ky      == kx - qxy*gamma + qyx *gamma
			//    ky - kx ==    - qxy*gamma + qyx *gamma
			//    ky - kx ==    (-qxy       + qyx)*gamma
			//
			//                   (ky - kx) / (-qxy + qyx)           == gamma
			//   (vpx - vca.x * ((ky - kx) / (-qxy + qyx))) / vba.x == alpha
			//
			const m_coor_type kxy = vp.x() * vba.y();
			const m_coor_type kyx = vp.y() * vba.x();
			const m_coor_type qyx = vca.y() * vba.x();
			const m_coor_type qxy = vca.x() * vba.y();

			const m_coor_type gamma = (kyx - kxy) / (-qxy + qyx);
			const m_coor_type alpha = (vp.x() - vca.x() * gamma) / vba.x();

			return (alpha >= m_coor_type(0) && gamma >= m_coor_type(0) && (alpha + gamma) <= m_coor_type(1));
		}
		#endif

		return false;
	}


	// test if <this> intersects <plane>
	// TODO: return intersection-points
	bool t_raw_triangle::intersect_plane(const t_plane& plane, const m_coor_type eps) const {
		if ((m_coor_type(1) - std::fabs(m_normal.inner(plane.get_normal()))) <= eps)
			return false;

		// we intersect <plane> if and only if not all three vertices
		// are on the same side (above or below) of it, ie. if the sum
		// of their distance signs is between 0 and 3
		const int sign_va = (plane.point_distance(m_verts[0]) >= -eps);
		const int sign_vb = (plane.point_distance(m_verts[1]) >= -eps);
		const int sign_vc = (plane.point_distance(m_verts[2]) >= -eps);

		return ((sign_va + sign_vb + sign_vc) > 0 && (sign_va + sign_vb + sign_vc) < 3);
	}


	unsigned int t_raw_triangle::intersect_triangle_cop(const t_raw_triangle& triangle, m_point_type* points, const m_coor_type eps) const {
		unsigned int n = 0;

		#if 1
		if (std::fabs((triangle.get_vertex(0) - m_verts[0]).inner(m_normal)) > MAX_TRIANGLE_SLACK_SPACE)
			return n;
		#endif

		// <this> might be fully inside <triangle>; the partially
		// inside and fully outside cases are covered by int_rays
		for (unsigned int k = 0; k < 3; k++) {
			points[n]     = m_verts[k];
			points[n].w() = triangle.point_in_triangle(m_verts[k]);

			// if point is valid, advance index
			n += (points[n].w() == m_coor_type(1));
		}

		// if all three vertices of <this> are inside <triangle>,
		// we do not need to perform any edge-intersection tests
		if (n == 3)
			return n;

		// if only one or two points of <this> are inside <triangle>
		// (but not zero), one or two edges of <this> must intersect
		// edges of <triangle>
		//
		// intersect each edge of <this> with each edge of <triangle>
		for (unsigned int i = 0; (i < 3 && n < get_max_intersect_points()); i++) {
			const t_ray r0 = t_ray(get_vertex(i), get_vertex((i + 1) % 3));

			for (unsigned int j = 0; (j < 3 && n < get_max_intersect_points()); j++) {
				const t_ray r1 = t_ray(triangle.get_vertex(j), triangle.get_vertex((j + 1) % 3));

				// parallel edges can not intersect
				if (std::fabs((r0.dir()).inner(r1.dir())) > (m_coor_type(1) - eps))
					continue;

				if (!r0.ray_intersect(r1, &points[n], eps))
					continue;

				// the intersection-point must lie inside *both* triangles
				// (in case the triangles do not overlap, a point might be
				// located on an edge for triangle A but not for B)
				n += (point_in_triangle(points[n]) && triangle.point_in_triangle(points[n]));
			}
		}

		return n;
	}

	// test if <this> intersects <triangle> (note that
	// triangles must be in the same coordinate-space)
	//
	// if so, exactly one or two edges of <this> will
	// cross the plane in which <triangle> is embedded
	// and their intersection point will lie inside it
	//
	// note that is NOT a symmetric test: triangle A
	// might intersect triangle B but not vice versa
	//
	unsigned int t_raw_triangle::intersect_triangle(const t_raw_triangle& triangle, m_point_type* points, const m_coor_type eps) const {
		unsigned int n = 0;

		// early-out #1; bail if bounding-spheres do not overlap
		if ((calc_midpoint() - triangle.calc_midpoint()).sq_len() > lib_math::square(m_radius + triangle.get_radius()))
			return n;
		// early-out #2; bail if <this> does not intersect <triangle>'s plane
		if (!intersect_plane(triangle.to_plane()))
			return n;

		// first test the parallel case; if triangles overlap and
		// are not vertically separated, we have four contacts at
		// most
		if ((m_coor_type(1) - std::fabs(m_normal.inner(triangle.get_normal()))) <= eps)
			return (n = intersect_triangle_cop(triangle, points, eps));

		const t_ray ray_ab = t_ray(m_verts[0], m_verts[1]), ray_ba = ray_ab.invert_pos();
		const t_ray ray_bc = t_ray(m_verts[1], m_verts[2]), ray_cb = ray_bc.invert_pos();
		const t_ray ray_ca = t_ray(m_verts[2], m_verts[0]), ray_ac = ray_ca.invert_pos();

		// triangles are one-sided so an edge can only intersect
		// if it runs ANTI-parallel to a triangle's normal vector
		// (if ray_xy intersects then ray_yx is guaranteed not to)
		//
		// ray_xy runs FROM vertex x TO vertex y (pos=x, dir=y-x)
		if ((ray_ab.dir()).inner(triangle.get_normal()) >= m_coor_type(0)) {
			points[n] = ray_ba.triangle_intersection_point(triangle);
		} else {
			points[n] = ray_ab.triangle_intersection_point(triangle);
		}

		// if point is valid, advance index
		n += (points[n].w() == m_coor_type(1));

		if ((ray_bc.dir()).inner(triangle.get_normal()) >= m_coor_type(0)) {
			points[n] = ray_cb.triangle_intersection_point(triangle);
		} else {
			points[n] = ray_bc.triangle_intersection_point(triangle);
		}

		n += (points[n].w() == m_coor_type(1));

		if ((ray_ca.dir()).inner(triangle.get_normal()) >= m_coor_type(0)) {
			points[n] = ray_ac.triangle_intersection_point(triangle);
		} else {
			points[n] = ray_ca.triangle_intersection_point(triangle);
		}

		n += (points[n].w() == m_coor_type(1));

		// return the number of valid intersection points
		return n;
	}

	unsigned int t_raw_triangle::intersect_triangle(
		const t_idx_triangle& idx_triangle,
		const m_point_type* idx_verts,
		      m_point_type* points,
		const m_coor_type eps
	) const {
		return (intersect_triangle(idx_triangle.to_raw_triangle(idx_verts), points, eps));
	}
};

