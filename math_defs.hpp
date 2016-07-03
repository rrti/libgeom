#ifndef LIBGEOM_DEFS_HDR
#define LIBGEOM_DEFS_HDR

#define LIBGEOM_POINT_SIZE  (4)
#define LIBGEOM_VECTOR_SIZE (4)
#define LIBGEOM_TUPLE_SIZE  (4)
#define LIBGEOM_MATRIX_SIZE (4 * 4)

#define M_FEPS 1e-4f
#define M_DEPS 1e-4
#define M_FINF std::numeric_limits< float>::infinity()
#define M_DINF std::numeric_limits<double>::infinity()
#define M_FNAN std::numeric_limits< float>::signaling_NaN()
#define M_DNAN std::numeric_limits<double>::signaling_NaN()

// slack-space allowed between (coplanar) triangle surfaces
#define MAX_TRIANGLE_SLACK_SPACE 0.05f

#endif

