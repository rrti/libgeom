#ifndef LIBGEOM_DEFS_HDR
#define LIBGEOM_DEFS_HDR

#define M_FEPS 1e-4f
#define M_DEPS 1e-4
#define M_FINF std::numeric_limits< float>::infinity()
#define M_DINF std::numeric_limits<double>::infinity()
#define M_FNAN std::numeric_limits< float>::signaling_NaN()
#define M_DNAN std::numeric_limits<double>::signaling_NaN()

// slack-space allowed between (coplanar) triangle surfaces
#define MAX_TRIANGLE_SLACK_SPACE 0.05f

#define X_AXIS_IDX 0
#define Z_AXIS_IDX 2

#define P_AXIS_IDX 0 // pitch (x-axis, [0])
#define Y_AXIS_IDX 1 //   yaw (y-axis, [1])
#define R_AXIS_IDX 2 //  roll (z-axis, [2])

#endif

