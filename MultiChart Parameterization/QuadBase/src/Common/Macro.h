/* ================== Library Information ================== */
// [Name] 
// MeshLib Library
//
// [Developer]
// Xu Dong
// State Key Lab of CAD&CG, Zhejiang University
// 
// [Date]
// 2005-08-05
//
// [Goal]
// A general, flexible, versatile and easy-to-use mesh library for research purpose.
// Supporting arbitrary polygonal meshes as input, but with a 
// primary focus on triangle meshes.

/* ================== File Information ================== */
// [Name]
// Macro.h
//
// [Developer]
// Xu Dong
// State Key Lab of CAD&CG, Zhejiang University
// 
// [Date]
// 2005-08-05
//
// [Goal]
// Defining the global macros and enumerations.



#pragma once


#define PLATFORM_WINDOWS 1

#ifdef PLATFORM_WINDOWS
#define MESHLIB_ASSERT(a) ASSERT(a)
#endif

/* ================== Enumerations ================== */

typedef enum {SOLID_SMOOTH, SOLID_FLAT, WIREFRAME, SOLID_AND_WIREFRAME, VERTICES, HIGHLIGHT_ONLY, TEXTURE_MAPPING,
	SELECT_NULL, SELECT_VERTEX, SELECT_EDGE, SELECT_FACE, } RENDER_MODE;
 
typedef enum {TRANSLATE, SCALE, ROTATE, NONE} TRANSFORM_MODE;

typedef enum {FACE, EDGE, VERTEX} PICKING_MODE;

typedef enum {AXIS_X, AXIS_Y, AXIS_Z} AXIS;

typedef enum {R=0, G, B} COLOR;



/* ================== Constants ================== */

#define PI 3.1415926535897932384626433832795
#define PI_12 1.5707963267948966192313216916395
#define EXP 2.7
#define LARGE_ZERO_EPSILON 1.0E-12
#define SMALL_ZERO_EPSILON 1.0E-16



/* ================== Math Macros ================== */

// Min & Max & ABS
#define MAX3(a, b, c) (((a) > (b)) ? (((a) > (c)) ? (a) : (c)) : (((b) > (c)) ? (b) : (c)))
#define MAX(a, b)	  (((a) > (b)) ? (a) : (b))

#define ABS(a) (((a) > -(a)) ? (a) : (-(a)))

#define MIN3(a, b, c) (((a) < (b)) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))
#define MIN(a, b)	  (((a) < (b)) ? (a) : (b))

// Equal within epsilon tolerance
#define ALMOST_EQUAL_LARGE(x,y) (fabs((x)-(y)) < LARGE_ZERO_EPSILON)
#define ALMOST_EQUAL_SMALL(x,y) (fabs((x)-(y)) < SMALL_ZERO_EPSILON)

// Degree & Radius
#define DEG2RAD(d) ((d)*0.017453292519943295769236907684886)
#define RAD2DEG(r) ((r)*57.29577951308232087679815481410500)



