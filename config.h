/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Parallel execution */
#define EXEC_PARALLEL 1

/* Serial execution */
#define EXEC_SERIAL 0

/* Choice of simulation execution type */
#define EXEC_TYPE EXEC_PARALLEL

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `memset' function. */
#define HAVE_MEMSET 1

/* Define to 1 if you have the `pow' function. */
#define HAVE_POW 1

/* Define to 1 if you have the `sqrt' function. */
#define HAVE_SQRT 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if the system has the type `_Bool'. */
/* #undef HAVE__BOOL */

/* Name of package */
#define PACKAGE "spectrum_batl"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT ""

/* Define to the full name of this package. */
#define PACKAGE_NAME "spectrum_batl"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "spectrum_batl 1.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "spectrum_batl"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.0"

/* Choice of the Runge-Kutta method */
#define RK_INTEGRATOR_TYPE 0

/* Cartesian server with AMR (part of BATS-R-US) */
#define SERVER_BATL 301

/* Cartesian server with uniform grid */
#define SERVER_CARTESIAN 300

/* Choice of background server interpolation order */
#define SERVER_INTERP_ORDER 1

/* Choice of background server number of ghost cells per side of block */
#define SERVER_NUM_GHOST_CELLS 1

/* No server */
#define SERVER_SELF 299

/* Choice of the background server */
#define SERVER_TYPE SERVER_BATL

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Field line tracer (-,-,-) */
#define TRAJ_FIELDLINE 199

/* Focused transport model (p,mu,-) */
#define TRAJ_FOCUSED 205

/* Guiding center model (p_para,-,p_perp) */
#define TRAJ_GUIDING 201

/* Guiding center model with perp. diffusion (p_para,-,p_perp) */
#define TRAJ_GUIDING_DIFF 203

/* Guiding center model with perp. diffusion and PA scattering
   (p_para,-,p_perp) */
#define TRAJ_GUIDING_DIFF_SCATT 204

/* Guiding center model with PA scattering (p_para,-,p_perp) */
#define TRAJ_GUIDING_SCATT 202

/* Newton-Lorentz model (px,py,pz) */
#define TRAJ_LORENTZ 200

/* Isotropic model (p,-,-) */
#define TRAJ_PARKER 206

/* Isotropic model with source (p,-,-) */
#define TRAJ_PARKER_SOURCE 207

/* Choice of trajectory time flow direction */
#define TRAJ_TIME_FLOW TRAJ_TIME_FLOW_BACKWARD

/* Backward-in-time simulation */
#define TRAJ_TIME_FLOW_BACKWARD 1

/* Forward-in-time simulation */
#define TRAJ_TIME_FLOW_FORWARD 0

/* Choice of the trajectory physics */
#define TRAJ_TYPE TRAJ_PARKER

/* Using CUDA */
/* #undef USE_CUDA */

/* Using GSL */
#define USE_GSL 1

/* Using MPI */
#define USE_MPI 1

/* Using SILO */
/* #undef USE_SILO */

/* Using SLURM */
#define USE_SLURM 1

/* Version number of package */
#define VERSION "1.0"

/* Define for Solaris 2.5.1 so the uint32_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef were allowed, the
   #define below would cause a syntax error. */
/* #undef _UINT32_T */

/* Define for Solaris 2.5.1 so the uint8_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef were allowed, the
   #define below would cause a syntax error. */
/* #undef _UINT8_T */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to the type of an unsigned integer type of width exactly 32 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint32_t */

/* Define to the type of an unsigned integer type of width exactly 8 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint8_t */
