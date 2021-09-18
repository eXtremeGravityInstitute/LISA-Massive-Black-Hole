/* ************************************************************************ */
/*                                                                          */
/* TITLE: Preprocessor Macros For The Galactic Background                   */
/*                                                                          */
/* ABSTRACT: This header file contains a number of constants used           */
/* throughout the making of the galactic backgound realizations.            */
/*                                                                          */
/* ************************************************************************ */



          /* --------------  MATHEMATICAL CONSTANTS  -------------- */

 /* Set the value of pi */
#define pi 3.141592653589793

#define sq2pi 2.506628274631

#define TPI 6.2831853071795865

 /* Square root of 3 */
#define sq3 1.73205080757
          /* ----------------  NATURAL CONSTANTS  ----------------- */

 /* Speed of light (m/s) */
#define clight 299792458.

 /* Number of seconds in a sidereal year */
#define year 3.15581498e7

/* Number of seconds in a day */
#define day 86400.0

/* Number of seconds in a week */
#define week 604800.0

/* Number of second in hour */
#define hour 3600.0

 /* Newton's gravitational constant (mks) */
#define G 6.67259e-11

 /* Astronomical unit (meters) */
#define au 1.49597870660e11

 /* Number of meters in a parsec */
#define pc 3.0856775807e16

 /* Number of meters in a kiloparsec */
#define kpc 3.0856775807e19

 /* Distance bewteen the Sun and the center of the galaxy (kpc) */
#define Rsun 8.5

 /* Mass of the Sun (kg) */
#define Msun 1.9889e30

/* Mass of the Sun (s) */
#define Tsun 4.92569043916e-6

#define h22fac  0.31539156525252   //  2.*sqrt(5./(64.*PI)) factor for h22 to h conversion

#define STtoSQD 3282.8063500117437948  // steradians to square degree

          /* ----------------  DETECTOR CONSTANTS  ---------------- */

/* how far back in time to start */
#define Tpad 1.0e6

/* sample cadence in seconds */
#define cadence 5.0

/* segemnet duration in seconds */
#define Tsegment 2621440.0


 /* Initial azimuthal position of the guiding center */
#define kappa0 0.0

 /* Initial orientation of the LISA constellation */
#define lambda0 0.0

 /* Orbital radius of the guiding center */
#define Rgc (1.0*au)

 /* Mean arm length of the LISA detector (meters) */
#define Larm 2.5e9

 /* Transfer frequency */ 
#define fstar 0.01908538064

 /* LISA orbital eccentricity */
#define ec 0.0048241852175

 /* LISA modulation frequency */
#define fm 3.168753575e-8

#define sq3 1.73205080757

#define GEOM 4.92569043916e-6

#define NParams 11

#define NI 6

#define NE 7

#define MF 4096

/* set to 1 to run on LDC, 0 for a cataliog source */
#define LDC 1

/* set to 1 for noise level marginalization, 0 otherwise */
#define nflag 0

/* set to zero for no cut to data */
#define cut 0

/* length of cosine filter used to cut data */
#define t_tuke 5.0e4

/* number of cold chains */
#define NCC 4

/* set to 1.0e-3 for v1, 0.0 for v2 */
#define PDfref 0.0

/* set to 1 for prior recovery test */
#define lhold 0

/* set to 0 for Massive Black holes, 1 for SOBBHs */
#define sflag 0

#define etamin 0.05     // lower limit on the symmetric mass ratio (18:1 mass ratio)


