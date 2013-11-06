#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846
#define NW 5

#define USAGE "lightcurve_v8 -robj 0.25 -d 40  -star 0.02 -f r -b 0 -out outfile "

double km_per_AU = 1.5e8;
double km_per_nm = 1.0e-12;
double radians_per_arcsec = 1.0/206265.;

int main(int argc, char **argv)
{
    int i, j, k;
    double x, y, xmax, omega;
    double rho, rho_star;
    double eta, eta_p;
    double eta_0, eta_1;
    double sep, sep_max, dsep;
    double g, *g_sum;
    double intensity(double rho, double eta);
    double I, *I_array;
    double I_0, I_1, I_eta;
    double *lightcurve;
    int idx, idx_min, idx_max;
    double area;
    double w;            // characteristic event width
    double tmax = 5;     // +/- tmax is lightcurve duration (sec)
    double sp = 40.0;    // sampling frequency (Hz)
    double as = 0.0;     // angular size of star (mas)
    double d = 40.0;     // distance to occulting object (AU)
    double b = 0.0;      // impact parameter in fraction of stellar radius
    double robj = 1.0;   // radius of occulting object (km)
    double rstar = 0.0;  // projected radius of star in occulter plane.
    double mw = 600.0;   // mean wavelength (nm)
    double vt = 25.0;    // transverse velocity (km/sec)
    double offset = 0.0; // time offset (sec)
    double F = 0.0;      // Fresnel scale
    double w1 = 552.0;   // short wavelength edge
    double w2 = 689.0;   // long wavelength edge
    double m = 2.0;      // number of event widths W to compute
    int    n_per_F = 100;// number of outputs per Fresnel scale
    int    n;            // number of outputs, calculated from other parameters
    double W;            // width of event, from Nihei et al 2007.
    double dt;

    FILE *outfile;

    /* Handle options */
    if(argc<4) { 
	printf("Usage: %s\n", USAGE);
	exit(0);
    }
    for(k=1;k<argc; k++) {
	if(strcmp(argv[k], "-h")==0) {
	    printf("Usage: %s \n", USAGE); 
	    exit(0);
	}
	if(strcmp(argv[k], "-robj")==0)  {robj   =atof(argv[k+1]);} 
	if(strcmp(argv[k], "-d")==0)     {d      =atof(argv[k+1]);}
	if(strcmp(argv[k], "-star")==0)  {as     =atof(argv[k+1]);}
	if(strcmp(argv[k], "-sp")==0)    {sp     =atof(argv[k+1]);}
	if(strcmp(argv[k], "-vt")==0)    {vt     =atof(argv[k+1]);}
	if(strcmp(argv[k], "-b")==0)     {b      =atof(argv[k+1]);}
	if(strcmp(argv[k], "-m")==0)     {m      =atof(argv[k+1]);}
	if(strcmp(argv[k], "-off")==0)   {offset =atof(argv[k+1]);}
	if(strcmp(argv[k], "-nf")==0)    {n_per_F=atof(argv[k+1]);}
	if(strcmp(argv[k], "-mw")==0)    {mw     =atof(argv[k+1]); w1 = mw; w2 = mw;}
	if(strcmp(argv[k], "-out")==0)   {outfile=fopen(argv[k+1], "w");}
	if(strcmp(argv[k], "-f")==0) { /* FILTERS */
	    if(strcmp(argv[k+1], "g")==0) {w1=405; w2=550;}
	    if(strcmp(argv[k+1], "r")==0) {w1=552; w2=689;}
	    if(strcmp(argv[k+1], "i")==0) {w1=691; w2=815;}
	    if(strcmp(argv[k+1], "z")==0) {w1=815; w2=915;} 
	    if(strcmp(argv[k+1], "y")==0) {w1=900; w2=1215;} 
	}
    }

    /* Internal units are km and sec. */

    mw = 0.5*(w1 + w2);                     // mean wavelength, for estimating Fresnel scale

    F = sqrt(0.5*mw*km_per_nm*d*km_per_AU); // Fresnel scale (km)

    rstar=0.5*as/1000*radians_per_arcsec*d*km_per_AU;   // projected stellar radius
    rho_star = rstar/F;                      // projected stellar radius in Fresnel units
    rho = robj/F;                            // object radius in Fresnel units

    /* omega is the characteristic width of the occultation event in Fresnel units, Nihei et al. 2007 */
    omega = 2.0*pow(pow(sqrt(3.0), 1.5) + pow(rho, 1.5), 2.0/3.0) + 2.0*rho_star;
    W = omega*F;                             // characteristic width of occultation event, physical units

    n = m*omega*n_per_F;                     // number of points in half-event 
    tmax = m*W/vt;                           // time span of half-event
    xmax = vt*tmax;                          // distance moved in time span of half-event

    // maximum separation between center of occulter and a point on the face of the star.
    sep_max = sqrt(xmax*xmax + b*W*b*W) + rstar;  
    dsep = sep_max/(n-1);

    dt = dsep/vt;

    /* Allocate arrays    */
    g_sum      = calloc(2*n + 3, sizeof(double));
    I_array    = calloc(n+1, sizeof(double));
    lightcurve = calloc(2*n + 3, sizeof(double));

    /* Initialize arrays */
    for(i=0; i<=n; i++){
	I_array[i] = 0.0;
	g_sum[i] = 0.0;
	g_sum[2*n-i] = 0.0;
    }

    if(2.0*rstar < dsep){
	area = 1.0;                          // considering start to be a point source.
    }else{
	area = PI*rstar*rstar;               // projected area of star (km^2)
    }

    fprintf(outfile, "# %lf %lf %lf %lf\n", 1.0/dt, W, F, rstar);

    /* Calculate lightcurve values for an array of uniformly spaced values. */
    /* These are averaged over a number of wavelength bins, assuming a neutral */
    /* color.  The values get used in the final lightcurve calculation below. */
    for(i=0; i<n; i++){
	sep = i*dsep;
	for(j=0; j<NW; j++){
	    w = w1 + j*((w2-w1)/NW);                // wavelength
	    F = sqrt(0.5*w*km_per_nm*d*km_per_AU); // Fresnel scale (km)
	    I_array[i] += intensity(robj/F, sep/F);
	}
	I_array[i] /= NW;
    }

    /* Calculate the final lightcurve, integrating over the face of the star, for the full time span. */
    for(i=0; i<=n; i++){
	x = -xmax + xmax*(float) i/n;
	y = b*W;
	eta = sqrt(x*x + y*y);
	I = 0.0;
	g_sum[i] = 0.0;

	if(2.0*rstar < dsep){ /* Very small star; linearly interpolate to a single best value */

	  idx = (int) floor(eta/dsep);
	  eta_0 = dsep*idx;
	  I_0 = I_array[idx];
	  eta_1 = eta_0 + dsep;
	  I_1 = I_array[idx+1];
	  I_eta = I_0 + ((I_1-I_0)/dsep)*(eta-eta_0);
	  g_sum[i] += 1.0;
	  /*I += intensity(robj/F, eta/F);*/
	  I += I_eta;

	}else{ /* Star is resolved; integrate over face of the star */

	  if(eta > rstar){ /* Center of object is outside of limb of star. */

	    idx_min = (int) ((eta - rstar)/dsep) + 1;
	    idx_max = (int) ((eta  + rstar)/dsep);

	  }else{ /* Center of object is inside of limb of star. */

	    idx_min = 1;
	    idx_max = (int) ((eta  + rstar)/dsep);

	  }

	  for(j=idx_min; j<idx_max; j++){
	    eta_p = dsep*j;
	    if(eta_p < (rstar - eta)){
	      g = 2.0*PI;
	    }else{
	      g = 2*acos((eta*eta + eta_p*eta_p - rstar*rstar)/(2.0*eta*eta_p));
	    }
	    g_sum[i] += g*eta_p*dsep;
	    I += g*eta_p*I_array[j]*dsep;
	  }
	}

	g_sum[2*n-i] = g_sum[i];
	lightcurve[i] = I/g_sum[i];
	lightcurve[2*n-i] = I/g_sum[2*n-i];

    }

    for(i=0; i<=2*n; i++){
      	x = -xmax + xmax*(float) i/n;
	fprintf(outfile, "%.8lf %.8lf %.8lf %.8lf\n", x/vt, g_sum[i]/area, lightcurve[i], x);
    }

    return 0;

}

double intensity(double rho, double eta)
{

  double u0, u1, u2;
  double x, val;
  double lommel(int n, double mu, double nu);

  if(eta <= rho){

    u0 = lommel(0, eta, rho);
    u1 = lommel(1, eta, rho);

    val = u0*u0 + u1*u1;

  }else{

    u1 = lommel(1, rho, eta);
    u2 = lommel(2, rho, eta);

    x = 0.5*PI*(rho*rho + eta*eta);

    val = 1.0 + u1*u1 + u2*u2 - 2.0*u1*sin(x) + 2.0*u2*cos(x);

  }

  return val;

}

double lommel(int n, double mu, double nu)
{
  
  double val, x;
  int k, kmax;
  double eps = 1e-15;

  kmax = 100;

   val = 0.0;

  for(k=0; k<=kmax; k++){

    x = pow(-1.0, k)*pow(mu/nu, n + 2*k) * jn(n+2*k, PI*mu*nu);

    if(fabs(x) < eps){
      val += x;
      break;
    }else{
      val += x;
    }
  }

  return val;

}

  

