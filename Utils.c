#include "Constants.h"
#include "IMRPhenomD.h"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Declarations.h"

int *int_vector(int N)
{
    return malloc( (N+1) * sizeof(int) );
}

void free_int_vector(int *v)
{
    free(v);
}

double *double_vector(int N)
{
    return malloc( (N+1) * sizeof(double) );
}

void free_double_vector(double *v)
{
    free(v);
}

double **double_matrix(int N, int M)
{
    int i;
    double **m = malloc( (N+1) * sizeof(double *));
    
    for(i=0; i<N+1; i++)
    {
        m[i] = malloc( (M+1) * sizeof(double));
    }
    
    return m;
}

void free_double_matrix(double **m, int N)
{
    int i;
    for(i=0; i<N+1; i++) free_double_vector(m[i]);
    free(m);
}


int **int_matrix(int N, int M)
{
    int i;
    int **m = malloc( (N+1) * sizeof(int *));
    
    for(i=0; i<N+1; i++)
    {
        m[i] = malloc( (M+1) * sizeof(int));
    }
    
    return m;
}

void free_int_matrix(int **m, int N)
{
    int i;
    for(i=0; i<N+1; i++) free_int_vector(m[i]);
    free(m);
}

double ***double_tensor(int N, int M, int L)
{
    int i,j;
    
    double ***t = malloc( (N+1) * sizeof(double **));
    for(i=0; i<N+1; i++)
    {
        t[i] = malloc( (M+1) * sizeof(double *));
        for(j=0; j<M+1; j++)
        {
            t[i][j] = malloc( (L+1) * sizeof(double));
        }
    }
    
    return t;
}

void free_double_tensor(double ***t, int N, int M)
{
    int i;
    
    for(i=0; i<N+1; i++) free_double_matrix(t[i],M);
    
    free(t);
}

void ang2pix_ring( const long nside, double theta, double phi, long *ipix) {
    /*
     c=======================================================================
     c     gives the pixel number ipix (RING)
     c     corresponding to angles theta and phi
     c=======================================================================
     */
    
    int nl2, nl4, ncap, npix, jp, jm, ipix1;
    double  z, za, tt, tp, tmp;
    int ir, ip, kshift;
    
    double piover2 = 0.5*M_PI;
    double twopi=2.0*M_PI;
    double z0=2.0/3.0;
    long ns_max=8192;
    
    if( nside<1 || nside>ns_max ) {
        fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
        exit(0);
    }
    
    if( theta<0. || theta>PI) {
        fprintf(stderr, "%s (%d): theta out of range: %f\n", __FILE__, __LINE__, theta);
        exit(0);
    }
    
    z = cos(theta);
    za = fabs(z);
    if( phi >= twopi)  phi = phi - twopi;
    if (phi < 0.)     phi = phi + twopi;
    tt = phi / piover2;//  ! in [0,4)
    
    nl2 = 2*nside;
    nl4 = 4*nside;
    ncap  = nl2*(nside-1);// ! number of pixels in the north polar cap
    npix  = 12*nside*nside;
    
    if( za <= z0 ) {
        
        jp = (int)floor(nside*(0.5 + tt - z*0.75)); /*index of ascending edge line*/
        jm = (int)floor(nside*(0.5 + tt + z*0.75)); /*index of descending edge line*/
        
        ir = nside + 1 + jp - jm;// ! in {1,2n+1} (ring number counted from z=2/3)
        kshift = 0;
        if (fmod(ir,2)==0.) kshift = 1;// ! kshift=1 if ir even, 0 otherwise
        
        ip = (int)floor( ( jp+jm - nside + kshift + 1 ) / 2 ) + 1;// ! in {1,4n}
        if( ip>nl4 ) ip = ip - nl4;
        
        ipix1 = ncap + nl4*(ir-1) + ip ;
    }
    else {
        
        tp = tt - floor(tt);//      !MOD(tt,1.d0)
        tmp = sqrt( 3.*(1. - za) );
        
        jp = (int)floor( nside * tp * tmp );// ! increasing edge line index
        jm = (int)floor( nside * (1. - tp) * tmp );// ! decreasing edge line index
        
        ir = jp + jm + 1;//        ! ring number counted from the closest pole
        ip = (int)floor( tt * ir ) + 1;// ! in {1,4*ir}
        if( ip>4*ir ) ip = ip - 4*ir;
        
        ipix1 = 2*ir*(ir-1) + ip;
        if( z<=0. ) {
            ipix1 = npix - 2*ir*(ir+1) + ip;
        }
    }
    *ipix = ipix1 - 1;// ! in {0, npix-1}
    
}

void pix2ang_ring( long nside, long ipix, double *theta, double *phi) {
    /*
     c=======================================================================
     c     gives theta and phi corresponding to pixel ipix (RING)
     c     for a parameter nside
     c=======================================================================
     */
    
    int nl2, nl4, npix, ncap, iring, iphi, ip, ipix1;
    double  fact1, fact2, fodd, hip, fihip;
    //      PARAMETER (pi     = 3.1415926535897932384626434d0)
    //      parameter (ns_max = 8192) ! 2^13 : largest nside available
    
    int ns_max=8192;
    
    if( nside<1 || nside>ns_max ) {
        fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
        exit(0);
    }
    npix = 12*nside*nside;      // ! total number of points
    if( ipix<0 || ipix>npix-1 ) {
        fprintf(stderr, "%s (%d): ipix out of range: %ld\n", __FILE__, __LINE__, ipix);
        exit(0);
    }
    
    ipix1 = ipix + 1; // in {1, npix}
    nl2 = 2*nside;
    nl4 = 4*nside;
    ncap = 2*nside*(nside-1);// ! points in each polar cap, =0 for nside =1
    fact1 = 1.5*nside;
    fact2 = 3.0*nside*nside;
    
    if( ipix1 <= ncap ) {  //! North Polar cap -------------
        
        hip   = ipix1/2.;
        fihip = floor(hip);
        iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;// ! counted from North pole
        iphi  = ipix1 - 2*iring*(iring - 1);
        
        *theta = acos( 1. - iring*iring / fact2 );
        *phi   = (1.*iphi - 0.5) * PI/(2.*iring);
    }
    else if( ipix1 <= nl2*(5*nside+1) ) {//then ! Equatorial region ------
        
        ip    = ipix1 - ncap - 1;
        iring = (int)floor( ip / nl4 ) + nside;// ! counted from North pole
        iphi  = (int)fmod(ip,nl4) + 1;
        
        fodd  = 0.5 * (1 + fmod((double)(iring+nside),2));//  ! 1 if iring+nside is odd, 1/2 otherwise
        *theta = acos( (nl2 - iring) / fact1 );
        *phi   = (1.*iphi - fodd) * PI /(2.*nside);
    }
    else {//! South Polar cap -----------------------------------
        
        ip    = npix - ipix1 + 1;
        hip   = ip/2.;
        /* bug corrige floor instead of 1.* */
        fihip = floor(hip);
        iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;//     ! counted from South pole
        iphi  = (int)(4.*iring + 1 - (ip - 2.*iring*(iring-1)));
        
        *theta = acos( -1. + iring*iring / fact2 );
        *phi   = (1.*iphi - 0.5) * PI/(2.*iring);
    }
}

void legendre_maker(int J, int U, double **P)
{
    double *gs, *a, *b;
    double c0, c1, c2, d0, d1, d2;
    int i, k;
    
    gs = double_vector(J+1);  // downsampled function
    a = double_vector(J+1);  // coefficients
    b = double_vector(J+1);  // data

    for (k = 0; k <= U; ++k)
    {
       
     c0 = 0.0;
     c1 = (double)(U-2*k);
     c2 = (double)(U);
     d0 = c2;
     d1 = c1+c1;
     d2 = d0;

     P[0][k] = 1.0;
     P[1][k] = c1/(double)(U);
 
    for (i = 1; i < J; ++i)
     {
        d0 += 2.0;
        d2 -= 2.0;
        c0 += d0;
        c1 += d1;
        c2 += d2;
        P[i+1][k] = (c1*P[i][k]-c0*P[i-1][k])/c2;
    }
        
    }
        
    free_double_vector(gs);
    free_double_vector(a);
    free_double_vector(b);

}

double fourier_nwip2(double *a, double *b, double *Sn, int imin, int imax, int N)
{
    int i, j, k;
    double arg, product;
    double test;
    double ReA, ReB, ImA, ImB;
    
    arg = 0.0;
    for(i=imin; i<imax; i++)
    {
        j = i;
        k = N-i;
        ReA = a[j]; ImA = a[k];
        ReB = b[j]; ImB = b[k];
        product = ReA*ReB + ImA*ImB;
        arg += product/(Sn[i]);
    }
    
    return(4.0*arg);
    
}

double fourier_nwip(double *a, double *b, double *Sn, int n)
{
    int i, j, k;
    double arg, product;
    double test;
    double ReA, ReB, ImA, ImB;
    
    arg = 0.0;
    for(i=1; i<n/2; i++)
    {
        j = i;
        k = n-i;
        ReA = a[j]; ImA = a[k];
        ReB = b[j]; ImB = b[k];
        product = ReA*ReB + ImA*ImB;
        arg += product/(Sn[i]);
    }
    
    return(4.0*arg);
    
}

void fourier_nwip_dual_time(double *abt, double *aA, double *bA, double *aE, double *bE, double *Sn, double Tobs, int n)
{
    int i, j, k;
    double arg, product;
    double test;
    double ReA, ReB, ImA, ImB;
    double scale;
    
    scale = 2.0*(double)(n)/Tobs;
    
    for(i=1; i<n/2; i++)
    {
        j = i;
        k = n-i;
        ReA = aA[j]; ImA = aA[k];
        ReB = bA[j]; ImB = bA[k];
        abt[j] = scale*(ReA*ReB + ImA*ImB)/Sn[i];
        abt[k] = scale*(ImA*ReB - ReA*ImB)/Sn[i];
        ReA = aE[j]; ImA = aE[k];
        ReB = bE[j]; ImB = bE[k];
        abt[j] += scale*(ReA*ReB + ImA*ImB)/Sn[i];
        abt[k] += scale*(ImA*ReB - ReA*ImB)/Sn[i];
    }
    
    abt[0]=0.0;
    abt[n/2]=0.0;
    
    gsl_fft_halfcomplex_radix2_inverse(abt, 1, n);
    
}

void fourier_nwip_time(double *abt, double *a, double *b, double *Sn, double Tobs, int n)
{
    int i, j, k;
    double arg, product;
    double test;
    double ReA, ReB, ImA, ImB;
    double scale;
    
    scale = 2.0*(double)(n)/Tobs;
    
    for(i=1; i<n/2; i++)
    {
        j = i;
        k = n-i;
        ReA = a[j]; ImA = a[k];
        ReB = b[j]; ImB = b[k];
        abt[j] = scale*(ReA*ReB + ImA*ImB)/Sn[i];
        abt[k] = scale*(ImA*ReB - ReA*ImB)/Sn[i];
    }
    
    abt[0]=0.0;
    abt[n/2]=0.0;
    
    gsl_fft_halfcomplex_radix2_inverse(abt, 1, n);
    
}



void pbt_shift(double *corr, double *corrf, double *data1, double *data2, double *Sn, int n)
{
    int nb2, i, l, k, j;
    int imax, imin;
    
    nb2 = n/2;
    
    for (i=1; i < nb2; i++)
    {
        l=i;
        k=n-i;
        
        corr[l]    = (data1[l]*data2[l] + data1[k]*data2[k])/Sn[i];
        corr[k]    = (data1[k]*data2[l] - data1[l]*data2[k])/Sn[i];
        corrf[l] = corr[k];
        corrf[k] = -corr[l];
    }
    
    corr[0] = 0.0;
    corrf[0] = 0.0;
    corr[nb2] = 0.0;
    corrf[nb2] = 0.0;
    
    
    gsl_fft_halfcomplex_radix2_inverse(corr, 1, n);
    gsl_fft_halfcomplex_radix2_inverse(corrf, 1, n);
    
    
}

void InChl(int ll, double *params, double **Fisher, double **iChl)
{
    int i, j, k, kk;
    double **CovI;
    double *scale;
    
    scale = (double*)malloc(sizeof(double)* (NI));
    
    CovI = double_matrix(NI,NI);

    Inverse(Fisher, CovI, NI);

   // Fisher uses log derivatives on m1, m2,  Want to switch to linear derivatives since using
   // uniform prior in m1, m2

    for(j = 0; j < NI; j++) scale[j] = 1.0;
    
    if(ll==0)
    {
    scale[0] = params[0];
    scale[1] = params[1];
    }
    else
    {
    scale[0] = exp(params[0]);
    scale[1] = exp(params[1]);
    }

    for(j = 0; j < NI; j++)
      {
       for(k=0; k < NI; k++)
       {
           CovI[j][k] *= scale[j]*scale[k];
       }
    }

    cholesky(CovI, iChl, NI);
  
    free(scale);
    free_double_matrix(CovI,NI);

    
}

void Ext_In(int ll, double *params, double **Fisher, double **eChl, double **iChl)
{
    int i, j, k, kk;
    double **FishE;
    double **CovE, **CovI;
    double *scale;
    
    kk = NP-NE;
    
    scale = (double*)malloc(sizeof(double)* (NP));
    
    FishE = double_matrix(NE,NE);
    CovE = double_matrix(NE,NE);
    CovI = double_matrix(NI,NI);
    
       for(j = 0; j < NE; j++)
         {
          for(k=0; k < NE; k++)
          {
              FishE[j][k] = Fisher[j+kk][k+kk];
          }
         }

    Inverse(Fisher, CovI, NI);
    Inverse(FishE, CovE, NE);

   // Fisher uses log derivatives on m1, m2, DL. Want to switch to linear derivatives since using
   // uniform prior in m1, m2 and DL

    for(j = 0; j < NP; j++) scale[j] = 1.0;
    
    if(ll==0)
    {
    scale[0] = params[0];
    scale[1] = params[1];
    scale[6] = params[6];
    }
    else
    {
    scale[0] = exp(params[0]);
    scale[1] = exp(params[1]);
    scale[6] = exp(params[6]);
    }

    for(j = 0; j < NI; j++)
      {
       for(k=0; k < NI; k++)
       {
           CovI[j][k] *= scale[j]*scale[k];
       }
    }
    
    for(j = 0; j < NE; j++)
      {
       for(k=0; k < NE; k++)
       {
           CovE[j][k] *= scale[j+kk]*scale[k+kk];
       }
    }


    cholesky(CovI, iChl, NI);
    cholesky(CovE, eChl, NE);
    
    free(scale);
    free_double_matrix(CovI,NI);
    free_double_matrix(CovE,NE);
    free_double_matrix(FishE,NE);
    
}

void de_jump(double *paramsx, double *paramsy, double **history, int m, int d, gsl_rng *r)
{
    int i, j, k;
    double alpha, beta;
    
    // pick two points from the history
    i = (int)((double)(m)*gsl_rng_uniform(r));
    k = 0;
    do
    {
        j = (int)((double)(m)*gsl_rng_uniform(r));
        k++;
    } while(i==j);
    
   // printf("%d\n", k);
    
    alpha = 1.0;
    beta = gsl_rng_uniform(r);
    if(beta < 0.9) alpha = gsl_ran_gaussian(r,0.5);
    
    beta = gsl_rng_uniform(r);
    if(beta > 0.5)
    {  // all parameters
    for(k=0; k< d; k++) paramsy[k] = paramsx[k]+alpha*(history[i][k]-history[j][k]);
    }
    else
    { // just the intrinsic (plus phase)
      for(k=0; k< NI; k++) paramsy[k] = paramsx[k]+alpha*(history[i][k]-history[j][k]);
      for(k=NI; k< d; k++) paramsy[k] = paramsx[k];
    }
    
}

void FisherEvec(double **fish, double *ej, double **ev, int d)
{
    int i, j, ecc, sc;
    double x, maxc;
 
    ecc = 0;
    for (i = 0 ; i < d ; i++) if(fabs(fish[i][i]) < 1.0e-16) ecc = 1;
    
    if(ecc == 0)
    {
        
        gsl_matrix *m = gsl_matrix_alloc (d, d);
        
        for (i = 0 ; i < d ; i++)
        {
            for (j = 0 ; j < d ; j++)
            {
                gsl_matrix_set(m, i, j, fish[i][j]);
            }
        }
        
        
        gsl_vector *eval = gsl_vector_alloc (d);
        gsl_matrix *evec = gsl_matrix_alloc (d, d);
        
        gsl_eigen_symmv_workspace * w =
        gsl_eigen_symmv_alloc (d);
        
        ecc = gsl_eigen_symmv (m, eval, evec, w);
        
        gsl_eigen_symmv_free (w);
        
        
        sc = gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
        
        for (i = 0; i < d; i++)
        {
            ej[i] = gsl_vector_get (eval, i);
            
           // printf("eigenvalue = %g\n", ej[i]);
            for (j = 0 ; j < d ; j++)
            {
                ev[i][j] = gsl_matrix_get(evec, j, i);
               // printf("%f ", ev[i][j]);
            }
           // printf("\n");
            
        }
        
        for (i = 0; i < d; i++)
        {
            // make sure no eigenvalue is too small
            //if(ej[i] < 100.0) ej[i] = 100.0;
            // turn into 1-sigma jump amplitudes
            ej[i] = 1.0/sqrt(fabs(ej[i]));
            //printf("jump %d = %g\n", i, ej[i]);
        }
        
        gsl_matrix_free (m);
        gsl_vector_free (eval);
        gsl_matrix_free (evec);
        
    }
    else
    {
        for (i = 0; i < d; i++)
        {
            ej[i] = 10000.0;
            for (j = 0 ; j < d ; j++)
            {
                ev[i][j] = 0.0;
                if(i==j) ev[i][j] = 1.0;
            }
        }
        
    }
    
    return;
    
}

void FisherEvecSVD(double **fish, double *ej, double **ev, int d)
{
    int i, j, ecc, sc;
    double x, maxc;
 
    ecc = 0;
    for (i = 0 ; i < d ; i++) if(fabs(fish[i][i]) < 1.0e-16) ecc = 1;
    
    if(ecc == 0)
    {
        
        gsl_matrix *m = gsl_matrix_alloc (d, d);
        
        for (i = 0 ; i < d ; i++)
        {
            for (j = 0 ; j < d ; j++)
            {
                gsl_matrix_set(m, i, j, fish[i][j]);
            }
        }
        
        
        gsl_vector *S = gsl_vector_alloc (d);
        gsl_matrix *V = gsl_matrix_alloc (d, d);
        
        gsl_linalg_SV_decomp_jacobi (m, V, S);
        
        
        for (i = 0; i < d; i++)
        {
            ej[i] = gsl_vector_get (S, i);
            
           // printf("eigenvalue = %g\n", ej[i]);
            for (j = 0 ; j < d ; j++)
            {
                ev[i][j] = gsl_matrix_get(m, j, i);
               // printf("%f ", ev[i][j]);
            }
           // printf("\n");
            
        }
        
        for (i = 0; i < d; i++)
        {
            // make sure no eigenvalue is too small
            //if(ej[i] < 100.0) ej[i] = 100.0;
            // turn into 1-sigma jump amplitudes
            //ej[i] = 1.0/sqrt(fabs(ej[i]));
            //printf("jump %d = %g\n", i, ej[i]);
        }
        
        gsl_matrix_free (m);
        gsl_vector_free (S);
        gsl_matrix_free (V);
        
    }
    else
    {
        for (i = 0; i < d; i++)
        {
            ej[i] = 10000.0;
            for (j = 0 ; j < d ; j++)
            {
                ev[i][j] = 0.0;
                if(i==j) ev[i][j] = 1.0;
            }
        }
        
    }
    
    return;
    
}

void FisherEvecSplit(double **fish, double *ej, double **ev, int d)
{
    int i, j, ecc, sc;
    double x, maxc;
    int NIN = 6;
    int NEX = 7;
    double **cov;
    
    
    ecc = 0;
    for (i = 0 ; i < d ; i++) if(fabs(fish[i][i]) < 1.0e-16) ecc = 1;
    
    
    
    if(ecc == 0)
    {
        
        for (i = 0; i < d; i++)
        {
            for (j = 0 ; j < d ; j++)
            {
                ev[i][j] = 0.0;
            }
        }
        
        cov = double_matrix(d,d);
        
        Inverse(fish,cov,d);
        
        // start with intrinsic
        
        gsl_matrix *m = gsl_matrix_alloc (NIN, NIN);
        
        for (i = 0 ; i < NIN ; i++)
        {
            for (j = 0 ; j < NIN ; j++)
            {
                gsl_matrix_set(m, i, j, cov[i][j]);
            }
        }
        
        
        gsl_vector *eval = gsl_vector_alloc (NIN);
        gsl_matrix *evec = gsl_matrix_alloc (NIN, NIN);
        
        gsl_eigen_symmv_workspace * w =
        gsl_eigen_symmv_alloc (NIN);
        
        ecc = gsl_eigen_symmv (m, eval, evec, w);
        
        gsl_eigen_symmv_free (w);
        
        
        sc = gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
        
        for (i = 0; i < NIN-1; i++)
        {
            ej[i] = gsl_vector_get (eval, i);
            
            // printf("eigenvalue = %g\n", ej[i]);
            for (j = 0 ; j < NIN ; j++)
            {
                ev[i][j] = gsl_matrix_get(evec, j, i);
                // printf("%f ", ev[i][j]);
            }
            // printf("\n");
            
        }
        
        for (i = 0; i < NIN-1; i++)
        {
            // make sure no eigenvalue is too small
            //if(ej[i] < 100.0) ej[i] = 100.0;
            // turn into 1-sigma jump amplitudes
            ej[i] = sqrt(fabs(ej[i]));
            //printf("jump %d = %g\n", i, ej[i]);
        }
        
        gsl_matrix_free (m);
        gsl_vector_free (eval);
        gsl_matrix_free (evec);
        
        // now the extrinsic
        
        gsl_matrix *mx = gsl_matrix_alloc (NEX, NEX);
        
        for (i = 0 ; i < NEX ; i++)
        {
            for (j = 0 ; j < NEX ; j++)
            {
                gsl_matrix_set(mx, i, j, cov[i+(NP-NEX)][j+(NP-NEX)]);
            }
        }
        
        
        gsl_vector *evalx = gsl_vector_alloc (NEX);
        gsl_matrix *evecx = gsl_matrix_alloc (NEX, NEX);
        
        gsl_eigen_symmv_workspace * wx =
        gsl_eigen_symmv_alloc (NEX);
        
        ecc = gsl_eigen_symmv (mx, evalx, evecx, wx);
        
        gsl_eigen_symmv_free (wx);
        
        
        sc = gsl_eigen_symmv_sort (evalx, evecx, GSL_EIGEN_SORT_ABS_ASC);
        
        for (i = 0; i < NEX-1; i++)
        {
            ej[i+NIN-1] = gsl_vector_get (evalx, i);
            
            // printf("eigenvalue = %g\n", ej[i]);
            for (j = 0 ; j < NEX ; j++)
            {
                ev[i+NIN-1][j+NP-NEX] = gsl_matrix_get(evecx, j, i);
                // printf("%f ", ev[i][j]);
            }
            // printf("\n");
            
        }
        
        for (i = 0; i < NEX; i++)
        {
            // make sure no eigenvalue is too small
            //if(ej[i] < 100.0) ej[i] = 100.0;
            // turn into 1-sigma jump amplitudes
            ej[i+NIN-1] = sqrt(fabs(ej[i]));
            //printf("jump %d = %g\n", i, ej[i]);
        }
        
        gsl_matrix_free (mx);
        gsl_vector_free (evalx);
        gsl_matrix_free (evecx);

        free_double_matrix(cov,d);
        
    }
    else
    {
        for (i = 0; i < d; i++)
        {
            ej[i] = 10000.0;
            for (j = 0 ; j < d ; j++)
            {
                ev[i][j] = 0.0;
                if(i==j) ev[i][j] = 1.0;
            }
        }
        
    }
    
    return;
    
}



void Inverse(double **M, double **IM, int d)
{
    int i, j;
    int s;
    double x, maxc;
    
    gsl_matrix *m = gsl_matrix_alloc (d, d);
    
    for (i = 0 ; i < d ; i++)
    {
        for (j = 0 ; j < d ; j++)
        {
            gsl_matrix_set(m, i, j, M[i][j]);
        }
    }
    
    gsl_permutation *p = gsl_permutation_alloc(d);
    
    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(m, p, &s);
    
    // Compute the  inverse of the LU decomposition
    gsl_matrix *inv = gsl_matrix_alloc(d, d);
    
    gsl_linalg_LU_invert(m, p, inv);
    
    gsl_permutation_free(p);
    
    
    for (i = 0; i < d; i++)
    {
        for (j = 0 ; j < d ; j++)
        {
            IM[i][j] = gsl_matrix_get(inv, j, i);
        }
    }
    
    
    gsl_matrix_free (inv);
    gsl_matrix_free (m);
    
    return;
    
}


void cholesky(double **A, double **C, int N)
{
    int *IPIV;
    int LWORK = N*N;
    int INFO;
    int i, j;
    double dx, dy;
    

    
    gsl_matrix *m = gsl_matrix_alloc (N, N);
    gsl_vector *s = gsl_vector_alloc (N);
    
    for (i = 0 ; i < N ; i++)
    {
        for (j = 0 ; j < N ; j++)
        {
            gsl_matrix_set(m, i, j, A[i][j]);
        }
    }
    
   /* gsl_linalg_cholesky_decomp1(m);
    for (i = 0; i < N; i++)
    {
        for (j = 0 ; j < N ; j++)
        {
            C[i][j] = gsl_matrix_get(m, i, j);
        }
    }
    */
    
    gsl_linalg_cholesky_decomp2(m, s);
    
    
    for (i = 0; i < N; i++)
    {
        for (j = 0 ; j < N ; j++)
        {
            C[i][j] = gsl_matrix_get(m, i, j)/gsl_vector_get(s, i);
        }
    }
    
    for (i = 0; i < N; i++)
    {
        for (j = i+1 ; j < N ; j++)
        {
            C[i][j] = 0.0;
        }
    }


    gsl_matrix_free(m);
    gsl_vector_free(s);
    
    return;
    
    
}



double det(double **A, int N)
{
    int *IPIV;
    int LWORK = N*N;
    int INFO;
    int i, j;
    double dx, dy;
    int s;
    
    gsl_permutation *p = gsl_permutation_alloc(N);
    
    gsl_matrix *m = gsl_matrix_alloc (N, N);
    
    for (i = 0 ; i < N ; i++)
    {
        for (j = 0 ; j < N ; j++)
        {
            gsl_matrix_set(m, i, j, A[i][j]);
        }
    }
    
    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(m, p, &s);
    
    dx = 1.0;
    for (i = 0; i < N; i++) dx *= gsl_matrix_get (m, i, i);
    dx = fabs(dx);
    
    //returns the absolute value of the determinant.
    
    gsl_permutation_free(p);
    gsl_matrix_free(m);
    
    
    return dx;
    
    
}
