/*******************************************************************************************
 
 Copyright (c) 2021 Neil Cornish & Tyson Littenberg
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 **********************************************************************************************/

#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "mbh.h"

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


void Inverse(double **M, double **IM, int d)
{
    int i, j;
    int s;
    
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



double det(double **A, int N)
{
    int i, j;
    double dx;
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

void get_component_masses(double *params, int flag, double *m1, double *m2)
{
    double Mc, Mtot, eta, dm;
    switch(flag)
    {
        case 0:
            *m1 = params[0];
            *m2 = params[1];
            break;
        case 1:
            *m1 = exp(params[0]);
            *m2 = exp(params[1]);
            break;
        case 2:
            Mc = exp(params[0]);
            Mtot = exp(params[1]);
            eta = pow(Mc/Mtot, 5./3.);
            dm = 0.0;
            if(eta < 0.25) dm = sqrt(1. - 4.*eta);
            *m1 = 0.5*Mtot*(1.+dm);
            *m2 = 0.5*Mtot*(1.-dm);
            break;
        default:
            break;
    }
}

void print_mbh_chain_file(struct MBH_Data *dat, struct Het *het, int *who, double **paramx, double *logLx, double **sx, int ll, int mc, FILE *chain)
{
    int i;
    int k;
    int q;
    double m1;
    double m2;
    double thetaL;
    double phiL;
    double DL;
    double x;
    
    for(k = 0; k < NCC; k++)
    {
        q = who[k];
        
        get_component_masses(paramx[q], ll, &m1, &m2);
        
        lisaskyloc(dat->Tend, paramx[q], &thetaL, &phiL);
        
        DL = exp(paramx[q][6]);
        
        x = logLx[q];
        
        if(nflag == 1)  // only want to record the reduced log likelihood
        {
            for (i=0; i< dat->Nch; i++)
            {
                x += 0.5*het->DD[i]/sx[q][i];
                x += (double)(het->MM-het->MN)*log(sx[q][i]);
            }
        }
        
        fprintf(chain,"%d %.12e %.12e %.12e ", mc+k, x, m1, m2);
        for(i = 2; i < NParams; i++)
        {
            if(i == 6)
            {
                fprintf(chain, "%.12e ", DL);
            }
            else
            {
                fprintf(chain, "%.15e ", paramx[q][i]);
            }
        }
        fprintf(chain,"%d %f %f %f %f\n", q, thetaL, phiL, sx[q][0], sx[q][1]);
    }
}

void set_mbh_priors(struct MBH_Data *dat, int massFlag, double *min, double *max)
{
    // ll = 0  [0] Mass1  [1] Mass2
    // ll = 1  [0] ln Mass1  [1] ln Mass2
    // ll = 2  [0] ln Mc  [1] ln Mtot
    // [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln distance
    // [7] cos EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] cos inclination

    switch (massFlag)
    {
        case 0:
            max[0] = 5.0e8;
            max[1] = 5.0e8;
            min[0] = 1.0e3;
            min[1] = 1.0e3;
            break;
        case 1:
            max[0] = log(5.0e8);
            max[1] = log(5.0e8);
            min[0] = log(1.0e3);
            min[1] = log(1.0e3);
            break;
        case 2:
            max[0] = log(0.44*5.0e8);
            max[1] = log(5.0e8);
            min[0] = log(1.0e2);
            min[1] = log(1.0e3);
            break;
        default:
            fprintf(stderr,"massFlag must be [0,1,2]\n");
            exit(1);
            break;
    }
    
    
    
    max[2] = 0.999;
    max[3] = 0.999;
    max[4] = PI;
    max[5] = 2.0*dat->Tend;
    max[6] = log(1.0e3);
    max[7] = 1.0;
    max[8] = 2.0*PI;
    max[9] = PI;
    max[10] = 1.0;
    
    
    min[2] = -0.999;
    min[3] = -0.999;
    min[4] = 0.0;
    min[5] = 1.01*dat->Tstart;
    min[6] = log(0.1);
    min[7] = -1.0;
    min[8] = 0.0;
    min[9] = 0.0;
    min[10] = -1.0;
}

