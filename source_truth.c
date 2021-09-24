/*******************************************************************************************

Copyright (c) 2021 Neil Cornish

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Constants.h"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>


int *int_vector(int N);
void free_int_vector(int *v);
double **double_matrix(int N, int M);
void free_double_matrix(double **m, int N);
int **int_matrix(int N, int M);
void free_int_matrix(int **m, int N);
double *double_vector(int N);
void free_double_vector(double *v);
double ***double_tensor(int N, int M, int L);
void free_double_tensor(double ***t, int N, int M);

double DL_fit(double z);
double z_DL(double DL);


//gcc -o source_truth source_truth.c -lgsl


int main(int argc,char **argv)
{
    char filename[50];
    int i, j, k, mc, N;
    double m1, m2, x, z;
    double chi1, chi2, elat, elong, DL, thetaL, phiL;
    FILE *in;
    FILE *out;

    
    if(argc<3)
      {
          printf("./source_truth file_in file_out\n");
          return 0;
      }
    
    in = fopen(argv[1],"r");
    out = fopen(argv[2],"w");
    
    //m1 m2 chi1 chi2 elat elong DL thetaL phiL
    
    fscanf(in,"%lf%lf%lf%lf%lf%lf%lf%lf%lf", &m1, &m2, &chi1, &chi2, &elat, &elong, &DL, &thetaL, &phiL);
    
    z = z_DL(DL*1000.0);
    
    printf("%f %e\n", z, DL);
         
    fprintf(out,"%e %e %e %e %e %e %e %e %e\n", m1/(1.0+z), m2/(1.0+z), chi1, chi2, elat, elong, DL, thetaL, phiL);
    
    fclose(in);
    fclose(out);

    
}

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

double z_DL(double DL)
{
    double zlow, zhigh;
    double zmid;
    double Dlow, Dhigh, Dmid;
    double u, v, w;
    int i;
    
    zlow = 1.0;
    zhigh = 1000.0;
    
    if(DL < 1.0)
    {
        zlow = 1.0e-10;
        zhigh = 3.0e-4;
    }
    else if (DL < 100.0)
    {
        zlow = 2.0e-4;
        zhigh = 5.0e-2;
    }
    else if (DL < 1000.0)
    {
        zlow = 1.0e-2;
        zhigh = 3.0e-1;
    }
    else if (DL < 10000.0)
    {
        zlow = 1.0e-1;
        zhigh = 2.0;
    }
    
    zmid = (zhigh+zlow)/2.0;
    
    Dlow = DL_fit(zlow);
    Dhigh = DL_fit(zhigh);
    Dmid = DL_fit(zmid);
    
    for (i=0; i< 30; i++)
    {
        
        u = Dlow-DL;
        v = Dmid-DL;
        
        if(u*v < 0.0)
        {
            zhigh = zmid;
            Dhigh = Dmid;
        }
        else
        {
            zlow = zmid;
            Dlow = Dmid;
        }
        zmid = (zhigh+zlow)/2.0;
        Dmid = DL_fit(zmid);
        
        
    }
    
    //  printf("%e %e\n", Dmid-DL, zmid);
    
    return(zmid);
    
    
}


double DL_fit(double z)
{
    double D;
    // Planck values https://arxiv.org/abs/1807.06209
    double Om = 0.315;
    double H0 = 67.4;
    double RH, x1, x2, Om1, Om2;
    
    // See http://arxiv.org/pdf/1111.6396v1.pdf
    
    RH = clight/H0*(1.0e-3*MPC);
    
    x1 = (1.0-Om)/(Om);
    
    x2 = (1.0-Om)/(Om*pow((1.0+z),3.0));
    
    Om1 = (1.0+1.32*x1+0.4415*x1*x1+0.02656*x1*x1*x1)/(1.0+1.392*x1+0.5121*x1*x1+0.03944*x1*x1*x1);
    
    Om2 = (1.0+1.32*x2+0.4415*x2*x2+0.02656*x2*x2*x2)/(1.0+1.392*x2+0.5121*x2*x2+0.03944*x2*x2*x2);
    
    D = 2.0*RH*(1.0+z)/sqrt(Om)*(Om1-Om2/sqrt(1.0+z));
    
    return(D/MPC);
}


    
