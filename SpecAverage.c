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

double logL(int M, int Nseg, double *cc, double *ss, double **CT, double **ST, double *SMT);

#define TPI 6.2831853071795865
#define year 3.15581498e7

//gcc -o SpecAverage SpecAverage.c -lgsl


int main(int argc,char **argv)
{
    char filename[50];
    int i, j, k, id, mc;
    int ii, jj, kk, ll;
    double x, y, z, u, v, f, dt;
    double Tobs, Tseg, Ttot;
    double fac;
    int N, Nseg, Nch, Ntot;
    double ***SM, ***SN;
    double **SNA, **SMA;
    double *fa;
    FILE *in;
    FILE *out;
    
    Nch = 2;
    Nseg = 12;
    N = 524288;
    Ntot = 8388608;   // padded data length for entire data set
    dt = 5.0;
    Tseg = dt*(double)(N);
    Ttot = dt*(double)(Ntot);
    
    fa = double_vector(N/2);
    
    SM = double_tensor(Nseg,Nch,N/2);
    SN = double_tensor(Nseg,Nch,N/2);
    
    SNA = double_matrix(Nch,N/2);
    SMA = double_matrix(Nch,N/2);
    
    for(id=0; id < Nch; id++)
         {
           for(i=0; i< N/2; i++)
           {
             SMA[id][i] = 0.0;
             SNA[id][i] = 0.0;
           }
         }
    
    for(j=0; j< Nseg; j++)
        {
     for(id=0; id < Nch; id++)
      {
        sprintf(filename, "specfit_%d_%d.dat", id, j);
        in = fopen(filename,"r");
        for(i=0; i< N/2; i++)
        {
          fscanf(in,"%lf%lf%lf%lf\n", &fa[i], &SM[j][id][i], &x, &SN[j][id][i]);
            SMA[id][i] += SM[j][id][i];
            SNA[id][i] += SN[j][id][i];
        }
        fclose(in);
      }
    }
    
    for(i=0; i< N/2; i++) fa[i] = (double)(i)/Tseg;
    
       for(id=0; id < Nch; id++)
            {
              for(i=0; i< N/2; i++)
              {
                SMA[id][i] /= (double)(Nseg);
                SNA[id][i] /= (double)(Nseg);
              }
            }
    
     
    for(id=0; id < Nch; id++)
      {
        sprintf(filename, "sav_%d.dat", id);
        out = fopen(filename,"w");
        for(i=0; i< N/2; i++)
         {
             fprintf(out,"%.14e %.14e %.14e\n", fa[i], SMA[id][i], SNA[id][i]);
         }
          fclose(out);
      }
    
    // 
    
    for(id=0; id < Nch; id++)
    {
      sprintf(filename, "specav_%d.dat", id);
      out = fopen(filename,"w");
      for(i=0; i< Ntot/2; i++)
      {
          f = (double)(i)/Ttot;
    
          ii = (int)(floor(f*Tseg)); // find lower bin
          
          //printf("%d %d\n", i, ii);
          
          x = 1.0;
          y = 1.0;
          
          if(ii+1 < N/2)
          {
          u = (fa[ii+1]-f)*Tseg;
          v = (f-fa[ii])*Tseg;
          x = SMA[id][ii]*u + SMA[id][ii+1]*v;
          y = SNA[id][ii]*u + SNA[id][ii+1]*v;
          }
          
          // The factor of 0.75 corrects for the ratio of filled to empty
          // time samples 12/16 = 3/4
          fprintf(out,"%.14e %.14e %.14e\n", f, 0.75*x, 0.75*y);
      }
      fclose(out);
    }
    

    
}

double logL(int M, int Nseg, double *cc, double *ss, double **CT, double **ST, double *SMT)
{
    int i, j;
    double LL;
    double x, y;
    
    LL = 0.0;
    
    for(j=0; j< Nseg; j++)
      {
          x = 0.0;
           for(i=0; i< M; i++)
            {
                x += cc[i]*CT[i][j]+ss[i]*ST[i][j];
            }
          y = (SMT[j]-x)/0.005;
          LL -= y*y;
      }
    
    return LL;
    
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

    
