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

double logL(int M, int Nseg, double *cc, double *ss, double **CT, double **ST, double *SMT);

#define TPI 6.2831853071795865
#define year 3.15581498e7

//gcc -o unique unique.c -lgsl


int main(int argc,char **argv)
{
    char filename[50];
    int i, j, k, id, mc;
    int ii, jj, kk, ll, flag;
    double x, y, z, f, dt;
    double Tobs, Tseg;
    int N, Nseg, Nch;
    double **params, **paramsU;
    double *SNR, *SNRU, *Mc, *McU;
    int *NS;
    int NST, NSU;
    FILE *in;
    FILE *out;
    
    Nch = 2;
    Nseg = 12;
    N = 524288;
    dt = 5.0;
    Tseg = dt*(double)(N);
    
    
    NS = int_vector(Nseg);
    
    NST = 0;
    for(j=0; j< Nseg; j++)
        {
         i = 0;
         sprintf(filename, "skymax_%d_%d.dat", j, i);
          while ((in = fopen(filename, "r")))
          {
            i++;
            fclose(in);
            sprintf(filename, "skymax_%d_%d.dat", j, i);
          }
           //printf("%d %d\n", j, i);
            NS[j] = i;
            NST += i;
        }
    
    params = double_matrix(NST,NParams);
    SNR = double_vector(NST);
    Mc = double_vector(NST);
    McU = double_vector(NST);
    SNRU = double_vector(NST);
    
    ii = 0;
    for(j=0; j< Nseg; j++)
        {
        for(i=0; i< NS[j]; i++)
          {
            sprintf(filename, "skymax_%d_%d.dat", j, i);
              in = fopen(filename,"r");
              fscanf(in,"%lf", &x);
              SNR[ii] = sqrt(2.0*x);
              for(k=0; k< NParams; k++) fscanf(in,"%lf", &params[ii][k]);
              fclose(in);
              Mc[ii] = pow((params[ii][0]*params[ii][1]),3.0/5.0)/pow((params[ii][0]+params[ii][1]),1.0/5.0);
              //printf("%f %e %e\n", SNR[ii], params[ii][5], (double)(j+1)*Tseg);
              ii++;
          }
        }
    
    printf("\n");
   
      paramsU = double_matrix(NST,NParams);
      ii = 0;
      jj = 0;
       for(j=0; j< Nseg-1; j++)
           {
           for(i=0; i< NS[j]; i++)
             {
                //printf("%d %d\n", ii, jj);
                // only keep sources that merge in a given segment
                // only keep sources with large enough SNR
                // skip sources that are in the boundary regions (impacted by Tukey filter)
            if(params[ii][5] < (double)(j+1)*Tseg && SNR[ii] > 12.0 && fabs(params[ii][5]-(double)(j)*Tseg) > 1.0e5 && fabs(params[ii][5]-(double)(j+1)*Tseg) > 1.0e5 && fabs(params[ii][5]-(double)(j+2)*Tseg) > 1.0e5)
               {
                   flag = 0;
                   
                 for(kk=0; kk< jj; kk++)
                 {
                     y = fabs(params[ii][5]-paramsU[kk][5])/3600.0;
                     z = fabs(Mc[ii]-McU[kk])/(Mc[ii]+McU[kk]);
                    // if(y < 50.0 && z < 0.2) printf("%f %f %e %e %e %e\n", y, z, params[ii][5], paramsU[kk][5], Mc[ii], McU[kk]);
                    if(y < 1.0 && z < 0.2 && SNR[ii] < 20.0)
                    {
                        flag = 1;
                        if(SNR[ii] > SNRU[kk]) // we will take the higher SNR copy
                        {
                            SNRU[kk] = SNR[ii];
                            McU[kk] = Mc[ii];
                            for(k=0; k< NParams; k++) paramsU[kk][k] = params[ii][k];
                        }
                    }
                 }
                if(flag == 0)
                {
                for(k=0; k< NParams; k++) paramsU[jj][k] = params[ii][k];
                McU[jj] = Mc[ii];
                SNRU[jj] = SNR[ii];
                jj++;
                }
               }
              ii++;
             }
           }
    
    // for the last segment keep all the sources
    j = Nseg-1;
    for(i=0; i< NS[j]; i++)
    {
       if(SNR[ii] > 12.0)
       {
         for(k=0; k< NParams; k++) paramsU[jj][k] = params[ii][k];
        // printf("%d %e\n", i, params[ii][5]);
         McU[jj] = Mc[ii];
         SNRU[jj] = SNR[ii];
         jj++;
       }
         ii++;
    }

    NSU = jj;
    
    printf("Total sources %d\nUnique sources %d\n", NST, NSU);
    
    out = fopen("search_sources.dat","w");
    for(i=0; i< NSU; i++)
    {
       fprintf(out,"%.15e %.15e ", SNRU[i], McU[i]);
       for(k=0; k< NParams; k++) fprintf(out,"%.15e ", paramsU[i][k]);
       fprintf(out,"\n");
    }
    fclose(out);
    
    printf("\n");
    
    out = fopen("source_info.dat","w");
    for(i=0; i< NSU; i++)
    {
       printf("Source %d merges in segment %d\n", i, (int)(paramsU[i][5]/Tseg));
       fprintf(out,"Source %d merges in segment %d\n", i, (int)(paramsU[i][5]/Tseg));
    }
    fclose(out);
    
    printf("\n");
 
    
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

    
