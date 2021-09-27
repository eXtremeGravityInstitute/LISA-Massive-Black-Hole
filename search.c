/*******************************************************************************************

Copyright (c) 2019 Neil Cornish

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


#include "IMRPhenomD.h"

#include <ctype.h>
#include <time.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#include "Constants_search.h"
#include "Detector.h"
#include "Declarations_search.h"

#define NR_END 1
#define FREE_ARG char*

#ifndef _OPENMP
#define omp ignore
#endif

// OSX
// clang -Xpreprocessor -fopenmp -lomp -w -o search search.c IMRPhenomD_internals.c IMRPhenomD.c -lgsl -lgslcblas  -lm

// Linux
// gcc -std=gnu99 -fopenmp -w -o search search.c IMRPhenomD_internals.c IMRPhenomD.c -lgsl -lgslcblas  -lm

//##############################################
//MT modifications

gsl_rng **rvec;
//##############################################

int main(int argc, char *argv[])
{
    double phi0;
    double deltaF;
    double m1_SI;
    double m2_SI;
    double chi1;
    double chi2;
    double distance;
    double tmerger;
    double zred;
    double Mtot, Mc, m1, m2, s1, s2;
    double f, cp, sp, p, tshift;
    double cv, sv, v, A, dtm, tlength, t0;
    double logL;
    int ret, seg;
    double *params, *pref, *zv;
    double *pmax;
    double HH, HD, DD;
    double thetaL, phiL, betaE, lambdaE, psi, iota;
    double fstart, fstop;
    double Tzero, beta;
    double **paramx;
    int ll, k, rep;
    int st, se;
    char filename[1024];

    FILE *in;
    FILE *out;
    
    double x, dt, Tobs;
    int i, j, N;
    
    const gsl_rng_type * P;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    
    P = gsl_rng_default;
    r = gsl_rng_alloc (P);
    
    //##############################################
    //open MP modifications
    omp_set_num_threads(NC);
    rvec = (gsl_rng **)malloc(sizeof(gsl_rng *) * (NC));
    for(i = 0 ; i< NC; i++){
        rvec[i] = gsl_rng_alloc(P);
        gsl_rng_set(rvec[i] , i);
    }
    //##############################################
    
    dt = 5.0;
    N = 524288;
    Tobs = (double)(N)*dt;

    pmax = double_vector(NP);
    params = (double*)malloc(sizeof(double)*NP);
    
    double *AC, *EC, *TC, *SN;
    double *AS, *ES;
    double *SA, *SE;
    double *SAS, *SES;
    AC = (double*)malloc(sizeof(double)* (N));
    EC = (double*)malloc(sizeof(double)* (N));
    AS = (double*)malloc(sizeof(double)* (N));
    ES = (double*)malloc(sizeof(double)* (N));
    TC = (double*)malloc(sizeof(double)* (N));
    SN = (double*)malloc(sizeof(double)* (N));  // goes out further in frequency than we need.
    SA = (double*)malloc(sizeof(double)* (N/2));
    SE = (double*)malloc(sizeof(double)* (N/2));
    SAS = (double*)malloc(sizeof(double)* (N/2));
    SES = (double*)malloc(sizeof(double)* (N/2));
    
    paramx = double_matrix(NC,NP);
    
    if(argc<3)
    {
        printf("./search start_seg end_seg\n");
        printf("segment numbers run from 0 to 22\n");
        return 0;
    }
    
    st = atoi(argv[1]);
    se = atoi(argv[2]);
    
   for(seg=st; seg<= se; seg++)
    {
        
    Tzero = (double)(seg)*Tobs/2.0;
        
    printf("%.0f %.0f %d\n", Tobs, Tzero, N);
    
    // Read in FFTed LDC data and A,E PSD from segmentAET.c
    sprintf(filename, "AET_seg%d_f.dat", seg);
    in = fopen(filename,"r");
    for(i=0; i< N; i++)
    {
        fscanf(in,"%lf%lf%lf%lf\n", &f, &AC[i], &EC[i], &TC[i]);
    }
    fclose(in);
    
    sprintf(filename, "specfit_0_%d.dat", seg);
    in = fopen(filename,"r");
    for(i=0; i< N/2; i++)
    {
        fscanf(in,"%lf%lf%lf%lf\n", &f, &SAS[i], &x, &SA[i]);
    }
    fclose(in);
    
    sprintf(filename, "specfit_1_%d.dat", seg);
    in = fopen(filename,"r");
    for(i=0; i< N/2; i++)
    {
        fscanf(in,"%lf%lf%lf%lf\n", &f, &SES[i], &x, &SE[i]);
    }
    fclose(in);
    
    // used in FisherMax
     for(i=0; i< N/2; i++)
       {
           SN[i] = 1.0/(1.0/SAS[i]+1.0/SES[i]);
       }
    
    
    rep = 0;
    
    
    sprintf(filename, "Spec_%d_%d.dat", seg, rep);
    out = fopen(filename,"w");
    for(i=1; i< N/2; i++)
     {
           fprintf(out,"%e %e %e\n", (double)(i)/Tobs, 2.0*(AC[i]*AC[i]+AC[N-i]*AC[N-i]), 2.0*(EC[i]*EC[i]+EC[N-i]*EC[N-i]));
     }
    fclose(out);
    
       // Used for testing against the PTMCMC code
       /*sprintf(filename, "skymax_%d_%d.dat", seg, rep);
       in = fopen(filename,"r");
       fscanf(in,"%lf", &x);
       for(i=0; i< NP; i++) fscanf(in,"%lf", &params[i]);
       fclose(in);
       logL =  Likelihood(0, params, N, AC, EC, SA, SE, Tobs, Tzero);
       printf("log L = %.12e\n", logL);
       return(0);
       */
        
        
    
    sprintf(filename, "Qscan_0_%d_%d.dat", seg, rep);
    qscanf(filename, AC, SA, Tobs, N);
    sprintf(filename, "Qscan_1_%d_%d.dat", seg, rep);
    qscanf(filename, EC, SE, Tobs, N);
     
    
    do
    {
        
      search(pmax, paramx, AC, EC, SN, SA, SE, Tobs, seg, N, rep);
    
       printf("\n Sky search \n");
        
        // Parameter order
        
    // [0] ln(Mc)  [1] ln(Mt)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
    // [7] cos(EclipticCoLatitude), [8] EclipticLongitude  [9] polarization, [10] inclination
        
     searchsky(pmax, paramx, AC, EC, SAS, SES, SA, SE, Tobs, seg, N, rep);
        
     logL = Likelihood(mtype, pmax, N, AC, EC, SA, SE, Tobs, Tzero);
        
    // Remove the best fit signal
     ResponseFreq(mtype, pmax, N, AS, ES, Tobs, Tzero);
     for(i=0; i< N; i++)
     {
        AC[i] -= AS[i];
        EC[i] -= ES[i];
     }
        
     rep++;
        
     sprintf(filename, "Qscan_0_%d_%d.dat", seg, rep);
     qscanf(filename, AC, SA, Tobs, N);
     sprintf(filename, "Qscan_1_%d_%d.dat", seg, rep);
     qscanf(filename, EC, SE, Tobs, N);
        
     sprintf(filename, "Spec_%d_%d.dat", seg, rep);
     out = fopen(filename,"w");
     for(i=1; i< N/2; i++)
      {
            fprintf(out,"%e %e %e\n", (double)(i)/Tobs, 2.0*(AC[i]*AC[i]+AC[N-i]*AC[N-i]), 2.0*(EC[i]*EC[i]+EC[N-i]*EC[N-i]));
      }
     fclose(out);
    
        
    }while(logL > 80.0 && rep < 5);
    
    // store the residual
    sprintf(filename, "Res_%d.dat", seg);
    out = fopen(filename,"w");
    for(i=0; i< N; i++)
     {
           fprintf(out,"%e %e\n", AC[i], EC[i]);
     }
    fclose(out);
        
    }
    
    
    return 1;
    
    
}


void search(double *pmax,  double **paramx, double *AC, double *EC, double *SN, double *SA, double *SE, double Tobs, int seg, int N, int rep)
{
    int i, j, k, jj, mc, hold;
    int *mh;
    double x, y, logLmax;
    double *logLx;
    double *min, *max;
    int *scount, *sacc;
    int **av, **cv;
    int ll=mtype;
    int *who;
    double *heat;
    int q, flag;
    double m1, m2, DL, Tzero, Mc, Mtot, dm, eta, leta;
    double alpha, beta, fny;
    double ***history, ***evec, ***Fisher, **eval;
    char filename[1024];
    
    FILE *chain;
    FILE *levels;
    
    const gsl_rng_type * P;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    
    P = gsl_rng_default;
    r = gsl_rng_alloc (P);
    
    mh = int_vector(NC);
    sacc = int_vector(NC);
    scount = int_vector(NC);
    who = int_vector(NC);
    max = double_vector(NP);
    min = double_vector(NP);
    heat = double_vector(NC);
    logLx = double_vector(NC);
    av = int_matrix(5,NC);
    cv = int_matrix(5,NC);
    
   
    Fisher = double_tensor(NC,NV,NV);
    history = double_tensor(NC,NH,NV);
    eval = double_matrix(NC,NV);
    evec = double_tensor(NC,NV,NV);
    
    Tzero = (double)(seg)*Tobs/2.0;
    
    for (i=0; i< NC; i++) who[i] = i;
    heat[0] = 1.0;
    for (i=1; i< NC; i++) heat[i] = heat[i-1]*1.1;
    
    
    for(k=0; k < NC; k++) mh[k] = 0;
    
    for(j = 0; j < 5; j++)
    {
        for(k=0; k < NC; k++)
        {
            av[j][k] = 0;
            cv[j][k] = 1;
        }
    }
    
    for(k=0; k < NC; k++)
    {
        sacc[k] = 0;
        scount[k] = 1;
    }
    
    max[0] = log(1.0e9);
    max[1] = log(2.3e9);  // (should be not more than 2.3 times Mcmax)
    max[2] = 0.95;
    max[3] = 0.95;
    max[4] = PI;
    max[5] = Tzero+2.0*Tobs-buffer*Tobs;  // merger could be up to 2 segments away
    max[6] = log(1.0e3);
    max[7] = 10.0;
    max[8] = 2.0*PI;
    max[9] = PI;
    max[10] = 1.0;
   
    min[0] = log(1.0e4);
    min[1] = log(1.0e4);
    min[2] = -0.95;
    min[3] = -0.95;
    min[4] = 0.0;
    min[5] = Tzero+buffer*Tobs;
    min[6] = log(0.1);
    min[7] = 0.1;
    min[8] = 0.0;
    min[9] = 0.0;
    min[10] = -1.0;
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
    // [7] E/A amplitude [8] PhaseA-PhaseE
    

    for (i=0; i< NC; i++)
       {
            
            for (j=0; j< NP; j++) paramx[i][j] = min[j]+(max[j]-min[j])*gsl_rng_uniform(r);
                    
             if(ll == 2)  // have to make sure the mass ratio is ok when using chirp,total mass
              {
                  
               do
               {
                   
               flag = 0;
               for (j=0; j< 2; j++) paramx[i][j] = min[j]+(max[j]-min[j])*gsl_rng_uniform(r);
               // eta cannot exceed 0.25
               leta = (5.0/3.0)*(paramx[i][0]-paramx[i][1]);
               if(leta > log(0.25)) paramx[i][1] = paramx[i][0] + 0.856; // a little below equal mass
               for(j = 0; j < 2; j++)
               {
                if(paramx[i][j] > max[j] || paramx[i][j] < min[j]) flag = 1;
               }
                   
                }while(flag == 1);
              }
           
           logLx[i] = log_likelihood_max_dual(ll, AC, EC, paramx[i], SA, SE, N, Tobs, Tzero);
           
       }
                
    
    sprintf(filename, "search_%d_%d.dat", seg, rep);
    chain = fopen(filename,"w");
    sprintf(filename, "likes_%d_%d.dat", seg, rep);
    levels = fopen(filename,"w");
    
    logLmax = -1.0e10;
    
    for(mc = 0; mc < MS; mc++)
    {
        
        if((mc+1)%(MS/8) == 0) // adapt the temperatures
        {
         // effective SNR of 7 for hottest chain
         alpha = pow((2.0*logLmax/50.0),1.0/(double)(NC));
         if(alpha > 1.3) alpha = 1.3;
         if(alpha < 1.0) alpha = 1.05;
         printf("reset alpha = %f\n", alpha);
         heat[0] = 1.0;
         for (i=1; i< NC; i++) heat[i] = heat[i-1]*alpha;
        }
        
        if(mc%50==0)
        {
        #pragma omp parallel for
        for(i = 0; i < NC; i++)
        {
            // printf("%d %e %e %e\n", i, exp(paramx[i][0]), exp(paramx[i][1]), paramx[i][5]);
            FisherMax(ll, paramx[i], Fisher[i], SN, Tobs, Tzero, N);
            FisherEvec(Fisher[i], eval[i], evec[i], NV);
        }
        }
    
      alpha = gsl_rng_uniform(r);
        
      if((NC > 1) && (alpha < 0.3))  // decide if we are doing a MCMC update of all the chains or a PT swap
       {
        // chain swap
       
        for(i = 0; i < NC/2; i++)
        {
        alpha = (double)(NC-1)*gsl_rng_uniform(r);
        j = (int)(alpha);
        beta = ((logLx[who[j]]-logLx[who[j+1]])/heat[j+1] - (logLx[who[j]]-logLx[who[j+1]])/heat[j]);
        alpha = log(gsl_rng_uniform(r));
        if(beta > alpha)
        {
            hold = who[j];
            who[j] = who[j+1];
            who[j+1] = hold;
            sacc[j]++;
        }
        
         scount[j]++;
        }
        
       }
      else      // MCMC update
       {
           #pragma omp parallel for
           for(k=0; k < NC; k++)
            {
             update(mc, k, ll, logLx, paramx, eval, evec, history, min, max, who, heat, av, cv, AC, EC, SA, SE, Tobs, Tzero, N, rvec[k]);
            }
       }
        
        for(k=0; k < NC; k++)
        {
          if(logLx[k] > logLmax)
          {
              logLmax = logLx[k];
              for(i = 0; i < NP; i++) pmax[i] = paramx[k][i];
          }
        }
        
        // add to the history file
           if(mc%2 == 0)
            {
             for(k=0; k < NC; k++)
             {
                q = who[k];
                i = mh[k]%NH;
                // the history file is kept for each temperature
                for(j=0; j<NV; j++) history[k][i][j] = paramx[q][j];
                mh[k]++;
             }
            }
        
        
        if(mc%50 == 0)
        {
            
            // clone the best solution into the worst
            x = -1.0e4;
            y = 1.0e20;
            j = jj = 0;
            for(k=0; k < NC; k++)
            {
              if(logLx[k] > x)
              {
                  x = logLx[k];
                  j = k;
              }
                
              if(logLx[k] < y)
              {
                 y = logLx[k];
                 jj = k;
              }
                
            }
            
            q = who[jj];
            logLx[q] = logLx[j];
            for(i = 0; i < NP; i++) paramx[q][i] = paramx[j][i];
            
        }
        
        if(mc%10 == 0)
        {
            
            q = who[0];
            
            if(ll == 0)
            {
             m1 = paramx[q][0];
             m2 = paramx[q][1];
             DL = paramx[q][6];
            }
            if(ll==1)
            {
             m1 = exp(paramx[q][0]);
             m2 = exp(paramx[q][1]);
             DL = exp(paramx[q][6]);
            }
            if(ll == 2)
            {
            DL = exp(paramx[q][6]);
            Mc = exp(paramx[q][0]);
            Mtot = exp(paramx[q][1]);
            eta = pow((Mc/Mtot), (5.0/3.0));
             if(eta > 0.25)
             {
                dm = 0.0;
             }
             else
             {
                dm = sqrt(1.0-4.0*eta);
             }
            m1 = Mtot*(1.0+dm)/2.0;
            m2 = Mtot*(1.0-dm)/2.0;
            }
            
               fprintf(chain,"%d %.12e %.12e %.12e ", mc, logLx[q], m1, m2);
               for(i = 2; i < NP; i++)
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
              fprintf(chain,"%d\n", q);
            
              fprintf(levels, "%d ", mc);
               for(k = 0; k < NC; k++)
                {
                fprintf(levels, "%.12e ", logLx[who[k]]);
                }
                 fprintf(levels, "\n");
          }
        
              if(mc%10 == 0)
                {
                    q = who[0];
                    
                    if(ll == 0)
                    {
                     m1 = paramx[q][0];
                     m2 = paramx[q][1];
                     DL = paramx[q][6];
                    }
                    if(ll==1)
                    {
                     m1 = exp(paramx[q][0]);
                     m2 = exp(paramx[q][1]);
                     DL = exp(paramx[q][6]);
                    }
                    if(ll == 2)
                    {
                    DL = exp(paramx[q][6]);
                    Mc = exp(paramx[q][0]);
                    Mtot = exp(paramx[q][1]);
                    eta = pow((Mc/Mtot), (5.0/3.0));
                     if(eta > 0.25)
                     {
                        dm = 0.0;
                     }
                     else
                     {
                        dm = sqrt(1.0-4.0*eta);
                     }
                    m1 = Mtot*(1.0+dm)/2.0;
                    m2 = Mtot*(1.0-dm)/2.0;
                    }
                    
                    printf("%d %e %e %e %e %f %f %f %f %f\n", mc, logLx[q], paramx[q][5], Mc, heat[NC-1],
                    (double)(sacc[0])/(double)(scount[0]),
                    (double)(av[0][q])/(double)(cv[0][q]),
                    (double)(av[1][q])/(double)(cv[1][q]),
                    (double)(av[2][q])/(double)(cv[2][q]),
                    (double)(av[3][q])/(double)(cv[3][q]));
                }
        
     }
    
    fclose(chain);
    fclose(levels);
    
    sprintf(filename, "start_%d_%d.dat", seg, rep);
    chain = fopen(filename,"w");
    
    if(ll == 0)
    {
     m1 = pmax[0];
     m2 = pmax[1];
    }
    
    if(ll==1)
    {
     m1 = exp(pmax[0]);
     m2 = exp(pmax[1]);
     DL = exp(pmax[6]);
    }
    
    if(ll == 2)
    {
    DL = exp(pmax[6]);
    Mc = exp(pmax[0]);
    Mtot = exp(pmax[1]);
    eta = pow((Mc/Mtot), (5.0/3.0));
     if(eta > 0.25)
     {
        dm = 0.0;
     }
     else
     {
        dm = sqrt(1.0-4.0*eta);
     }
    m1 = Mtot*(1.0+dm)/2.0;
    m2 = Mtot*(1.0-dm)/2.0;
    }
    
    
    
    
    fprintf(chain,"%.12e %.12e %.12e ", logLmax, m1, m2);
    for(i = 2; i < NP; i++)
     {
     if(i == 6)
      {
      fprintf(chain, "%.12e ", DL);
      }
      else
      {
       fprintf(chain, "%.15e ", pmax[i]);
      }
    }
    fclose(chain);
    
    // use the last entry to pass the log likelihood over to the sky search
    pmax[NP-1] = logLmax;
    
    free_double_tensor(history,NC,NH);
    free_double_tensor(Fisher,NC,NV);
    free_double_tensor(evec,NC,NV);
    free_double_matrix(eval,NC);
    free_int_vector(mh);
    free_int_vector(sacc);
    free_int_vector(scount);
    free_int_vector(who);
    free_double_vector(max);
    free_double_vector(min);
    free_double_vector(heat);
    free_double_vector(logLx);
    free_int_matrix(av,5);
    free_int_matrix(cv,5);
   
}

void update(int mc, int k, int ll, double *logLx, double **paramx, double **eval, double ***evec, double ***history, double *min, double *max, int *who, double *heat, int **av, int **cv, double *AC, double *EC, double *SA, double *SE, double Tobs, double Tzero, int N, gsl_rng *r)
{
    int q, i, j;
    int typ, flag;
    double lm1, lm2, s1, s2;
    double alpha, beta;
    double logLy, H, x;
    double SNR;
    double *paramy;
    double a, b, c, leta;

    paramy = double_vector(NP);
    
    q = who[k];
    
    for(i = 0; i < NP; i++) paramy[i] = paramx[q][i];
    
    alpha = gsl_rng_uniform(r);
    
    if(mc < MS/8)
    {
        a = 0.5;
        b = 0.3;
        c = 0.0;
    }
    else if (mc < MS/4)
    {
        a = 0.8;
        b = 0.4;
        c = 0.1;
    }
    else if (mc < MS/2)
    {
        a = 0.95;
        b = 0.5;
        c = 0.3;
    }
    else
    {
        a = 1.0;
        b = 0.5;
        c = 0.3;
    }
    
    if(alpha > a ) // uniform draw
    {
        typ = 0;
        
        for (j=0; j< NP; j++) paramy[j] = min[j]+(max[j]-min[j])*gsl_rng_uniform(r);
        
        // keep tc at current value part of the time
        beta = gsl_rng_uniform(r);
        if(beta > 0.3) paramy[5] = paramx[q][5];
        
    }
    else if (alpha > b && logLx[q] > 20.0)// fisher jump
    {
        typ = 1;
        
        beta = gsl_rng_uniform(r);
        i = (int)(beta*(NV));
        beta = sqrt(heat[k])*eval[q][i]*gsl_ran_gaussian(r,1.0);
        for(j = 0; j < NV; j++) paramy[j] += beta*evec[q][i][j];
        
    }
    else if (alpha > c)  // uncorrelated Gaussian jumps
    {
        typ = 2;
        
        SNR = sqrt(2.0*logLx[q]);
        x = 1.0/SNR;
        
        beta = gsl_rng_uniform(r);
        
        if(beta > 0.6)
        {
        paramy[0] += gsl_ran_gaussian(r,1.0e-2*x);
        paramy[1] += gsl_ran_gaussian(r,1.0e-2*x);
        paramy[2] += gsl_ran_gaussian(r,0.1*x);
        paramy[3] += gsl_ran_gaussian(r,0.1*x);
        paramy[5] += gsl_ran_gaussian(r,1.0e2*x);
        }
        else if(beta > 0.3)
        {
            paramy[0] += gsl_ran_gaussian(r,1.0e-1*x);
            paramy[1] += gsl_ran_gaussian(r,1.0e-1*x);
            paramy[2] += gsl_ran_gaussian(r,x);
            paramy[3] += gsl_ran_gaussian(r,x);
            paramy[5] += gsl_ran_gaussian(r,1.0e3*x);
        }
        else
        {
            paramy[0] += gsl_ran_gaussian(r,1.0e-3*x);
            paramy[1] += gsl_ran_gaussian(r,1.0e-3*x);
            paramy[2] += gsl_ran_gaussian(r,0.01*x);
            paramy[3] += gsl_ran_gaussian(r,0.01*x);
            paramy[5] += gsl_ran_gaussian(r,10.0*x);
        }
        
    }
    else  // DE
    {
        typ = 3;
        
        // the history file is kept for each temperature
        
        de_jump(paramx[q], paramy, history[k], NH, NV, r);
    }
    

    
    cv[typ][k]++;
    
    if(ll < 2)
    {
    lm1 = paramy[0];
    lm2 = paramy[1];
    s1 = paramy[2];
    s2 = paramy[3];
    if(lm2 > lm1)
    {
        paramy[0] = lm2;
        paramy[1] = lm1;
        paramy[2] = s2;
        paramy[3] = s1;
    }
    }
    
    if(ll == 2)
    {
     // eta cannot exceed 0.25
     leta = (5.0/3.0)*(paramy[0]-paramy[1]);
     if(leta > log(0.25)) paramy[1] = paramy[0] + 0.856; // a little below equal mass
    }
    
    logLy = -1.0e10;
    flag = 0;
    for(i = 0; i < NV; i++)
    {
        if(paramy[i] > max[i] || paramy[i] < min[i]) paramy[i] = min[i]+(max[i]-min[i])*gsl_rng_uniform(r);
        if(paramy[i] != paramy[i]) flag = 1; // crazier things have happened...
    }
    
    x = (paramy[5]-Tzero)/Tobs;
    if( fabs(x-rint(x)) < buffer ) flag = 1;  // keep merger away from boundaries

    
    if(flag == 0)
     {
     if (mc < MS/8)  // time, phase and amplitude maximized
      {
      logLy = log_likelihood_max_dual(ll, AC, EC, paramy, SA, SE, N, Tobs, Tzero);
      }
      else // phase and amplitude maximized
      {
       logLy = log_likelihood_APmax_dual(ll, AC, EC, paramy, SA, SE, N, Tobs, Tzero);
      }
     }
    
    H = (logLy-logLx[q])/heat[k];
     
    alpha = log(gsl_rng_uniform(r));
        
     if(H > alpha)
     {
         // copy over new state if accepted
         logLx[q] = logLy;
         for(i = 0; i < NP; i++) paramx[q][i] = paramy[i];
         av[typ][k]++;
     }
    
    free_double_vector(paramy);
    
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
    
    for(k=0; k< d; k++) paramsy[k] = paramsx[k]+alpha*(history[i][k]-history[j][k]);
    
}


double findsky(int ll, double *params, double *min, double *max, double Tobs, double Tzero, int N,  double *AC, double *EC, double *SA, double *SE, gsl_rng *r)
{
    int i, j, k;
   double logL, logLmax, x, ts, tref, alpha;
   double *paramy, *pmax;
    
    pmax = double_vector(NP);
    paramy = double_vector(NP);
    
    tref = params[5];
    
    // copy over parameters
    for (j=0; j< NP; j++) paramy[j] = params[j];
    
    for (j=7; j< 9; j++) paramy[j] = min[j]+(max[j]-min[j])*gsl_rng_uniform(r);
    for (j=0; j< NP; j++) pmax[j] = paramy[j];
    
    logLmax = -1.0;
    
    k = 0;
    do
    {
        
     paramy[6] = 0.0;
     for (j=7; j< 9; j++) paramy[j] = min[j]+(max[j]-min[j])*gsl_rng_uniform(r);
  
     ts = Tmap(paramy, params[5]);
     paramy[5] = params[5]-ts;
    
     x = (paramy[5]-Tzero)/Tobs;
     logL = -100.0;
        
    if(fabs(x-rint(x)) > buffer)
    {
     logL = likelihoodFstatTmax(ll, paramy, N, AC, EC, SA, SE, Tobs, Tzero);
    //logL = likelihoodFstat(ll, paramy, N, AC, EC, SA, SE, Tobs, Tzero);
    }
     
    if(logL > logLmax)
    {
        logLmax = logL;
        for (j=0; j< NP; j++) pmax[j] = paramy[j];
    }
        
   // printf("%d %e %e\n", k, logL, logLmax);
        
    k++;
        
    }while(k < 10);
    
    for (j=0; j< NP; j++) params[j] = pmax[j];
    
    free_double_vector(paramy);
    free_double_vector(pmax);
    
    return logLmax;
    
}

// maps tc so that merger time at detector is held fixed
double Tmap(double *params, double tdet)
{
    
    /*   Indicies   */
    int i,j, k, n, m, a, M;
    
    /*   Gravitational Wave basis vectors   */
    double *kv;
    
    /*   Spacecraft position and separation vector   */
    double *x, *y, *z;

    /*   Dot products   */
    double kdotx;
    
    /*   GW Source data   */
    double Mc, theta, phi, psi, D, iota, A, Aplus, Across, f0, fdot, phio;
    double costh, sinth, cosph, sinph, cosi, cosps, sinps;
    
    /*   Time and distance variables   */
    double xa, ya, za;
    
    /*   Miscellaneous  */
    double xm, td;
    
    /*   Allocating Arrays   */
    
    kv = double_vector(4);
    x = double_vector(4); y = double_vector(4); z = double_vector(4);

    phi = params[8];   // EclipticLongitude
    //Calculate cos and sin of sky position
    costh = params[7];   sinth = sqrt(1.0-costh*costh);
    cosph = cos(phi);     sinph = sin(phi);

    kv[1] = -sinth*cosph;
    kv[2] = -sinth*sinph;
    kv[3] = -costh;
    
    spacecraft(tdet, x, y, z);
        
        // guiding center
        xa = (x[1]+x[2]+x[3])/3.0;
        ya = (y[1]+y[2]+y[3])/3.0;
        za = (z[1]+z[2]+z[3])/3.0;
        
    kdotx = (xa*kv[1]+ya*kv[2]+za*kv[3])/clight;
    
    free_double_vector(kv);
    free_double_vector(x); free_double_vector(y); free_double_vector(z);
    
    return kdotx;
}

void searchsky(double *pmax, double **paramx, double *AC, double *EC, double *SAS, double *SES, double *SA, double *SE, double Tobs, int seg, int N, int rep)
{
    int i, j, k, mc, hold;
    double x, logLmax, ts;
    double Mc, Mtot, eta, dm;
    double *logLx;
    double *min, *max;
    int *scount, *sacc;
    double ***Fisher, **eval, ***evec;
    int **av, **cv;
    int ll=mtype;
    int *who;
    double *heat;
    int q;
    double m1, m2, DL, Tzero;
    double alpha, beta;
    char filename[1024];
    
    FILE *chain;
    FILE *levels;
    
    const gsl_rng_type * P;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    
    P = gsl_rng_default;
    r = gsl_rng_alloc (P);
    
    sacc = int_vector(NC);
    scount = int_vector(NC);
    who = int_vector(NC);
    max = double_vector(NP);
    min = double_vector(NP);
    heat = double_vector(NC);
    logLx = double_vector(NC);
    av = int_matrix(4,NC);
    cv = int_matrix(4,NC);
    
    Fisher = double_tensor(NC,NP,NP);
    eval = double_matrix(NC,NP);
    evec = double_tensor(NC,NP,NP);
    
    Tzero = (double)(seg)*Tobs/2.0;
    
    for (i=0; i< NC; i++) who[i] = i;
    // effective SNR of 7 for hottest chain
    logLmax = pmax[NP-1];
    alpha = pow((2.0*logLmax/50.0),1.0/(double)(NC));
    if(alpha > 1.3) alpha = 1.3;
    heat[0] = 1.0;
    for (i=1; i< NC; i++) heat[i] = heat[i-1]*alpha;
  
    
    for(j = 0; j < 4; j++)
    {
        for(k=0; k < NC; k++)
        {
            av[j][k] = 0;
            cv[j][k] = 1;
        }
    }
    
    for(k=0; k < NC; k++)
    {
        sacc[k] = 0;
        scount[k] = 1;
    }
   
    max[0] = log(1.0e9);
    max[1] = log(2.3e9);  // (should be not more than 2.3 times Mcmax)
    max[2] = 0.95;
    max[3] = 0.95;
    max[4] = PI;
    max[5] = Tzero+2.0*Tobs-buffer*Tobs;
    max[6] = log(1.0e3);
    max[7] = 1.0;
    max[8] = 2.0*PI;
    max[9] = PI;
    max[10] = 1.0;
   
    min[0] = log(1.0e4);
    min[1] = log(1.0e4);
    min[2] = -0.95;
    min[3] = -0.95;
    min[4] = 0.0;
    min[5] = Tzero+buffer*Tobs;
    min[6] = log(0.1);
    min[7] = -1.0;
    min[8] = 0.0;
    min[9] = 0.0;
    min[10] = -1.0;
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
    // [7] cos(EclipticCoLatitude), [8] EclipticLongitude  [9] polarization, [10] inclination
    
     #pragma omp parallel for
     for (i=0; i< NC; i++)
       {
           logLx[i] = findsky(ll, paramx[i], min, max, Tobs, Tzero, N, AC, EC, SA, SE, rvec[i]);
       }
    
    sprintf(filename, "searchsky_%d_%d.dat", seg, rep);
    chain = fopen(filename,"w");
    sprintf(filename, "likesky_%d_%d.dat", seg, rep);
    levels = fopen(filename,"w");
    
    logLmax = -1.0e10;
    
    for(mc = 0; mc < MSS; mc++)
    {
        
       // printf("mc = %d \n", mc);
    
     alpha = gsl_rng_uniform(r);
        
        if(mc%50==0)
        {
         #pragma omp parallel for
         for(i = 0; i < NC; i++)
         {
            FisherFast(ll, paramx[i], Fisher[i], SAS, SES, Tobs, Tzero, N);
            FisherEvec(Fisher[i], eval[i], evec[i], NP);
         }
        }
    
      if((NC > 1) && (alpha < 0.3))  // decide if we are doing a MCMC update of all the chains or a PT swap
       {
        // chain swap
       
        for(i = 0; i < NC/2; i++)
        {
        alpha = (double)(NC-1)*gsl_rng_uniform(r);
        j = (int)(alpha);
        beta = ((logLx[who[j]]-logLx[who[j+1]])/heat[j+1] - (logLx[who[j]]-logLx[who[j+1]])/heat[j]);
        alpha = log(gsl_rng_uniform(r));
        if(beta > alpha)
        {
            hold = who[j];
            who[j] = who[j+1];
            who[j+1] = hold;
            sacc[j]++;
        }
        
         scount[j]++;
        }
        
       }
      else      // MCMC update
       {
           #pragma omp parallel for
           for(k=0; k < NC; k++)
            {
             updatesky(mc, k, ll, logLx, paramx, eval, evec, min, max, who, heat, av, cv, AC, EC, SA, SE, Tobs, Tzero, N, rvec[k]);
            }
       }
        
        for(k=0; k < NC; k++)
        {
          if(logLx[k] > logLmax)
          {
              logLmax = logLx[k];
              for(i = 0; i < NP; i++) pmax[i] = paramx[k][i];
          }
        }
        
        /*
        if(mc%100 == 0)
        {
            
            // clone the best solution into the hottest chain
            x = -1.0e4;
            for(k=0; k < NC; k++)
            {
              if(logLx[k] > x)
              {
                  x = logLx[k];
                  j = k;
              }
            }
            
            q = who[NC-1];
            logLx[q] = logLx[j];
            for(i = 0; i < NP; i++) paramx[q][i] = paramx[j][i];
            
        } */
        
        if(mc%10 == 0)
        {
            
                q = who[0];
            
                      if(ll == 0)
                       {
                        m1 = paramx[q][0];
                        m2 = paramx[q][1];
                        DL = paramx[q][6];
                       }
                       if(ll==1)
                       {
                        m1 = exp(paramx[q][0]);
                        m2 = exp(paramx[q][1]);
                        DL = exp(paramx[q][6]);
                       }
                       if(ll == 2)
                       {
                       DL = exp(paramx[q][6]);
                       Mc = exp(paramx[q][0]);
                       Mtot = exp(paramx[q][1]);
                       eta = pow((Mc/Mtot), (5.0/3.0));
                        if(eta > 0.25)
                        {
                           dm = 0.0;
                        }
                        else
                        {
                           dm = sqrt(1.0-4.0*eta);
                        }
                       m1 = Mtot*(1.0+dm)/2.0;
                       m2 = Mtot*(1.0-dm)/2.0;
                       }
            
               fprintf(chain,"%d %.12e %.12e %.12e ", mc, logLx[q], m1, m2);
               for(i = 2; i < NP; i++)
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
              fprintf(chain,"%d\n", q);
            
              fprintf(levels, "%d ", mc);
               for(k = 0; k < NC; k++)
                {
                fprintf(levels, "%.12e ", logLx[who[k]]);
                }
                 fprintf(levels, "\n");
          }
        
              if(mc%10 == 0)
                {
                    q = who[0];
                    //x =  Likelihood(ll, paramx[q], N, AC, EC, SA, SE, Tobs, Tzero);
                    printf("%d %e %e %e %e %f %f %f %f\n", mc, logLx[q], paramx[q][5], paramx[q][7], paramx[q][8],
                    (double)(sacc[0])/(double)(scount[0]),
                    (double)(av[0][q])/(double)(cv[0][q]),
                    (double)(av[1][q])/(double)(cv[1][q]),
                    (double)(av[2][q])/(double)(cv[2][q]));
                }
        
     }
    
    fclose(chain);
    fclose(levels);
    
    sprintf(filename, "skymax_%d_%d.dat", seg, rep);
    chain = fopen(filename,"w");
              if(ll == 0)
              {
                m1 = pmax[0];
                m2 = pmax[1];
                DL = pmax[6];
               }
               if(ll==1)
               {
                m1 = exp(pmax[0]);
                m2 = exp(pmax[1]);
                DL = exp(pmax[6]);
               }
               if(ll == 2)
               {
               DL = exp(pmax[6]);
               Mc = exp(pmax[0]);
               Mtot = exp(pmax[1]);
               eta = pow((Mc/Mtot), (5.0/3.0));
                if(eta > 0.25)
                {
                   dm = 0.0;
                }
                else
                {
                   dm = sqrt(1.0-4.0*eta);
                }
               m1 = Mtot*(1.0+dm)/2.0;
               m2 = Mtot*(1.0-dm)/2.0;
               }
    
    fprintf(chain,"%.12e %.15e %.15e ", logLmax, m1, m2);
    for(i = 2; i < NP; i++)
     {
     if(i == 6)
      {
      fprintf(chain, "%.15e ", DL);
      }
      else
      {
       fprintf(chain, "%.15e ", pmax[i]);
      }
    }
    fclose(chain);


    free_double_tensor(Fisher,NC,NP);
    free_double_tensor(evec,NC,NP);
    free_double_matrix(eval,NC);
    free_int_vector(sacc);
    free_int_vector(scount);
    free_int_vector(who);
    free_double_vector(max);
    free_double_vector(min);
    free_double_vector(heat);
    free_double_vector(logLx);
    free_int_matrix(av,4);
    free_int_matrix(cv,4);
   
}

void updatesky(int mc, int k, int ll, double *logLx, double **paramx, double **eval, double ***evec, double *min, double *max, int *who, double *heat, int **av, int **cv, double *AC, double *EC, double *SA, double *SE, double Tobs, double Tzero, int N, gsl_rng *r)
{
    int q, i, j;
    int typ, flag;
    double lm1, lm2;
    double alpha, beta;
    double tx, ty;
    double logLy, H, x;
    double s1, s2, leta;
    double *paramy;
    double a, b;
    double pDx, pDy;
    double **Chl, **Cov;
    double *zv;
    
    paramy = double_vector(NP);
   
    q = who[k];
    
    for(i = 0; i < NP; i++) paramy[i] = paramx[q][i];
    
    alpha = gsl_rng_uniform(r);
    
    if(mc < MSS/4)
    {
        a = 0.6;
        b = 0.2;
    }
    else if (mc < MSS/2)
    {
        a = 0.8;
        b = 0.4;
    }
    else
    {
        a = 0.9;
        b = 0.4;
    }
    
    
    if(alpha > a ) // uniform draw on sky
    {
        typ = 0;
        
        for (j=7; j< 9; j++) paramy[j] = min[j]+(max[j]-min[j])*gsl_rng_uniform(r);
    
        // time delay for current sky location
        tx = Tmap(paramx[q], paramx[q][5]);
        // time delay for new current sky location
        ty = Tmap(paramy, paramy[5]);
        
        paramy[5] = paramx[q][5]+tx-ty;
        
        
    }
    else if (alpha > b) //  sky jump
    {
        typ = 1;
     
       beta = gsl_rng_uniform(r);
        
        if(beta > 0.7)
        {
           paramy[5] += gsl_ran_gaussian(r,1.0);
           paramy[7] += gsl_ran_gaussian(r,0.001);
           paramy[8] += gsl_ran_gaussian(r,0.002);
        }
        else if(beta > 0.3)
        {
           paramy[5] += gsl_ran_gaussian(r,10.0);
           paramy[7] += gsl_ran_gaussian(r,0.01);
           paramy[8] += gsl_ran_gaussian(r,0.02);
        }
        else
        {
            paramy[7] += gsl_ran_gaussian(r,0.05);
            paramy[8] += gsl_ran_gaussian(r,0.1);
            // time delay for current sky location
            tx = Tmap(paramx[q], paramx[q][5]);
            // time delay for new current sky location
            ty = Tmap(paramy, paramy[5]);
            paramy[5] = paramx[q][5]+tx-ty;
        }
        
        
    }
    else  // Fisher update for all parameters
    {
        typ = 2;
        
        // pick an eigendirection to jump in
        beta = gsl_rng_uniform(r);
        i = (int)(beta*(double)(NP));
        // draw the jump size
        beta = sqrt(heat[k])*eval[q][i]*gsl_ran_gaussian(r,1.0);
        for(j = 0; j < NP; j++) paramy[j] = paramx[q][j]+beta*evec[q][i][j];
        
        if(ll < 2)
        {
        lm1 = paramy[0];
        lm2 = paramy[1];
        s1 = paramy[2];
        s2 = paramy[3];
        if(lm2 > lm1)
        {
            paramy[0] = lm2;
            paramy[1] = lm1;
            paramy[2] = s2;
            paramy[3] = s1;
        }
        }
        
        if(ll == 2)
        {
         // eta cannot exceed 0.25
         leta = (5.0/3.0)*(paramy[0]-paramy[1]);
         if(leta > log(0.25)) paramy[1] = paramy[0] + 0.856; // a little below equal mass
        }
        
    }
    
    if(paramy[7] > 1.0)  paramy[7] = 0.99;
    if(paramy[7] < -1.0)  paramy[7] = -0.99;
    if(paramy[8] > TPI)  paramy[8] -= TPI;
    if(paramy[8] < 0.0)  paramy[8] += TPI;
    
    for(i = 0; i < NP; i++)
    {
        if(paramy[i] > max[i] || paramy[i] < min[i])
        {
            paramy[i] = min[i]+(max[i]-min[i])*gsl_rng_uniform(r);
        }
    }
    
    cv[typ][k]++;
    
    
       logLy = -1.0e10;
       flag = 0;
       for(i = 0; i < NV; i++)
       {
           if(paramy[i] != paramy[i]) flag = 1; // crazier things have happened...
       }
    
       x = (paramy[5]-Tzero)/Tobs;
       if( fabs(x-rint(x)) < buffer ) flag = 1;  // keep merger away from boundaries
    
      if(flag == 0)
      {
         logLy = likelihoodFstat(ll, paramy, N, AC, EC, SA, SE, Tobs, Tzero);
      }
    
 
    // catch failures in the maximization
    if(paramy[4] != paramy[4]) paramy[4] = PI;
    if(paramy[6] != paramy[6]) paramy[6] = 1.0;
    if(paramy[9] != paramy[9]) paramy[9] = PI/4.0;
    if(paramy[10] != paramy[10]) paramy[10] = 0.0;
    

    pDx = paramx[q][6];   // uniform in distance prior
    pDy = paramy[6];   // uniform in distance prior
    
    H = (logLy-logLx[q])/heat[k] + pDy - pDx;
     
    alpha = log(gsl_rng_uniform(r));
        
     if(H > alpha)
     {
         // copy over new state if accepted
         logLx[q] = logLy;
         for(i = 0; i < NP; i++) paramx[q][i] = paramy[i];
         av[typ][k]++;
     }
    
    free_double_vector(paramy);
 
    
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

double fourier_nwip_Tshift(double delT, double Tobs, double *a, double *b, double *Sn, int n)
{
    int i, j, k;
    double arg, product;
    double test;
    double ReA, ReB, ImA, ImB;
    double f, c, s;
    
    arg = 0.0;
    for(i=1; i<n/2; i++)
    {
        j = i;
        k = n-i;
        
        f = (double)(i)/Tobs;
        c = cos(TPI*f*delT);
        s = sin(TPI*f*delT);
        
        ReA = a[j]*c+a[k]*s;
        ImA = -a[j]*s+a[k]*c;
        
        ReB = b[j]; ImB = b[k];
        product = ReA*ReB + ImA*ImB;
        arg += product/(Sn[i]);
    }
    
    return(4.0*arg);
    
}



double f_at_t(double m1, double m2, double chi1, double chi2, double tc, double t)
{
    // 3PN f(t)
    int i;
    double f, fr, fny, af, M, eta, chi, theta;
    double gamma_E=0.5772156649; //Euler's Constant-- shows up in 3PN term
    double PN1, PN15, PN2, PN25, PN3, PN35;
    double theta2, theta3, theta4, theta5, theta6, theta7;
    
    M = m1+m2;
    eta = m1*m2/(M*M);
    chi = (m1*chi1+m2*chi2)/M;
    
    af = FinalSpin0815(eta, chi1, chi2);
    fr = fring(eta, chi1, chi2, af)/M;
    
    // Taylor T3
    
    if(tc > t)
    {
        
        theta = pow(eta*(tc-t)/(5.0*M),-1.0/8.0);
        theta2 = theta*theta;
        theta3 = theta2*theta;
        theta4 = theta2*theta2;
        theta5 = theta2*theta3;
        theta6 = theta3*theta3;
        theta7 = theta3*theta4;
        
        PN1 = (11.0/32.0*eta+743.0/2688.0)*theta2;
        PN15 = -3.0*PI/10.0*theta3 + (1.0/160.0)*(113.0*chi-38.0*eta*(chi1+chi2))*theta3;
        PN2 = (1855099.0/14450688.0+56975.0/258048.0*eta+371.0/2048.0*eta*eta)*theta4 + (1.0/14450688.0)*(-3386880.0*chi*chi+1512.0*chi1*chi2)*theta4;
        PN25 = -(7729.0/21504.0-13.0/256.0*eta)*PI*theta5;
        PN3 = (-720817631400877.0/288412611379200.0+53.0/200.0*PI*PI+107.0/280.0*gamma_E
               +(25302017977.0/4161798144.0-451.0/2048.0*PI*PI)*eta-30913.0/1835008.0*eta*eta+235925.0/1769472.0*eta*eta*eta +107.0/280.0*log(2.0*theta))*theta6;
        PN35 = (141769.0/1290240.0*eta*eta - 97765.0/258048.0*eta - 188516689.0/433520640.0)*PI*theta7;
        
        //printf("%f %f %e %e %e %e %e\n", t, theta3/(8.0*M*PI), PN1, PN15, PN2, PN25, PN3);
        
        f = theta3/(8.0*M*PI)*(1.0 + PN1 + PN15 + PN2 + PN25 + PN3 + PN35);
        
        if(f > fr) f = 2.0*fr;
    }
    else
    {
        f = 2.0*fr;
    }
    
    return(f);
    
}




//uses frequency domain data
void qscanf(char filename[], double *data, double *Sn, double Tobs, int N)
{
    double t_rise, alpha;
    double dt, fmx, fmn;
    double *dcopy;
    double *freqs;
    double x, dx, dlnf;
    double t, f;
    int subscale, octaves, Nf, i, j;
    double **tfDR, **tfDI, **tfD;
    double fac;
    double *SX;
    
    FILE *out;
    
    dt = Tobs/(double)(N);
    fmx = 1.0/(2.0*dt);
    fmn = 8.0/Tobs;
    
    dcopy = (double*)malloc(sizeof(double)* (N));
    SX = (double*)malloc(sizeof(double)* (N/2));
    
    // Qscan uses different scaling convention
   // fac = Tobs/((double)(N)*(double)(N));
    for (i = 0; i < N/2; ++i) SX[i] = Sn[i];

    // make copy
    for (i = 0; i < N; ++i) dcopy[i] = data[i];

    whiten(dcopy, SX, N);
    
    // logarithmic frequency spacing
    subscale = 20;  // number of semi-tones per octave
    octaves = (int)(rint(log(fmx/fmn)/log(2.0))); // number of octaves
    Nf = subscale*octaves+1;
    freqs = (double*)malloc(sizeof(double)* (Nf));   // frequencies used in the analysis
    dx = log(2.0)/(double)(subscale);
    dlnf = dx;
    x = log(fmn);
    for(i=0; i< Nf; i++)
    {
        freqs[i] = exp(x);
        x += dx;
    }
    
    printf("%d freqs %d times\n", Nf, N);
    
    tfDR = double_matrix(Nf,N);
    tfDI = double_matrix(Nf,N);
    tfD = double_matrix(Nf,N);
    
    // Wavelet transform
    TransformC(dcopy, freqs, tfD, tfDR, tfDI, 8.0, Tobs, N, Nf);
    
    
    out = fopen(filename,"w");
    for(j = 0; j < Nf; j+=2)
    {
        f = freqs[j];
        
        for(i = 0; i < N; i+=16)
        {
            t = (double)(i)*dt;
            fprintf(out,"%e %e %e\n", t, f, tfD[j][i]);
        }
        
        fprintf(out,"\n");
    }
    fclose(out);
    
    free(SX);
    free(dcopy);
    free(freqs);
    free_double_matrix(tfDR,Nf);
    free_double_matrix(tfDI,Nf);
    free_double_matrix(tfD,Nf);
    
}

void whiten(double *data, double *Sn, int N)
{
    double f, x, y, fix;
    int i;
    
    data[0] = 0.0;
    data[N/2] = 0.0;
    
    for(i=1; i< N/2; i++)
    {
        x = 1.0/sqrt(Sn[i]);
        data[i] *= x;
        data[N-i] *= x;
    }
    
}

void TransformC(double *a, double *freqs, double **tf, double **tfR, double **tfI, double Q, double Tobs, int n, int m)
{
    int j;
    double fix;
    
    // [0] t0 [1] f0 [2] Q [3] Amp [4] phi
    
    fix = sqrt((double)(n/2));
    
    #pragma omp parallel for
    for(j = 0; j < m; j++)
    {
        layerC(a, freqs[j], tf[j], tfR[j], tfI[j], Q, Tobs, fix, n);
    }
    
}

void layerC(double *a, double f, double *tf, double *tfR, double *tfI, double Q, double Tobs, double fix, int n)
{
    int i;
    double *AC, *AF;
    double *b;
    double *params;
    double bmag;
    
    params= double_vector(6);
    AC=double_vector(n);  AF=double_vector(n);
    b = double_vector(n);
    
    params[0] = 0.0;
    params[1] = f;
    params[2] = Q;
    params[3] = 1.0;
    params[4] = 0.0;
    
    SineGaussianC(b, params, Tobs, n);
    
    bmag = sqrt(f_nwip(b, b, n)/(double)n);
    
    bmag /= fix;
    
    phase_blind_time_shift(AC, AF, a, b, n);
    
    for(i = 0; i < n; i++)
    {
        tfR[i] = AC[i]/bmag;
        tfI[i] = AF[i]/bmag;
        tf[i] = tfR[i]*tfR[i]+tfI[i]*tfI[i];
    }
    
    free_double_vector(AC);  free_double_vector(AF);
    free_double_vector(b);  free_double_vector(params);
    
}

void SineGaussianC(double *hs, double *sigpar, double Tobs, int N)
{
    double f0, t0, Q, sf, sx, Amp;
    double fmx, fmn;//, fac;
    double phi, f;
    double tau;
    double re,im;
    double q, p, u;
    double A, B, C;
    double Am, Bm, Cm;
    double a, b, c;
    double dt, fac;
    double cosPhase_m, sinPhase_m, cosPhase_p, sinPhase_p;
    
    int i, imid, istart,istop,imin,imax,iend,even,odd;
    
    t0  = sigpar[0];
    f0  = sigpar[1];
    Q   = sigpar[2];
    Amp = sigpar[3];
    phi = sigpar[4];
    
    tau = Q/(TPI*f0);
    
    fmx = f0 + 3.0/tau;  // no point evaluating waveform past this time (many efolds down)
    fmn = f0 - 3.0/tau;  // no point evaluating waveform before this time (many efolds down)
    
    i = (int)(f0*Tobs);
    imin = (int)(fmn*Tobs);
    imax = (int)(fmx*Tobs);
    if(imax - imin < 10)
    {
        imin = i-5;
        imax = i+5;
    }
    
    if(imin < 0) imin = 1;
    if(imax > N/2) imax = N/2;
    
    hs[0] = 0.0;
    hs[N/2] = 0.0;
    
    for(i = 1; i < N/2; i++)
    {
        hs[i] = 0.0;
        hs[N-i] = 0.0;
    }
    
    dt = Tobs/(double)(N);
    fac = sqrt(sqrt(2.0)*PI*tau/dt);
    
    /* Use recursion relationship  */
    
    imid = (int)(f0*Tobs);
    
    p = PI*PI*tau*tau/(Tobs*Tobs);
    
    Bm = exp(-p*(((double)(imid)-f0*Tobs)*((double)(imid)-f0*Tobs)));
    Cm = 1.0;
    
    b = exp(-p*(1.0+2.0*((double)(imid)-f0*Tobs)));
    c = exp(-2.0*p);
    
    // start in the middle and work outwards
    
    B = Bm;
    C = Cm;

    
    for(i = imid; i < imax; i++)
    {
 
        f = (double)(i)/Tobs;
        
        sf = fac*B;
        
        //  printf("%e\n", exp(-pi2tau2*(f-f0)*(f-f0))-B);
        
        hs[i] = sf;
        hs[N-i] = 0.0;
        
        B *= (C*b);
        C *= c;
        
    }
    
    // reset to midpoint
    
    B = Bm;
    C = Cm;
    
    b = exp(p*(-1.0+2.0*((double)(imid)-f0*Tobs)));
    // c unchanged
    
    for(i = imid; i > imin; i--)
    {

        f = (double)(i)/Tobs;
        
        sf = fac*B;
        
        // printf("%e\n", exp(-pi2tau2*(f-f0)*(f-f0))-B);
        
        hs[i] = sf;
        hs[N-i] = 0.0;
        
        B *= (C*b);
        C *= c;
        

    }
    
    
}

double f_nwip(double *a, double *b, int n)
{
    int i, j, k;
    double arg, product;
    double test;
    double ReA, ReB, ImA, ImB;
    
    arg = 0.0;
    for(i=1; i<n/2; i++)
    {
        j = i;
        k = n-j;
        ReA = a[j]; ImA = a[k];
        ReB = b[j]; ImB = b[k];
        product = ReA*ReB + ImA*ImB;
        arg += product;
    }
    
    return(arg);
    
}

void phase_blind_time_shift(double *corr, double *corrf, double *data1, double *data2, int n)
{
    int nb2, i, l, k, j;
    int imax, imin;
    
    nb2 = n / 2;
    
    corr[0] = 0.0;
    corrf[0] = 0.0;
    corr[nb2] = 0.0;
    corrf[nb2] = 0.0;
    
    for (i=1; i < nb2; i++)
    {
        l=i;
        k=n-i;
        
        corr[l]    = (data1[l]*data2[l] + data1[k]*data2[k]);
        corr[k]    = (data1[k]*data2[l] - data1[l]*data2[k]);
        corrf[l] = corr[k];
        corrf[k] = -corr[l];
    }
    
    gsl_fft_halfcomplex_radix2_inverse(corr, 1, n);
    gsl_fft_halfcomplex_radix2_inverse(corrf, 1, n);
    
    
}





void instrument_noise(double f, double *SAE, double *SXYZ)
{
    //Power spectral density of the detector noise and transfer frequency
    double Sn, red, confusion_noise;
    double Sloc, fonfs;
    double f1, f2;
    double A1, A2, slope, LC;
    FILE *outfile;
    
    fonfs = f/fstar;
    
    LC = 4.0*fonfs*fonfs;
    
    red = 16.0*((1.0e-4/f)*(1.0e-4/f));
    
    // Calculate the power spectral density of the detector noise at the given frequency
    
    *SAE = LC*16.0/3.0*pow(sin(fonfs),2.0)*( (2.0+cos(fonfs))*(Sps) + 2.0*(3.0+2.0*cos(fonfs)+cos(2.0*fonfs))*(Sacc/pow(2.0*PI*f,4.0)*(1.0+red)) ) / pow(2.0*L,2.0);
    
    *SXYZ = LC*4.0*pow(sin(fonfs),2.0)*( 4.0*(Sps) + 8.0*(1.0+pow(cos(fonfs),2.0))*(Sacc/pow(2.0*PI*f,4.0)*(1.0+red)) ) / pow(2.0*L,2.0);
    
}


void PDwave(int ll, double *wavef, double *params, int N, double Tobs, double Tzero)
{
    
    AmpPhaseFDWaveform *ap = NULL;
    
    double m1, m2, chi1, chi2, Mtot, Mc, eta, dm;
    double m1_SI, m2_SI, distance, phic;
    double phi0, q, f_min, f_max;
    int ret, flag, i, k;
    double p, cp, sp, Amp;
    double f, x, y, deltaF, ts, fr;
    double tc, fstart, fstop, sqrtTobs;
    int imin, imax;
    
    sqrtTobs = sqrt(Tobs);
    
    RealVector *freq;
    
    chi1 = params[2];
    chi2 = params[3];
    phi0 = params[4];
    tc = params[5];
    ts = Tzero+Tobs-tc;
    
    if(ll == 0)  // linear in m1, m2, DL
    {
    m1 = params[0]*TSUN;    // Mass1
    m2 = params[1]*TSUN;    // Mass2
    distance = params[6]*1.0e9*PC_SI; // distance
    Mtot = (m1+m2);
    eta = m1*m2/((m1+m2)*(m1+m2));
    Mc = pow(m1*m2,3.0/5.0)/pow(m1+m2,1.0/5.0);
    }
    else if(ll == 1)  // log in m1, m2, DL
    {
     m1 = exp(params[0])*TSUN;    // Mass1
     m2 = exp(params[1])*TSUN;    // Mass2
     distance = exp(params[6])*1.0e9*PC_SI; // distance
     Mtot = (m1+m2);
     eta = m1*m2/((m1+m2)*(m1+m2));
     Mc = pow(m1*m2,3.0/5.0)/pow(m1+m2,1.0/5.0);
    }
    else // log in Mc, Mt, DL
    {
    distance = exp(params[6])*1.0e9*PC_SI; // distance
    Mc = exp(params[0])*TSUN;
    Mtot = exp(params[1])*TSUN;
    eta = pow((Mc/Mtot), (5.0/3.0));
     if(eta > 0.25)
     {
        dm = 0.0;
     }
     else
     {
        dm = sqrt(1.0-4.0*eta);
     }
    m1 = Mtot*(1.0+dm)/2.0;
    m2 = Mtot*(1.0-dm)/2.0;
    }
    
    m1_SI = m1*MSUN_SI/TSUN;
    m2_SI = m2*MSUN_SI/TSUN;
    
    StartStop(ll, params, Tobs, Tzero, Tzero+Tobs, &fstart, &fstop, &fr);
    
    imin = (int)(fstart*Tobs);
    imax = (int)(fstop*Tobs)+1;
    
    if(imin < 1) imin = 1;
    if(imax <= imin) imax = imin+2;
    if(imax > N/2)
    {
        imax = N/2;
        if(imin >= imax) imin = imax-2;
    }
    
    //printf("%e %e %d %d\n", fstart, fstop, imin, imax);
    
    freq = CreateRealVector((imax-imin));
    for (i=0; i< (imax-imin); i++) freq->data[i] = (double)(i+imin)/Tobs;
    
    if(distance != distance) printf("nan %e\n", params[6]);
    if(distance <= 0.0) printf("<= 0 %e\n", params[6]);

    
    ret = IMRPhenomDGenerateh22FDAmpPhase(&ap, freq, 0.0, PDfref, m1_SI, m2_SI, chi1, chi2, distance);
    
    for (i=0; i< N; i++) wavef[i] = 0.0;
    
    for (i=imin; i< imax; i++)
    {
        k = i-imin;
        f = freq->data[k];
        x = f/fstar;
        Amp = 4.0*x*sin(x)*h22fac*ap->amp[k]/sqrtTobs;
        p = TPI*f*ts-ap->phase[k]-2.0*phi0;
        wavef[i] = Amp*cos(p);
        wavef[N-i] = Amp*sin(p);
    }
   
    
    DestroyAmpPhaseFDWaveform(ap);
    DestroyRealVector(freq);
    
}

void Extrinsic(double *params, double Tend, int NF, double *FF, double *TF, double *PF, double *AF, double *AAmp, double *EAmp, double *APhase, double *EPhase, double *kxm)
{
    
    /*   Indicies   */
    int i,j, k, n, m;
    
    /*   Time and distance variables   */
    double *xi;
    
    double pa, pe;
    double paold, peold;
    double ja, je;
    
    /*   Miscellaneous  */
    double xm, fstep, power, om, mx, tx, x;
    
    double *ta, *xia;
    double *FpAR, *FcAR, *FpER, *FcER;
    double *FpAI, *FcAI, *FpEI, *FcEI;
    
    double t, f, kdotx, Amp, Phase, RR, II, Tstart;
    
    int NA, flag;
    
    double iota, cosi;
    double Aplus, Across;
    
    double Tcut;
    
    clock_t start, end;
    double cpu_time_used;
    
    FILE *out;
    
    
    xi = (double*)malloc(sizeof(double)* (NF));
    FpAR = (double*)malloc(sizeof(double)* (NF));
    FcAR = (double*)malloc(sizeof(double)* (NF));
    FpER = (double*)malloc(sizeof(double)* (NF));
    FcER = (double*)malloc(sizeof(double)* (NF));
    FpAI = (double*)malloc(sizeof(double)* (NF));
    FcAI = (double*)malloc(sizeof(double)* (NF));
    FpEI = (double*)malloc(sizeof(double)* (NF));
    FcEI = (double*)malloc(sizeof(double)* (NF));
    
    Tstart = TF[0];

    RAantenna(params, NF, TF, FF, xi, FpAR, FpAI, FcAR, FcAI, FpER, FpEI, FcER, FcEI);
    
    for(n=0; n< NF; n++)
    {
        AAmp[n] = 0.0;
        EAmp[n] = 0.0;
        APhase[n] = 0.0;
        EPhase[n] = 0.0;
    }
    
    cosi = params[10];  // cos of inclination
    
    Aplus = 0.5*(1.+cosi*cosi);
    Across = -cosi;
    
    // Merger kdotx
    *kxm = (TF[NF-1]-xi[NF-1]);
    
    n = 0;
    RR = FpAR[n]*Aplus - FcAI[n]*Across;
    II = FcAR[n]*Across + FpAI[n]*Aplus;
    paold = atan2(II,RR);
    if(paold < 0.0) paold += TPI;
    RR = FpER[n]*Aplus - FcEI[n]*Across;
    II = FcER[n]*Across + FpEI[n]*Aplus;
    peold = atan2(II,RR);
    if(peold < 0.0) peold += TPI;
    
    ja = 0.0;
    je = 0.0;

    //out = fopen("map_fast.dat","w");
    for(n=0; n< NF; n++)
    {
        // Barycenter time and frequency
        t = TF[n];
        f = FF[n];
        
        // Tukey filter to match what is done in time domain
        x = 1.0;
        if(t < Tstart+t_tuke) x = 0.5*(1.0+cos(PI*((t-Tstart)/t_tuke-1.0)));
        if(t > Tend-t_tuke && t < Tend)
        {
        x = 0.5*(1.0-cos(PI*(t-Tend)/t_tuke));
        }
        if(t > Tend) x = 0.0;
        
        kdotx = t-xi[n];
        
        Amp = x*AF[n];
        
        //  - FF[n]/fstar  Derivation says this is needed in the phase. Doesn't seem to be.
        
        Phase = -2.0*PI*f*kdotx-PF[n];
        
        RR = FpAR[n]*Aplus - FcAI[n]*Across;
        II = FcAR[n]*Across + FpAI[n]*Aplus;
        
        pa = atan2(II,RR);
        if(pa < 0.0) pa += TPI;
        
        if(pa-paold > 6.0) ja -= TPI;
        if(paold-pa > 6.0) ja += TPI;
        paold = pa;
        
        AAmp[n] = Amp*sqrt(RR*RR+II*II);
        APhase[n] = Phase+pa+ja;
        
        RR = FpER[n]*Aplus - FcEI[n]*Across;
        II = FcER[n]*Across + FpEI[n]*Aplus;
        
        pe = atan2(II,RR);
        if(pe < 0.0) pe += TPI;
        
        if(pe-peold > 6.0) je -= TPI;
        if(peold-pe > 6.0) je += TPI;
        peold = pe;
        
        EAmp[n] = Amp*sqrt(RR*RR+II*II);
        EPhase[n] = Phase+pe+je;
        
        //printf("%d %e %e\n", n, t, f);
        
        //fprintf(out,"%d %e %e %e %e %e\n", n, t, f, kdotx, Amp, Phase);
        
    }
    //fclose(out);
    
    //end = clock();
    //cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    //printf("response loop took %f seconds\n", cpu_time_used);
    
    
    /*   Deallocate Arrays   */
    
    free(xi);
    
    free(FcAR);
    free(FpAR);
    free(FcER);
    free(FpER);
    free(FcAI);
    free(FpAI);
    free(FcEI);
    free(FpEI);
    
    return;
}



void FisherFast(int ll, double *params, double **Fisher, double *SAS, double *SES, double Tobs, double Tzero, int N)
{
    double *AAmp, *EAmp, *APhase, *EPhase, *TF, *FF, *PF, *AF;
    double *PFref, *AFref, *TFref;
    double af, fr, df, DT, fac, deltaF, f, fmin, x, t;
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double dt, fny;
    double Amp, Phase, tf, fguess;
    double m1_SI, m2_SI, distance, tc, phic;
    double *paramsP, *paramsM;
    int i, j, k, ii, NF, NW, lx;
    double A, P, Tend;
    double px, fnew, kxm;
    double *Aint, *Eint, *SA, *SE;
    double *AphaseP, *AphaseM, *AampP, *AampM;
    double *EphaseP, *EphaseM, *EampP, *EampM;
    double **DphaseA, **DampA;
    double **DphaseE, **DampE;
    double *epsilon;
    
    int NFmax = 5000;
    
    int kmin, kmax, flag;
    
    dt = Tobs/(double)(N);
    fny = 1.0/(2.0*dt);
    Tend = Tzero+Tobs;

    epsilon = (double*)malloc(sizeof(double)* (NP));
    paramsP = (double*)malloc(sizeof(double)* (NP));
    paramsM = (double*)malloc(sizeof(double)* (NP));
    
    AF = (double*)malloc(sizeof(double)* (NFmax));
    PF = (double*)malloc(sizeof(double)* (NFmax));
    FF = (double*)malloc(sizeof(double)* (NFmax));
    TF = (double*)malloc(sizeof(double)* (NFmax));
    
    // Note: This call sets the FF array, which then gets held and used even when the parameters change a little
    SetUp(ll, params, fny, NFmax, &NF, FF, TF, PF, AF, Tobs, Tzero);
    
    TFref = (double*)malloc(sizeof(double)* (NF));
    PFref = (double*)malloc(sizeof(double)* (NF));
    AFref = (double*)malloc(sizeof(double)* (NF));
    
    for (i=0; i< NF; i++)
     {
      TFref[i] = TF[i];
      PFref[i] = PF[i];
      AFref[i] = AF[i];
      //printf("%d %e %e %e\n", i, FF[i], TF[i], AF[i]);
     }
    
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
    // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
    
    // will take log derivatives for m1, m2, DL
    
    epsilon[0] = 1.0e-6;
    epsilon[1] = 1.0e-6;
    epsilon[2] = 1.0e-4;
    epsilon[3] = 1.0e-4;
    epsilon[4] = 1.0e-4;
    epsilon[5] = 1.0e-1;
    epsilon[6] = 1.0e-4;
    epsilon[7] = 1.0e-4;
    epsilon[8] = 1.0e-4;
    epsilon[9] = 1.0e-4;
    epsilon[10] = 1.0e-4;
    
    AphaseP = (double*)malloc(sizeof(double)* (NF));
    AphaseM = (double*)malloc(sizeof(double)* (NF));
    
    AampP = (double*)malloc(sizeof(double)* (NF));
    AampM = (double*)malloc(sizeof(double)* (NF));
    
    EphaseP = (double*)malloc(sizeof(double)* (NF));
    EphaseM = (double*)malloc(sizeof(double)* (NF));
    
    EampP = (double*)malloc(sizeof(double)* (NF));
    EampM = (double*)malloc(sizeof(double)* (NF));
    
    DphaseA = double_matrix(NP, NF);
    DampA = double_matrix(NP, NF);
    
    DphaseE = double_matrix(NP, NF);
    DampE = double_matrix(NP, NF);
 
    
    AAmp = (double*)malloc(sizeof(double)* (NF));
    EAmp = (double*)malloc(sizeof(double)* (NF));
    APhase = (double*)malloc(sizeof(double)* (NF));
    EPhase = (double*)malloc(sizeof(double)* (NF));
    
    Aint = (double*)malloc(sizeof(double)* (NF));
    Eint = (double*)malloc(sizeof(double)* (NF));
    SA = (double*)malloc(sizeof(double)* (NF));
    SE = (double*)malloc(sizeof(double)* (NF));
    
    // reference amplitude and phase
    Extrinsic(params, Tend, NF, FF, TF, PF, AF, AAmp, EAmp, APhase, EPhase, &kxm);
    
    for (k=0; k< NF; k++)
    {
        i = (int)(FF[k]*Tobs);
        SA[k] = 1.0e-40;
        SE[k] = 1.0e-40;
        if(i > -1 && i < N/2)
        {
        SA[k] = SAS[i];
        SE[k] = SES[i];
        }
      // printf("%d %e %e %e\n", k, FF[k], SA[k], SE[k]);
    }
    
    // masses and spins
    for (i=0; i< 4; i++)
    {
        
        for (j=0; j< NP; j++)
        {
            paramsP[j] = params[j];
            paramsM[j] = params[j];
        }
        
        
        if(ll==0)
        {
            if(i==0 || i==1) // log derivatives
            {
            x = exp(epsilon[i]);
            paramsP[i] *= x;
            paramsM[i] /= x;
            }
            else
            {
                paramsP[i] += epsilon[i];
                paramsM[i] -= epsilon[i];
            }
            
        }
        else
        {
            paramsP[i] += epsilon[i];
            paramsM[i] -= epsilon[i];
        }

            Intrinsic(ll,paramsP, NF, FF, TF, PF, AF, Tobs, Tzero);  // intrinsic
            Extrinsic(paramsP, Tend, NF, FF, TF, PF, AF, AampP, EampP, AphaseP, EphaseP, &kxm); // extrinsic
                
            Intrinsic(ll,paramsM, NF, FF, TF, PF, AF, Tobs, Tzero); // intrinsic
            Extrinsic(paramsM, Tend, NF, FF, TF, PF, AF, AampM, EampM, AphaseM, EphaseM, &kxm); // extrinsic

            for (k=0; k< NF; k++)
            {
                x = (AphaseP[k]-AphaseM[k]);
                if(x > 0.0) while(x > PI) x -= PI;
                if(x < 0.0) while(x < -PI) x += PI;
                DphaseA[i][k] = x;
                x = (EphaseP[k]-EphaseM[k]);
                if(x > 0.0) while(x > PI) x -= PI;
                if(x < 0.0) while(x < -PI) x += PI;
                DphaseE[i][k] = x;
                DampA[i][k] = (AampP[k]-AampM[k]);
                DampE[i][k] = (EampP[k]-EampM[k]);
            }
        
    }

    

   for (i=4; i< 7; i++)
    {
    
      if(i == 4)  // Put in coallesence phase manually. Note this is orbital phase, so times two in GW phase
       {
           
        x = -4.0*epsilon[4];
        for (k=0; k< NF; k++)
        {
            DphaseA[i][k] = x;
            DphaseE[i][k] = x;
            DampA[i][k] = 0.0;
            DampE[i][k] = 0.0;
        }
        
      }
     else if(i == 5)  // tc is not included in the phase subroutine
      {
        for (k=0; k< NF; k++)
        {
            x = -4.0*PI*FF[k]*epsilon[5];
            DphaseA[i][k] = x;
            DphaseE[i][k] = x;
            DampA[i][k] = 0.0;
            DampE[i][k] = 0.0;
        }
        
       }
      else if(i == 6)  // ln DL derivative
       {
           for (k=0; k< NF; k++)
           {
               DphaseA[i][k] = 0;
               DphaseE[i][k] = 0;
               DampA[i][k] = -2.0*epsilon[6]*AAmp[k];
               DampE[i][k] = -2.0*epsilon[6]*EAmp[k];
           }

       }
        

  }

    for (i=7; i< NP; i++)
    {
    
       for (j=0; j< NP; j++)
        {
        paramsP[j] = params[j];
        paramsM[j] = params[j];
        }
    
        paramsP[i] += epsilon[i];
        paramsM[i] -= epsilon[i];

        Extrinsic(paramsP, Tend, NF, FF, TFref, PFref, AFref, AampP, EampP, AphaseP, EphaseP, &kxm); // extrinsic
        Extrinsic(paramsM, Tend, NF, FF, TFref, PFref, AFref, AampM, EampM, AphaseM, EphaseM, &kxm); // extrinsic
        
        for (k=0; k< NF; k++)
        {
            x = (AphaseP[k]-AphaseM[k]);
            if(x > 0.0) while(x > PI) x -= PI;
            if(x < 0.0) while(x < -PI) x += PI;
            DphaseA[i][k] = x;
            
            x = (EphaseP[k]-EphaseM[k]);
            if(x > 0.0) while(x > PI) x -= PI;
            if(x < 0.0) while(x < -PI) x += PI;
            DphaseE[i][k] = x;
            
            DampA[i][k] = (AampP[k]-AampM[k]);
            DampE[i][k] = (EampP[k]-EampM[k]);
        }
        
    }

   


    gsl_interp_accel *IAacc = gsl_interp_accel_alloc();
    gsl_spline *IAspline = gsl_spline_alloc (gsl_interp_cspline, NF);


    gsl_interp_accel *IEacc = gsl_interp_accel_alloc();
    gsl_spline *IEspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    
    
   
    for (i=0; i< NP; i++)
    {
        for (j=i; j< NP; j++)
        {
  
           // printf("%d %d\n", i, j);
            
            for (k=0; k< NF; k++)
            {
               // printf("%d %d %e %e %e\n", i, k, DampA[i][k], AAmp[k], DphaseA[i][k]);
               // printf("%d %d %e %e %e\n", i, k, DampE[i][k], EAmp[k], DphaseE[i][k]);
            Aint[k] = 4.0*(DampA[i][k]*DampA[j][k]+AAmp[k]*AAmp[k]*DphaseA[i][k]*DphaseA[j][k])/SA[k];
            Eint[k] = 4.0*(DampE[i][k]*DampE[j][k]+EAmp[k]*EAmp[k]*DphaseE[i][k]*DphaseE[j][k])/SE[k];
                //printf("%d %d %e\n", i, k, Aint[k]);
               // printf("%d %d %e\n", i, k, Eint[k]);
            }
            
            gsl_spline_init(IAspline, FF, Aint, NF);
            gsl_spline_init(IEspline, FF, Eint, NF);
            
            Fisher[i][j] = gsl_spline_eval_integ(IAspline, FF[0], FF[NF-1], IAacc)+gsl_spline_eval_integ(IEspline, FF[0], FF[NF-1], IEacc);
            
            Fisher[i][j] /= (4.0*epsilon[i]*epsilon[j]/Tobs);
        
        }
        
    }
    
    for (i=0; i< NP; i++)
    {
        for (j=i+1; j< NP; j++)
        {
            Fisher[j][i] =  Fisher[i][j];
        }
    }
    
    
    free_double_matrix(DphaseA,NP);
    free_double_matrix(DampA,NP);
    free_double_matrix(DphaseE,NP);
    free_double_matrix(DampE,NP);
 
    free(epsilon);
    free(SA);
    free(SE);
    free(Aint);
    free(Eint);
    free(AAmp);
    free(EAmp);
    free(APhase);
    free(EPhase);
    free(FF);
    free(TF);
    free(AF);
    free(PF);
    free(TFref);
    free(AFref);
    free(PFref);
    
    free(paramsP);
    free(paramsM);
    
    free(AphaseP);
    free(AphaseM);
    free(AampP);
    free(AampM);
    
    free(EphaseP);
    free(EphaseM);
    free(EampP);
    free(EampM);

    gsl_spline_free(IAspline);
    gsl_spline_free(IEspline);

    gsl_interp_accel_free(IAacc);
    gsl_interp_accel_free(IEacc);
 
    
    
}


void FisherMax(int ll, double *params, double **Fisher, double *SN, double Tobs, double Tzero, int N)
{
    double *TF, *FF, *PF, *AF;
    double af, fr, df, DT, fac, deltaF, f, fmax, fmin, x, y, t;
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double Amp, Phase, tf, fguess;
    double m1_SI, m2_SI, distance, tc, phic;
    double *paramsP, *paramsM;
    int i, j, k, ii, NF, NW, lx;
    double A, P;
    double px, fnew, kxm;
    double *Xint, *SAE;
    double *phaseP, *phaseM, *ampP, *ampM;
    double **Dphase, **Damp;
    double *epsilon;
    double dt, fny;
    
    dt = Tobs/(double)(N);
    fny = 1.0/(2.0*dt);
    
    int NFmax = 10000;
    
    int kmin, kmax, flag;
    
    epsilon = (double*)malloc(sizeof(double)* (NP));
    paramsP = (double*)malloc(sizeof(double)* (NP));
    paramsM = (double*)malloc(sizeof(double)* (NP));
    
    AF = (double*)malloc(sizeof(double)* (NFmax));
    PF = (double*)malloc(sizeof(double)* (NFmax));
    FF = (double*)malloc(sizeof(double)* (NFmax));
    TF = (double*)malloc(sizeof(double)* (NFmax));
    
    //for (j=0; j< NV; j++) printf("%e ", params[j]);
   // printf("eta = %f\n", pow(exp(params[0])/exp(params[1]), (5.0/3.0)));
    
    // Note: This call sets the FF array, which then gets held and used even when the parameters change a little
    SetUp(ll, params, fny, NFmax, &NF, FF, TF, PF, AF, Tobs, Tzero);
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
    // [7] E/A amplitude [8] PhaseA-PhaseE
    
    // combines A annd E channel amplitudes
    fac = (1.0+params[7]*params[7]);
    
    // will take log derivatives for Mc, Mt
    
    epsilon[0] = 1.0e-6;
    epsilon[1] = 1.0e-6;
    epsilon[2] = 1.0e-4;
    epsilon[3] = 1.0e-4;
    epsilon[4] = 1.0e-4;
    epsilon[5] = 1.0e-1;
    
    phaseP = (double*)malloc(sizeof(double)* (NF));
    phaseM = (double*)malloc(sizeof(double)* (NF));
    
    ampP = (double*)malloc(sizeof(double)* (NF));
    ampM = (double*)malloc(sizeof(double)* (NF));
    
    Dphase = double_matrix(NV, NF);
    Damp = double_matrix(NV, NF);
    
    Xint = (double*)malloc(sizeof(double)* (NF));
    SAE = (double*)malloc(sizeof(double)* (NF));
    
    
    // masses and spins
    for (i=0; i< NV; i++)
    {
        
        for (j=0; j< NP; j++)
        {
            paramsP[j] = params[j];
            paramsM[j] = params[j];
        }
        
        
        if(ll==0)
        {
            if(i==0 || i==1) // log derivatives
            {
            x = exp(epsilon[i]);
            paramsP[i] *= x;
            paramsM[i] /= x;
            }
            else
            {
                paramsP[i] += epsilon[i];
                paramsM[i] -= epsilon[i];
            }
            
        }
        else
        {
            paramsP[i] += epsilon[i];
            paramsM[i] -= epsilon[i];
        }


        
            if(i == 4)  // Put in coallesence phase manually. Note this is orbital phase, so times two in GW phase
              {
                  
               x = -4.0*epsilon[4];
               for (k=0; k< NF; k++)
               {
                   Dphase[i][k] = x;
                   Damp[i][k] = 0.0;
               }
               
             }
            else if(i == 5)  // tc is not included in the phase subroutine
             {
               for (k=0; k< NF; k++)
               {
                   x = -4.0*PI*FF[k]*epsilon[5];
                   Dphase[i][k] = x;
                   Damp[i][k] = 0.0;
               }
             }
            else
            {
                
            Intrinsic(ll,paramsP, NF, FF, TF, phaseP, ampP, Tobs, Tzero);
            Intrinsic(ll,paramsM, NF, FF, TF, phaseM, ampM, Tobs, Tzero);
            
            for (k=0; k< NF; k++)
            {
                Dphase[i][k] = (phaseP[k]-phaseM[k]);
                Damp[i][k] = (ampP[k]-ampM[k]);
            }
                
            }
        
    }

    
    //printf("%d %d\n", (int)(FF[0]*Tobs), (int)(FF[NF-1]*Tobs));
    for (k=0; k< NF; k++)
    {
        i = (int)(FF[k]*Tobs);
        SAE[k] = 1.0e-40;
        if(i > -1 && i < N/2) SAE[k] = SN[i];
    }

    gsl_interp_accel *Iacc = gsl_interp_accel_alloc();
    gsl_spline *Ispline = gsl_spline_alloc (gsl_interp_cspline, NF);

   
    for (i=0; i< NV; i++)
    {
        for (j=i; j< NV; j++)
        {
  
            for (k=0; k< NF; k++)
            {
            Xint[k] = 4.0*fac*(Damp[i][k]*Damp[j][k]+AF[k]*AF[k]*Dphase[i][k]*Dphase[j][k])/SAE[k];
            }
            
            gsl_spline_init(Ispline, FF, Xint, NF);
    
            Fisher[i][j] = gsl_spline_eval_integ(Ispline, FF[0], FF[NF-1], Iacc);
            
            Fisher[i][j] /= (4.0*epsilon[i]*epsilon[j]/Tobs);
        }
        
    }
    
    for (i=0; i< NV; i++)
    {
        for (j=i+1; j< NV; j++)
        {
            Fisher[j][i] =  Fisher[i][j];
        }
    }
    
    /*
    printf("\n");
    for (i=0; i< NV; i++)
    {
        for (j=0; j< NV; j++)
        {
            printf("%e ", Fisher[i][j]);
        }
        printf("\n");
    }
    printf("\n"); */
    
    // stabilizer
    for (i=0; i< NV; i++) Fisher[i][i] += 4.0;

    free_double_matrix(Dphase,NV);
    free_double_matrix(Damp,NV);
 
    free(epsilon);
    free(SAE);
    free(Xint);
    free(FF);
    free(TF);
    free(AF);
    free(PF);
    
    free(paramsP);
    free(paramsM);
    
    free(phaseP);
    free(phaseM);
    free(ampP);
    free(ampM);

    gsl_spline_free(Ispline);
    gsl_interp_accel_free(Iacc);
   
    
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
        gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (d);
        
        ecc = gsl_eigen_symmv (m, eval, evec, w);
        
        gsl_eigen_symmv_free (w);
        
        
        sc = gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
        
        for (i = 0; i < d; i++)
        {
            ej[i] = gsl_vector_get (eval, i);
        
            for (j = 0 ; j < d ; j++)
            {
                ev[i][j] = gsl_matrix_get(evec, j, i);
            }
            
        }
        
        for (i = 0; i < d; i++)
        {
           // printf("%d %e ", i, ej[i]);
            // turn into 1-sigma jump amplitudes
            ej[i] = 1.0/sqrt(fabs(ej[i]));
          //  printf("%e\n", ej[i]);
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

void StartStop(int ll, double *params, double Tseg, double tstart, double tstop, double *fstart, double *fstop, double *frg)
{
    
    double m1, m2, m1_SI, m2_SI, chi1, chi2, tc, dm;
    double Mtot, eta, Mc, af, fr;
    double Amp, Phase, distance;
    double fmin, fmax;
    double fnew, tf, fny;
    int i;
    
    if(ll == 0)  // linear in m1, m2, DL
    {
    m1 = params[0]*TSUN;    // Mass1
    m2 = params[1]*TSUN;    // Mass2
    distance = params[6]*1.0e9*PC_SI; // distance
    Mtot = (m1+m2);
    eta = m1*m2/((m1+m2)*(m1+m2));
    Mc = pow(m1*m2,3.0/5.0)/pow(m1+m2,1.0/5.0);
    }
    else if(ll == 1)  // log in m1, m2, DL
    {
     m1 = exp(params[0])*TSUN;    // Mass1
     m2 = exp(params[1])*TSUN;    // Mass2
     distance = exp(params[6])*1.0e9*PC_SI; // distance
     Mtot = (m1+m2);
     eta = m1*m2/((m1+m2)*(m1+m2));
     Mc = pow(m1*m2,3.0/5.0)/pow(m1+m2,1.0/5.0);
    }
    else // log in Mc, Mt, DL
    {
    distance = exp(params[6])*1.0e9*PC_SI; // distance
    Mc = exp(params[0])*TSUN;
    Mtot = exp(params[1])*TSUN;
    eta = pow((Mc/Mtot), (5.0/3.0));
     if(eta > 0.25)
     {
        dm = 0.0;
     }
     else
     {
        dm = sqrt(1.0-4.0*eta);
     }
    m1 = Mtot*(1.0+dm)/2.0;
    m2 = Mtot*(1.0-dm)/2.0;
    }
    
    m1_SI = m1*MSUN_SI/TSUN;
    m2_SI = m2*MSUN_SI/TSUN;
    
    chi1 = params[2];
    chi2 = params[3];
    tc = params[5];
    
    af = FinalSpin0815(eta, chi1, chi2);
    fr = fring(eta, chi1, chi2, af)/Mtot;
    
    if(fr < 0.0) fr = 1.0/Mtot;  // should not happen ....
    
    *frg = fr;
    
    // Here we space the frequency array to give approximately equal spacing in time
    // The dynamic frequency spacing is capped to be between 1e-6 and 1e-4 Hz

    fmin = f_at_t(m1, m2, chi1, chi2, tc, tstart);

    // nan catcher
    if(fmin != fmin) fmin = 1.0/Tseg;
    if(fmin < 0.0) fmin = 1.0/Tseg;
    if(fmin > fr) fmin = 0.9*fr;  // shouldn't happen
    
    fmax = f_at_t(m1, m2, chi1, chi2, tc, tstop);
    
    // nan catcher
    if(fmax != fmax) fmax = 2.0*fr;
    if(fmax > 2.0*fr) fmax = 2.0*fr; // shouldn't happen
    if(fmax < fmin) fmax = 2.0*fmin;
    
    //printf("%e %e %e\n", fmin, fmax, fr);
    
    *fstart = fmin;
    *fstop = fmax;
    
}


void SetUp(int ll, double *params, double fny, int NFmax, int *NFS, double *FF, double *TF, double *PF, double *AF, double Tobs, double Tzero)
{
    double af, fr, df, DT, fac, deltaF, f, fmax, fmin, x;
    double m1, m2, chi1, chi2, Mtot, Mc, eta, dm;
    double Amp, Phase, tf, fguess;
    double m1_SI, m2_SI, distance, tc, phic, told;
    int i, ii, NF;
    double A, P;
    double px, fnew, fonfs;
    double t, tx;
    double dfmin, dfmax;
    
    double fRef_in=PDfref;
    
    int ret, flag, flag1, flag2;
    
    // This subroutine sets up the frequency sample array, the time-frequency map and finds the PhenomD amplitude and phase
    
    // NFmax is the size of the holder arrays. NFS is the actual size.
    
    if(ll == 0)  // linear in m1, m2
       {
       m1 = params[0]*TSUN;    // Mass1
       m2 = params[1]*TSUN;    // Mass2
       distance = params[6]*1.0e9*PC_SI; // distance
       Mtot = (m1+m2);
       eta = m1*m2/((m1+m2)*(m1+m2));
       Mc = pow(m1*m2,3.0/5.0)/pow(m1+m2,1.0/5.0);
       }
       else if(ll == 1)  // log in m1, m2
       {
        m1 = exp(params[0])*TSUN;    // Mass1
        m2 = exp(params[1])*TSUN;    // Mass2
        distance = exp(params[6])*1.0e9*PC_SI; // distance
        Mtot = (m1+m2);
        eta = m1*m2/((m1+m2)*(m1+m2));
        Mc = pow(m1*m2,3.0/5.0)/pow(m1+m2,1.0/5.0);
       }
       else // log in Mc, Mt
       {
       distance = exp(params[6])*1.0e9*PC_SI; // distance
       Mc = exp(params[0])*TSUN;
       Mtot = exp(params[1])*TSUN;
       eta = pow((Mc/Mtot), (5.0/3.0));
        if(eta > 0.25)
        {
           dm = 0.0;
        }
        else
        {
           dm = sqrt(1.0-4.0*eta);
        }
       m1 = Mtot*(1.0+dm)/2.0;
       m2 = Mtot*(1.0-dm)/2.0;
       }
       
       m1_SI = m1*MSUN_SI/TSUN;
       m2_SI = m2*MSUN_SI/TSUN;
    
    StartStop(ll, params, Tobs, Tzero, Tzero+Tobs, &fmin, &fmax, &fr);
    
     // printf("%e %e\n", fmin, fmax);
    
    if(fmin > fny) fmin = 0.9*fny;
    if(fmax > fny) fmax = fny;
    
    // this can happen when tc is really small and the masses are large
    if(fmax < fmin)
    {
        fmin = 0.5*fmax;
    }

 
    dfmin = 1.0/Tobs;
    dfmax = fmax/100.0;
    DT = 3.0e5;
    
    fac = DT*pow(8.0*PI, 8.0/3.0)*3.0/40.0*pow(Mc,5.0/3.0);
    
    f = fmin;
    NF = 1;
    do
    {
        df = fac*pow(f,11.0/3.0);
        if(df < dfmin) df = dfmin;
        if(df > dfmax) df = dfmax;
        f += df;
        NF++;
    }while(f < fmax);
    
    //printf("%d\n", NF);
    
     // Need to catch is NF > NFmax
    
    if(NF > NFmax)
    {
        NF = NFmax;
        df = (fmax-fmin)/(double)(NFmax-1);
        for (i=0; i< NF; i++)
        {
            FF[i] = fmin + df*(double)(i);
        }
    }
    else
    {
    
    f = fmin;
    FF[0] = fmin;
    for (i=1; i< NF; i++)
    {
        df = fac*pow(f,11.0/3.0);
        if(df < dfmin) df = dfmin;
        if(df > dfmax) df = dfmax;
        f += df;
        FF[i] = f;
    }
        
    }
    

    
    if(NF < 4)
    {
        NF = 4;
        df = (fmax-fmin)/3.0;
        FF[0] = fmin;
        for (i=1; i< NF; i++) FF[i] = fmin+df*(double)(i);
    }
    
    flag = 0;
    if(FF[0] != FF[0]) flag = 1;
    for (i=1; i< NF; i++)
    {
        if(FF[i] != FF[i]) flag = 1;
        if(FF[i] <= FF[i-1]) flag = 1;
    }
    
    // to stop the code from crashing - during the search wierd stuff can happen
    if(flag == 1)
    {
        NF = 10;
        df = fny/10.0;
        if(df < dfmin) df = dfmin;
        FF[0] = 1.0/Tobs;
        for (i=1; i< NF; i++) FF[i] = FF[i-1]+df;
    }
    
    //printf("%d %d\n", flag, NF);
    

    Intrinsic(ll, params, NF, FF, TF, PF, AF, Tobs, Tzero);

    *NFS = NF;
    
}

void Intrinsic(int ll, double *params, int NF, double *FF, double *TF, double *PF, double *AF, double Tobs, double Tzero)
{
    
    AmpPhaseFDWaveform *ap = NULL;
    RealVector *freq;
    
    double m1, m2, chi1, chi2, Mtot, Mc, eta, dm;
    double m1_SI, m2_SI, distance, tc, phic;
    
    double af, fr, fx, dtdf;
    
    int i, ii, ret, flag;
    double fonfs, t, told, tx;
    
    double sqrtTobs;
    
    double fRef_in=PDfref;
    
    sqrtTobs = sqrt(Tobs);
    
    if(ll == 0)  // linear in m1, m2, DL
    {
    m1 = params[0]*TSUN;    // Mass1
    m2 = params[1]*TSUN;    // Mass2
    distance = params[6]*1.0e9*PC_SI; // distance
    Mtot = (m1+m2);
    eta = m1*m2/((m1+m2)*(m1+m2));
    Mc = pow(m1*m2,3.0/5.0)/pow(m1+m2,1.0/5.0);
    }
    else if(ll == 1)  // log in m1, m2, DL
    {
     m1 = exp(params[0])*TSUN;    // Mass1
     m2 = exp(params[1])*TSUN;    // Mass2
     distance = exp(params[6])*1.0e9*PC_SI; // distance
     Mtot = (m1+m2);
     eta = m1*m2/((m1+m2)*(m1+m2));
     Mc = pow(m1*m2,3.0/5.0)/pow(m1+m2,1.0/5.0);
    }
    else // log in Mc, Mt, DL
    {
    distance = exp(params[6])*1.0e9*PC_SI; // distance
    Mc = exp(params[0])*TSUN;
    Mtot = exp(params[1])*TSUN;
    eta = pow((Mc/Mtot), (5.0/3.0));
     if(eta > 0.25)
     {
        dm = 0.0;
     }
     else
     {
        dm = sqrt(1.0-4.0*eta);
     }
    m1 = Mtot*(1.0+dm)/2.0;
    m2 = Mtot*(1.0-dm)/2.0;
    }
    
    m1_SI = m1*MSUN_SI/TSUN;
    m2_SI = m2*MSUN_SI/TSUN;
    
    chi1 = params[2];  // Spin1
    chi2 = params[3];  // Spin2
    tc = params[5];    // merger time
    

    
    af = FinalSpin0815(eta, chi1, chi2);
    fr = fring(eta, chi1, chi2, af)/Mtot;
    
    dtdf = 0.1*Mtot/fr; // how we advance time after PhenomD freezes
    
    freq = CreateRealVector(NF);
    
    for (i=0; i< NF; i++) freq->data[i] = FF[i];
    
    // ther maximization can mess with the distances
    if(distance != distance)
    {
       // printf("nan %e\n", params[6]);
        distance = 1.0e9*PC_SI;
    }
    if(distance <= 0.0)
    {
        printf("<= 0 %e\n", params[6]);
        distance = 1.0e9*PC_SI;
    }

    
    ret = IMRPhenomDGenerateh22FDAmpPhase(
                                          &ap,      /**< [out] FD waveform */
                                          freq,    /**< Input: frequencies (Hz) on which to evaluate h22 FD */
                                          0.0,                  /**< Orbital phase at fRef (rad) */
                                          fRef_in,               /**< reference frequency (Hz) */
                                          m1_SI,                 /**< Mass of companion 1 (kg) */
                                          m2_SI,                 /**< Mass of companion 2 (kg) */
                                          chi1,                  /**< Aligned-spin parameter of companion 1 */
                                          chi2,                  /**< Aligned-spin parameter of companion 2 */
                                          distance               /**< Distance of source (m) */
                                          );
    
    
    flag = 0;
    told = ap->time[0]+tc;
    
    for (i=0; i< NF; i++)
    {
        PF[i] = ap->phase[i];
        
        AF[i] =  h22fac*ap->amp[i]/sqrtTobs;
        fonfs = freq->data[i]/fstar;
        AF[i] *= (4.0*fonfs*sin(fonfs));   // conversion to fractional frequency and leading order TDI transfer function
        
        t = ap->time[i]+tc;
        
        // catch the turn over in f(t)
        if(t < told && flag == 0)
        {
            flag = 1;
            ii = i-1;
            tx = told;
            fx = freq->data[ii];
        }
        
        TF[i] = t;
        
        // keep the clock ticking forward
        if(flag == 1) TF[i] = tx + (freq->data[i]-fx)*dtdf;
        told = t;
    }
    
    // printf("Tzero %e %e %e\n", TF[0], tc, tc+ap->time[0]);
    
    DestroyAmpPhaseFDWaveform(ap);
    DestroyRealVector(freq);

    
}

double log_likelihood_APmax_dual(int ll, double *A, double *E, double *params, double *SA, double *SE, int N, double Tobs, double Tzero)
{
    int i, j, k, NMAX;
    int ii, jj, Nend;
    double sum;
    double HHA, HHE, LL;
    double logL;
    double logLfs;
    double fmax;
    double HA, HE, LD, dt, x;
    double normA, normE, deltH, pshiftA, pshiftE;
    double cA, sA, cE, sE;
    double *h;
    
    //Chi^2 = (d-h|d-h)/Sn
    
    dt = Tobs/(double)(N);
    
    h = (double*)malloc(sizeof(double)*(N));
    
    PDwave(ll, h, params, N, Tobs, Tzero);
    
    
    HHA = fourier_nwip(h, h, SA, N);
    HHE = fourier_nwip(h, h, SE, N);
    
    logL = 0.0;
    
    if(HHA > 1.0e-8 && HHE > 1.0e-8)
    {

    quadratures(&cA, &sA, A, h, SA, N);
    quadratures(&cE, &sE, E, h, SE, N);
    
    HA = sqrt(cA*cA+sA*sA);
    HE = sqrt(cE*cE+sE*sE);

   // printf("new %e %e %e %e\n", HA, HE, HHA, HHE);
    
    pshiftA = atan2(sA,cA);
    pshiftE = atan2(sE,cA);
    normA = HA/HHA;
    normE = HE/HHE;
    
    // reference to channel A
    if(ll==0)
    {
     params[6] /= normA;
    }
    else
    {
     params[6] -= log(normA);
    }
    
    params[4] -= pshiftA/2.0;  //
    
    // [7] E/A amplitude [8] PhaseA-PhaseE
    
    params[7] = normE/normA;
    params[8] = 0.5*(pshiftA-pshiftE);
    
    logL = (HA*HA)/(2.0*HHA)+(HE*HE)/(2.0*HHE);
        
    // catch failed maximization
    if(params[6] != params[6] || params[4] != params[4] || params[7] != params[7] || params[8] != params[8]) logL = -1.0e20;
        
    }
    
    free(h);
    
    return logL;
}




double log_likelihood_max_dual(int ll, double *A, double *E, double *params, double *SA, double *SE, int N, double Tobs, double Tzero)
{
    int i, j, k, NMAX;
    int ii, jj, Nend;
    double sum;
    double HHA, HHE, LL;
    double logL;
    double logLfs;
    double fmax;
    double HA, HE, LD, dt, x, y, t;
    double normA, normE, deltH, pshiftA, pshiftE;
    double *AS, *ES;
    double *AC, *AF;
    double *EC, *EF;
    double *h;
    
    //Chi^2 = (d-h|d-h)/Sn
    
    dt = Tobs/(double)(N);
    
    h = (double*)malloc(sizeof(double)*(N));
    
    PDwave(ll, h, params, N, Tobs, Tzero);
    
    AS = double_vector(N); ES = double_vector(N);
    AC=double_vector(N);  AF=double_vector(N);
    EC=double_vector(N);  EF=double_vector(N);
    
    HHA = fourier_nwip(h, h, SA, N);
    HHE = fourier_nwip(h, h, SE, N);
    
    logL = 0.0;
    
    if(HHA > 1.0e-8 && HHE > 1.0e-8)
    {
        
    pbt_shift(AC, AF, A, h, SA, N);
    pbt_shift(EC, EF, E, h, SE, N);
    
    // only allow time shifts up to +/- Tobs/8 (otherwise frequency range can be way off)
    Nend = N/8;
    
    for(i = 0; i < Nend/2; i++) AS[i+N/2] = sqrt(AC[i]*AC[i]+AF[i]*AF[i]);
    for(i = -Nend/2; i < 0; i++) AS[i+N/2] = sqrt(AC[N+i]*AC[N+i]+AF[N+i]*AF[N+i]);
    
    for(i = 0; i < Nend/2; i++) ES[i+N/2] = sqrt(EC[i]*EC[i]+EF[i]*EF[i]);
    for(i = -Nend/2; i < 0; i++) ES[i+N/2] = sqrt(EC[N+i]*EC[N+i]+EF[N+i]*EF[N+i]);
    
    x = 0.0;
    k = 0;
    for (i = -Nend/2; i < Nend/2; ++i)
    {
        j = i+N/2;
        t = params[5]+(double)(i)*dt;
        
        y = (t-Tzero)/Tobs;
        
        if(fabs(y-rint(y)) > buffer)
        {
         if((AS[j]+ES[j]) > x)
         {
            x = AS[j]+ES[j];
            k = i;
         }
        }
    }
        
    HA = 0.0;
    HE = 0.0;
    
        if(x > 0.0)
        {
         HA = 2.0*(double)(N)*(AS[k+N/2]);
         HE = 2.0*(double)(N)*(ES[k+N/2]);
        }
    
    deltH = dt*(double)(k);
    
    if(k < 0) k += N;
    

     // printf("old %e %e %e %e\n", HA, HE, HHA, HHE);

    pshiftA = atan2(AF[k],AC[k]);
    pshiftE = atan2(EF[k],EC[k]);
    normA = HA/HHA;
    normE = HE/HHE;
    
    // reference to channel A
    if(ll==0)
    {
     params[6] /= normA;
    }
    else
    {
     params[6] -= log(normA);
    }
    params[5] += deltH;
    params[4] -= pshiftA/2.0;  //
    
    // [7] E/A amplitude [8] PhaseA-PhaseE
    
    params[7] = normE/normA;
    params[8] = 0.5*(pshiftA-pshiftE);
    
    logL = (HA*HA)/(2.0*HHA)+(HE*HE)/(2.0*HHE);
        
    // catch failed maximization
    if(params[6] != params[6] || params[4] != params[4] || params[7] != params[7] || params[8] != params[8]) logL = -1.0e20;
        
    }
    
    free_double_vector(AC);  free_double_vector(AF);
    free_double_vector(EC);  free_double_vector(EF);

    free_double_vector(AS);
    free_double_vector(ES);
    
    free(h);

    
    return logL;
}

double log_likelihood_max(int ll, double *H, double *params, double *SN, int N, double Tobs, double Tzero)
{
    int i, j, k, NMAX;
    int ii, jj;
    double sum;
    double HH, LL;
    double logL;
    double logLfs;
    double fmax;
    double HD, LD, dt, x;
    double normH, deltH, pshiftH;
    double *HS;
    double *HC, *HF;
    double *h;
    
    //Chi^2 = (d-h|d-h)/Sn
    
    dt = Tobs/(double)(N);
    
    h = (double*)malloc(sizeof(double)*(N));
    
    PDwave(ll, h, params, N, Tobs, Tzero);
    
    HS = double_vector(N);
    HC=double_vector(N);  HF=double_vector(N);
    
    pbt_shift(HC, HF, H, h, SN, N);
    
    for(i = 0; i < N/2; i++) HS[i+N/2] = sqrt(HC[i]*HC[i]+HF[i]*HF[i]);
    for(i = -N/2; i < 0; i++) HS[i+N/2] = sqrt(HC[N+i]*HC[N+i]+HF[N+i]*HF[N+i]);
    
    x = 0;
    for (i = -N/2; i < N/2; ++i)
    {
        if(HS[i+N/2] > x)
        {
            x = HS[i+N/2];
            k = i;
        }
    }
    
    HD = 2.0*(double)(N)*HS[k+N/2];
    
    deltH = dt*(double)(k);
    
    if(k < 0) k += N;
    
    pshiftH = atan2(HF[k],HC[k]);
    
    HH = fourier_nwip(h, h, SN, N);
    
    // Inverse FFTs in fourier_nwip are un-normlaized
   /* HH /= Tobs;
    HD /= Tobs; */

    normH = HD/HH;

    logL = (HD*HD/HH)/2.0;
    
    // Shift the waverform to match
    
    
    if(ll==0)
    {
     params[6] /= normH;
    }
    else
    {
     params[6] -= log(normH);
    }
    params[5] += deltH;
    params[4] -= pshiftH/2.0;  // PhenomD uses orbital phase, while here we have GW phase
    
    free_double_vector(HC);  free_double_vector(HF);

    free_double_vector(HS);
    
    free(h);

    
    return logL;
}

double likelihoodFstat(int ll, double *params, int N, double *AC, double *EC, double *SA, double *SE, double Tobs, double Tzero)
{
    
    /*   Indicies   */
    int i,j, k, n, m, imax;
    
    /*   Miscellaneous  */
    double xm, fstep, power, om, mx, tx;
    
    double *pfilt;
    
    double **filtA, **filtE;
    
    double **MM, **MI;
    
    double *KV, *aV;
    
    double AR, AI, ER, EI;
    
    double t, f, kdotx, Amp, Phase, Fp, Fc;
    
    double A, P, px, tc, tcs, pxi, pxj;
    
    double x, y, u, v, xold, yold;
    
    double cp, sp;
    
    int NA, flag;
    
    double cx, sx, logL, logLmax;
    
    int nn;
    
    double iota, cosi;
    double Ap, Ac, scale;
    
    double a1, a2, a3, a4;
    double cpsi, spsi;
    double psi, phic;
    double cpc, spc;
    
    FILE *out;
    
    pfilt = (double*)malloc(sizeof(double)* (NP));
    
    filtA = double_matrix(4,N);
    filtE = double_matrix(4,N);
    
    KV = (double*)malloc(sizeof(double)* (4));
    aV = (double*)malloc(sizeof(double)* (4));
    MM = double_matrix(4,4);
    MI = double_matrix(4,4);
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
    // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination

    ResponseFstat(ll, params, N, filtA, filtE, Tobs, Tzero);
    
        for (i=0; i< 4; i++)
        {
            KV[i] = 2.0*(fourier_nwip(filtA[i], AC, SA, N) + fourier_nwip(filtE[i], EC, SE, N));
            for (j=i; j< 4; j++)
            {
                MM[i][j] =  4.0*(fourier_nwip(filtA[i], filtA[j], SA, N) + fourier_nwip(filtE[i], filtE[j], SE, N));
               
            }
        }
        
        for (i=0; i< 4; i++)
        {
            for (j=i; j< 4; j++)
            {
                MM[j][i] = MM[i][j];
            }
        }
        
        // a small stabilizer
        for (i=0; i< 4; i++) MM[i][i] += 0.1;
        
        Inverse(MM, MI, 4);
        
        logL = 0.0;
        
        for (i=0; i< 4; i++)
        {
            aV[i] = 0.0;
            for (j=0; j< 4; j++)
            {
                aV[i] += MI[i][j]*KV[j];
                logL += 0.5*MI[i][j]*KV[i]*KV[j];
            }
        }
    
    x = (aV[0]+aV[3]);
    x = x*x;
    y = (aV[1]-aV[2]);
    y = y*y;
    u = (aV[0]-aV[3]);
    u = u*u;
    v = (aV[1]+aV[2]);
    v = v*v;
        
        Ap = sqrt(x+y)+sqrt(u+v);
        Ac = sqrt(x+y)-sqrt(u+v);
        A = Ap + sqrt(Ap*Ap-Ac*Ac);
        
        x = atan2((aV[1]-aV[2]),(aV[0]+aV[3]));
        
        y = atan2(-(aV[1]+aV[2]),(aV[0]-aV[3]));
        
        psi = 0.25*(y-x);
        while(psi < 0.0) psi += PI;
        while(psi > PI) psi -= PI;
        
        
        phic = 0.25*(x+y);
        while(phic < 0.0) phic += PI;
        while(phic > PI) phic -= PI;
        
        cosi = Ac/A;
        
        scale = 2.0/A;
        

      if(ll == 0)
      {
        params[6] *= scale;
      }
      else
      {
       params[6] += log(scale);
      }
    
       params[10] = cosi;
       params[9] = psi;
       params[4] = phic;
    
    // catch failed maximization
    if(params[6] != params[6] || params[4] != params[4] || params[9] != params[9] || params[10] != params[10]) logL = -1.0e20;

    /*   Deallocate Arrays   */
    
     free(KV);
     free(aV);
     free_double_matrix(MM,4);
     free_double_matrix(MI,4);
     free_double_matrix(filtA,4);
     free_double_matrix(filtE,4);
    
    
    return logL;
}

double likelihoodFstatOld(int ll, double *params, int N, double *AC, double *EC, double *SA, double *SE, double Tobs, double Tzero)
{
    
    /*   Indicies   */
    int i,j, k, n, m, imax;
    
    /*   Miscellaneous  */
    double xm, fstep, power, om, mx, tx;
    
    double *pfilt;
    
    double **filtA, **filtE;
    
    double **MM, **MI;
    
    double *KV, *aV;
    
    double AR, AI, ER, EI;
    
    double t, f, kdotx, Amp, Phase, Fp, Fc;
    
    double A, P, px, tc, tcs, pxi, pxj;
    
    double x, y, u, v, xold, yold;
    
    double cp, sp;
    
    int NA, flag;
    
    double cx, sx, logL, logLmax;
    
    int nn;
    
    double iota, cosi;
    double Ap, Ac, scale;
    
    double a1, a2, a3, a4;
    double cpsi, spsi;
    double psi, phic;
    double cpc, spc;
    
    FILE *out;
    
    pfilt = (double*)malloc(sizeof(double)* (NP));
    
    filtA = double_matrix(4,N);
    filtE = double_matrix(4,N);
    
    KV = (double*)malloc(sizeof(double)* (4));
    aV = (double*)malloc(sizeof(double)* (4));
    MM = double_matrix(4,4);
    MI = double_matrix(4,4);
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
    // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
    
    
    for (i=0; i< NP; i++) pfilt[i] = params[i];

    pfilt[10] = 0.0;
    
    pfilt[4] = 0.0;
    pfilt[9] = 0.0;
    ResponseFreq(ll, pfilt, N, filtA[0], filtE[0], Tobs, Tzero);
    
    pfilt[4] = PI/2.0;
    pfilt[9] = PI/4.0;
    ResponseFreq(ll, pfilt, N, filtA[1], filtE[1], Tobs, Tzero);
    
    pfilt[4] = 3.0*PI/4.0;
    pfilt[9] = 0.0;
    ResponseFreq(ll, pfilt, N, filtA[2], filtE[2], Tobs, Tzero);
    
    pfilt[4] = PI/4.0;
    pfilt[9] = PI/4.0;
    ResponseFreq(ll, pfilt, N, filtA[3], filtE[3], Tobs, Tzero);
     
    
        for (i=0; i< 4; i++)
        {
            KV[i] = 2.0*(fourier_nwip(filtA[i], AC, SA, N) + fourier_nwip(filtE[i], EC, SE, N));
            for (j=i; j< 4; j++)
            {
                MM[i][j] =  4.0*(fourier_nwip(filtA[i], filtA[j], SA, N) + fourier_nwip(filtE[i], filtE[j], SE, N));
               
            }
        }
        
        for (i=0; i< 4; i++)
        {
            for (j=i; j< 4; j++)
            {
                MM[j][i] = MM[i][j];
            }
        }
        
        // a small stabilizer
        for (i=0; i< 4; i++) MM[i][i] += 0.1;
        
        Inverse(MM, MI, 4);
        
        logL = 0.0;
        
        for (i=0; i< 4; i++)
        {
            aV[i] = 0.0;
            for (j=0; j< 4; j++)
            {
                aV[i] += MI[i][j]*KV[j];
                logL += 0.5*MI[i][j]*KV[i]*KV[j];
            }
        }
    
    x = (aV[0]+aV[3]);
    x = x*x;
    y = (aV[1]-aV[2]);
    y = y*y;
    u = (aV[0]-aV[3]);
    u = u*u;
    v = (aV[1]+aV[2]);
    v = v*v;
        
        Ap = sqrt(x+y)+sqrt(u+v);
        Ac = sqrt(x+y)-sqrt(u+v);
        A = Ap + sqrt(Ap*Ap-Ac*Ac);
        
        x = atan2((aV[1]-aV[2]),(aV[0]+aV[3]));
        
        y = atan2(-(aV[1]+aV[2]),(aV[0]-aV[3]));
        
        psi = 0.25*(y-x);
        while(psi < 0.0) psi += PI;
        while(psi > PI) psi -= PI;
        
        
        phic = 0.25*(x+y);
        while(phic < 0.0) phic += PI;
        while(phic > PI) phic -= PI;
        
        cosi = Ac/A;
        
        scale = 2.0/A;
        

      if(ll == 0)
      {
        params[6] *= scale;
      }
      else
      {
       params[6] += log(scale);
      }
    
       params[10] = cosi;
       params[9] = psi;
       params[4] = phic;
    
    // catch failed maximization
    if(params[6] != params[6] || params[4] != params[4] || params[9] != params[9] || params[10] != params[10]) logL = -1.0e20;

    /*   Deallocate Arrays   */
    
     free(KV);
     free(aV);
     free_double_matrix(MM,4);
     free_double_matrix(MI,4);
     free_double_matrix(filtA,4);
     free_double_matrix(filtE,4);
    
    
    return logL;
}

// Fstat likelihood that allows for a small (roughly +/- 30 s) time maximization
// Not working very well...
double likelihoodFstatTmax(int ll, double *params, int N, double *AC, double *EC, double *SA, double *SE, double Tobs, double Tzero)
{
    
    /*   Indicies   */
    int i,j, k, n, m, imax;
    
    /*   Miscellaneous  */
    double xm, fstep, power, om, mx, tx;
    
    double *pfilt;
    
    double **filtA, **filtE;
    
    double **MM, **MI;
    
    double *KV, *aV;
    
    double *KP, *KM;
    
    double AR, AI, ER, EI;
    
    double t, f, kdotx, Amp, Phase, Fp, Fc;
    
    double A, P, px, tc, tcs, pxi, pxj;
    
    double x, y, u, v, xold, yold;
    
    double cp, sp;
    
    int NA, flag;
    
    double cx, sx, logL, logLmax;
    
    double delT = 50.0;
    
    int nn;
    
    double iota, cosi;
    double Ap, Ac, scale;
    
    double a1, a2, a3, a4;
    double cpsi, spsi;
    double psi, phic;
    double cpc, spc;
    
    double T0, T1, T2, a, b, c, dtx;
    
    FILE *out;
    
    pfilt = (double*)malloc(sizeof(double)* (NP));
    
    filtA = double_matrix(4,N);
    filtE = double_matrix(4,N);
    
    KV = (double*)malloc(sizeof(double)* (4));
    KP = (double*)malloc(sizeof(double)* (4));
    KM = (double*)malloc(sizeof(double)* (4));
    
    aV = (double*)malloc(sizeof(double)* (4));
    MM = double_matrix(4,4);
    MI = double_matrix(4,4);
    
   
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
       // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
    

    ResponseFstat(ll, params, N, filtA, filtE, Tobs, Tzero);
    
   //printf("\n");
        for (i=0; i< 4; i++)
        {
            KV[i] = 2.0*(fourier_nwip(filtA[i], AC, SA, N) + fourier_nwip(filtE[i], EC, SE, N));
            KP[i] = 2.0*(fourier_nwip_Tshift(delT, Tobs, filtA[i], AC, SA, N) + fourier_nwip_Tshift(delT, Tobs, filtE[i], EC, SE, N));
            KM[i] = 2.0*(fourier_nwip_Tshift(-delT, Tobs, filtA[i], AC, SA, N) + fourier_nwip_Tshift(-delT, Tobs, filtE[i], EC, SE, N));
            for (j=i; j< 4; j++)
            {
                MM[i][j] =  4.0*(fourier_nwip(filtA[i], filtA[j], SA, N) + fourier_nwip(filtE[i], filtE[j], SE, N));
                //printf("%e ", MM[i][j]);
            }
            //printf("\n");
        }
        
        for (i=0; i< 4; i++)
        {
            for (j=i; j< 4; j++)
            {
                MM[j][i] = MM[i][j];
            }
        }
        
        // a small stabilizer
        for (i=0; i< 4; i++) MM[i][i] += 0.1;
        
        //x = det(MM,4);
        
        //printf("%e\n", x);
        
       // printf("%e %e %e %e\n", MM[0][0], MM[1][1], MM[2][2], MM[1][2]);
       //printf("%e %e %e %e\n", NV[0], NV[1], NV[2], NV[3]);
        
        Inverse(MM, MI, 4);
    
    
    T0 = 0.0;
    T1 = 0.0;
    T2 = 0.0;
    
    for (i=0; i< 4; i++)
    {
        for (j=0; j< 4; j++)
        {
            T0 += 0.5*MI[i][j]*KV[i]*KV[j];
            T1 += 0.5*MI[i][j]*KP[i]*KP[j];
            T2 += 0.5*MI[i][j]*KM[i]*KM[j];
        }
    }
    
    c = T0;
    b = (T1-T2)/(2.0*delT);
    a = (T1-b*delT-c)/(delT*delT);
    
    dtx = -b/(2.0*a);
                        
    //printf("%f\n", dtx);
    
    // This method is not well define for shifts of more than a few delT
    if(dtx < -2.0*delT) dtx = -2.0*delT;
    if(dtx > 2.0*delT) dtx = 2.0*delT;
    
    params[5] += dtx;
    
    for (i=0; i< 4; i++)
    {
        KV[i] = 2.0*(fourier_nwip_Tshift(dtx, Tobs, filtA[i], AC, SA, N) + fourier_nwip_Tshift(dtx, Tobs, filtE[i], EC, SE, N));
    }
        
        logL = 0.0;
        
        for (i=0; i< 4; i++)
        {
            aV[i] = 0.0;
            for (j=0; j< 4; j++)
            {
                aV[i] += MI[i][j]*KV[j];
                logL += 0.5*MI[i][j]*KV[i]*KV[j];
            }
        }
        
       //printf("%f %f %f %f\n", logL, T0, T1, T2);
        
        //printf("%e %e %e %e\n", aV[0], aV[1], aV[2], aV[3]);
    
    x = (aV[0]+aV[3]);
    x = x*x;
    y = (aV[1]-aV[2]);
    y = y*y;
    u = (aV[0]-aV[3]);
    u = u*u;
    v = (aV[1]+aV[2]);
    v = v*v;
        
        Ap = sqrt(x+y)+sqrt(u+v);
        Ac = sqrt(x+y)-sqrt(u+v);
        A = Ap + sqrt(Ap*Ap-Ac*Ac);
        
        x = atan2((aV[1]-aV[2]),(aV[0]+aV[3]));
        
        y = atan2(-(aV[1]+aV[2]),(aV[0]-aV[3]));
        
        psi = 0.25*(y-x);
        while(psi < 0.0) psi += PI;
        while(psi > PI) psi -= PI;
        
        
        phic = 0.25*(x+y);
        while(phic < 0.0) phic += PI;
        while(phic > PI) phic -= PI;
        
        cosi = Ac/A;
        
        scale = 2.0/A;
        
      // printf("%f %f %f %f %f %f %f %f %f\n", scale, x, y, phic, psi, cosi, params[4], params[9], params[10]);

      if(ll == 0)
      {
        params[6] *= scale;
      }
      else
      {
       params[6] += log(scale);
      }
    
       params[10] = cosi;
       params[9] = psi;
       params[4] = phic;
    
    // catch failed maximization
    if(params[6] != params[6] || params[4] != params[4] || params[9] != params[9] || params[10] != params[10]) logL = -1.0e20;

    /*   Deallocate Arrays   */
    
     free(KV);
     free(KP);
     free(KM);
     free(aV);
     free_double_matrix(MM,4);
     free_double_matrix(MI,4);
     free_double_matrix(filtA,4);
     free_double_matrix(filtE,4);
    
    
    return logL;
}

double Likelihood(int ll, double *params, long N, double *AD, double *ED, double *SA, double *SE, double Tobs, double Tzero)
{
    double *AS, *ES;
    double HH, HD;
    double logL;
    
    
    AS = (double*)malloc(sizeof(double)* (N));
    ES = (double*)malloc(sizeof(double)* (N));
    
    ResponseFreq(ll, params, N, AS, ES, Tobs, Tzero);
    
    HH = (fourier_nwip(AS, AS, SA, N)+fourier_nwip(ES, ES, SE, N));
    HD = (fourier_nwip(AS, AD, SA, N)+fourier_nwip(ES, ED, SE, N));
    
    logL = HD-0.5*HH;
    
    /*
     printf("%.14e %.14e %.14e\n", logL, HH, HD);
    
    FILE *out;
    int i;
    double AA, EE;
    out = fopen("int_full_search.dat","w");
          for (i=1; i< N/2; i++)
          {
              AA = 4.0*(AS[i]*AS[i]+AS[N-i]*AS[N-i])/SA[i];
              EE = 4.0*(ES[i]*ES[i]+ES[N-i]*ES[N-i])/SE[i];
              fprintf(out,"%e %e %e\n", (double)(i)/Tobs, AA, EE);
          }
          fclose(out);
     */
    
    free(AS);
    free(ES);
    
    return(logL);
    
}




void ResponseFreq(int ll, double *params, long N, double *AS, double *ES, double Tobs, double Tzero)
{
    
    /*   Indicies   */
    int i,j, k, ii, n, m, a, M, nn, nmin, nmax;
    
    /*   GW Source data   */
    double Mc, theta, phi, psi, D, iota, A, Aplus, Across, f0, fdot, phio;
    double costh, sinth, cosph, sinph, cosi, cosps, sinps;
    
    /*   Time and distance variables   */
    double xi, t, Tend, told, tx;
    
    /*   Miscellaneous  */
    double xm, fstep, power, om, mx, x;
    
    double Amp, Phase, fonfs, f;
    
    double Aprime, Pprime, fi;
    
    double HC, HS, hp, hc;
    
    double m1, m2, chi1, chi2, phic, tc, distance, dm, Mtot, eta, fr, af;
    
    double *ta, *xia, *FF;
    
    double Fp, Fc, kdotx, delt, fmin, fmax, dt;
    
    double XR, XI, YR, YI, ZR, ZI;
    
    int nfmin, nfmax, nf, flag;
    
    int NA;
    
    double fstart, fstop, fny;
    
    double m1_SI, m2_SI, deltaF;
    
     double fx, sqrtTobs;
    
 
    double *FpAR, *FpAI, *FcAR, *FcAI;
    double *FpER, *FpEI, *FcER, *FcEI;
    
    FILE *out;
    
    dt = Tobs/(double)(N);
    fny = 1.0/(2.0*dt);
    Tend = Tzero+Tobs;
    sqrtTobs = sqrt(Tobs);
    
    if(ll == 0)  // linear in m1, m2, DL
    {
    m1 = params[0]*TSUN;    // Mass1
    m2 = params[1]*TSUN;    // Mass2
    distance = params[6]*1.0e9*PC_SI; // distance
    Mtot = (m1+m2);
    eta = m1*m2/((m1+m2)*(m1+m2));
    Mc = pow(m1*m2,3.0/5.0)/pow(m1+m2,1.0/5.0);
    }
    else if(ll == 1)  // log in m1, m2, DL
    {
     m1 = exp(params[0])*TSUN;    // Mass1
     m2 = exp(params[1])*TSUN;    // Mass2
     distance = exp(params[6])*1.0e9*PC_SI; // distance
     Mtot = (m1+m2);
     eta = m1*m2/((m1+m2)*(m1+m2));
     Mc = pow(m1*m2,3.0/5.0)/pow(m1+m2,1.0/5.0);
    }
    else // log in Mc, Mt, DL
    {
    distance = exp(params[6])*1.0e9*PC_SI; // distance
    Mc = exp(params[0])*TSUN;
    Mtot = exp(params[1])*TSUN;
    eta = pow((Mc/Mtot), (5.0/3.0));
     if(eta > 0.25)
     {
        dm = 0.0;
     }
     else
     {
        dm = sqrt(1.0-4.0*eta);
     }
    m1 = Mtot*(1.0+dm)/2.0;
    m2 = Mtot*(1.0-dm)/2.0;
    }
    
    //printf("%e %e %e\n", m1, m2, Mtot);
    
    m1_SI = m1*MSUN_SI/TSUN;
    m2_SI = m2*MSUN_SI/TSUN;
    
    chi1 = params[2];  // Spin1
    chi2 = params[3];  // Spin2
    phic = params[4];  // merger phase
    tc = params[5];    // merger time
    cosi = params[10];  // inclination
    
    Aplus = 0.5*(1.+cosi*cosi);
    Across = -cosi;
    
    AmpPhaseFDWaveform *ap = NULL;
    double fRef_in=PDfref;
    double *AF, *TF;
    int ret, flag1, flag2;
    
    AF = (double*)malloc(sizeof(double)* (N/2));
    TF = (double*)malloc(sizeof(double)* (N/2));
    
    StartStop(ll, params, Tobs, Tzero, Tend, &fstart, &fstop, &fr);
    
    if(fstop > fny) fstop = fny;
    
    //printf("%.15e %e %e\n", fstart, fstop, fr);
    
    nfmin = (int)(fstart*Tobs);
    if(nfmin < 1) nfmin = 1;
    nfmax = (int)(fstop*Tobs);
    if(nfmax <= nfmin) nfmax = nfmin+1;
    if(nfmax > N/2) nfmax = N/2;
    nf = nfmax-nfmin;
    
    if(nf < 1)
    {
        nf = 1;
        if(nfmax < 2)
        {
            nfmin = 1;
            nfmax = 2;
        }
        if(nfmax == N/2)
        {
            nfmin = N/2-1;
        }
    }
    
    fmin = (double)(nfmin)/Tobs;
    fmax = (double)(nfmax)/Tobs;
    
    deltaF = 1.0/Tobs;
    
    //printf("%e %e %d %ld\n", fmin, fmax, nfmin, nf);
    
    
    ta = (double*)malloc(sizeof(double)* (nf));
    xia = (double*)malloc(sizeof(double)* (nf));
    FpAR = (double*)malloc(sizeof(double)* (nf));
    FcAR = (double*)malloc(sizeof(double)* (nf));
    FpER = (double*)malloc(sizeof(double)* (nf));
    FcER = (double*)malloc(sizeof(double)* (nf));
    FpAI = (double*)malloc(sizeof(double)* (nf));
    FcAI = (double*)malloc(sizeof(double)* (nf));
    FpEI = (double*)malloc(sizeof(double)* (nf));
    FcEI = (double*)malloc(sizeof(double)* (nf));
    FF = (double*)malloc(sizeof(double)* (nf));
    
    RealVector *freq;
    freq = CreateRealVector(nf);
    for (i=0; i< nf; i++) freq->data[i] = fmin+(double)(i)*deltaF;
    for (i=0; i< nf; i++) FF[i] = freq->data[i];
    
    // start = clock();
    
    if(distance != distance) printf("nan %e\n", params[6]);
    if(distance <= 0.0) printf("<= 0 %e\n", params[6]);

    
    ret = IMRPhenomDGenerateh22FDAmpPhase(
                                          &ap,      /**< [out] FD waveform */
                                          freq, /**< Input: frequencies (Hz) on which to evaluate h22 FD - will be copied in the output AmpPhaseFDWaveform. Frequencies exceeding max freq covered by PhenomD will be given 0 amplitude and phase. */
                                          0.0,                  /**< Orbital phase at fRef (rad) */
                                          fRef_in,               /**< reference frequency (Hz) */
                                          m1_SI,                 /**< Mass of companion 1 (kg) */
                                          m2_SI,                 /**< Mass of companion 2 (kg) */
                                          chi1,                  /**< Aligned-spin parameter of companion 1 */
                                          chi2,                  /**< Aligned-spin parameter of companion 2 */
                                          distance               /**< Distance of source (m) */
                                          );
    
    //end = clock();
    //cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    //printf("PhenomD call took %f seconds\n", cpu_time_used);
    
    // start = clock();
    
    
    // fill the Time - Frequency series
    
    flag = 0;
    told = ap->time[0]+tc;
    
    for (i=0; i< nf; i++)
    {
        t = ap->time[i]+tc;
        
        if(t < told && flag == 0)
        {
            flag = 1;
            ii = i-1;
            tx = told;
            fx = freq->data[ii];
        }
        
        TF[i] = t;
        if(flag == 1) TF[i] = tx + (freq->data[i]-fx)/fr*Mtot;
        told = t;
    }
    
    RAantenna(params, nf, TF, FF, xia, FpAR, FpAI, FcAR, FcAI, FpER, FpEI, FcER, FcEI);
    
    /*
     out = fopen("check.dat","w");
     for (i=0; i< nf; i++)
     {
     fprintf(out,"%d %e %e %e\n", i, freq->data[i], TF[i], ap->amp[i]);
     }
     fclose(out);
     */
    
    
    //end = clock();
    //cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    // printf("time array call took %f seconds\n", cpu_time_used);
    
    // start = clock();
    
    for (i=0; i< nf; i++)
    {
        AF[i] =  h22fac*ap->amp[i]/sqrtTobs;
        fonfs = freq->data[i]/fstar;
        AF[i] *= (4.0*fonfs*sin(fonfs));   // conversion to fractional frequency and leading order TDI transfer function
    }
    
    
    for(n=0; n< N; n++)
    {
        AS[n] = 0.0;
        ES[n] = 0.0;
    }
    
    /*   Main Loop   */
    
    // start = clock();
    
    //printf("%d %d\n", nfmin, nfmax);
    
    nn = 0;
    for(n=nfmin; n< nfmax; n++)
    {
        // Barycenter time and frequency
        
        // The frequency and time arrays start at nfmin
        m = n-nfmin;
        
        if(m > -1 && m < nf)
        {
            t = TF[m];
            f = FF[m];
            xi = xia[m];
            
            kdotx = t-xi;
            
            // Tukey filter to match what is done in time domain
            x = 1.0;
            if(t < Tzero+t_tuke) x = 0.5*(1.0+cos(PI*((t-Tzero)/t_tuke-1.0)));
            if(t > Tend-t_tuke && t < Tend)
            {
            x = 0.5*(1.0-cos(PI*(t-Tend)/t_tuke));
            }
            if(t > Tend) x = 0.0;
            
            Amp = x*AF[m];
            Phase = ap->phase[m]+2.0*phic;
            
            HC = Amp*cos(2.0*PI*f*(Tend-tc+dt/2.0-kdotx)-Phase);
            HS = Amp*sin(2.0*PI*f*(Tend-tc+dt/2.0-kdotx)-Phase);
            

            AS[n] = FpAR[m]*Aplus*HC - FpAI[m]*Aplus*HS - FcAR[m]*Across*HS - FcAI[m]*Across*HC;
            AS[N-n] = FpAI[m]*Aplus*HC + FpAR[m]*Aplus*HS - FcAI[m]*Across*HS + FcAR[m]*Across*HC;
            
            ES[n] = FpER[m]*Aplus*HC - FpEI[m]*Aplus*HS - FcER[m]*Across*HS - FcEI[m]*Across*HC;
            ES[N-n] = FpEI[m]*Aplus*HC + FpER[m]*Aplus*HS - FcEI[m]*Across*HS + FcER[m]*Across*HC;
           
            
        }
        
        
    }
    
    /*   Deallocate Arrays   */
    
    free(AF);
    free(TF);
    
    free(ta);
    free(xia);
    free(FcAR);
    free(FpAR);
    free(FcER);
    free(FpER);
    free(FcAI);
    free(FpAI);
    free(FcEI);
    free(FpEI);
    free(FF);
    
    DestroyAmpPhaseFDWaveform(ap);
    DestroyRealVector(freq);
    
    return;
}

void ResponseFstat(int ll, double *params, long N, double **AS, double **ES, double Tobs, double Tzero)
{
    
    /*   Indicies   */
    int i,j, k, ii, n, m, a, M, nn, nmin, nmax;
    
    /*   GW Source data   */
    double Mc, theta, phi, psi, D, iota, A, Aplus, Across, f0, fdot, phio;
    double costh, sinth, cosph, sinph, cosi, cosps, sinps;
    
    /*   Time and distance variables   */
    double xi, t, Tend, told, tx;
    
    /*   Miscellaneous  */
    double xm, fstep, power, om, mx, x;
    
    double Amp, Phase, fonfs, f;
    
    double Aprime, Pprime, fi;
    
    double HC, HS, hp, hc;
    
    double m1, m2, chi1, chi2, phic, tc, distance, dm, Mtot, eta, fr, af;
    
    double *ta, *xia, *FF;
    
    double Fp, Fc, kdotx, delt, fmin, fmax, dt;
    
    double XR, XI, YR, YI, ZR, ZI;
    
    int nfmin, nfmax, nf, flag;
    
    int NA;
    
    double fstart, fstop, fny;
    
    double m1_SI, m2_SI, deltaF;
    
     double fx, sqrtTobs;
    
    double *pfstat;
    
 
    double *FpAR, *FpAI, *FcAR, *FcAI;
    double *FpER, *FpEI, *FcER, *FcEI;
    
    FILE *out;
    
    dt = Tobs/(double)(N);
    fny = 1.0/(2.0*dt);
    Tend = Tzero+Tobs;
    sqrtTobs = sqrt(Tobs);
    
    if(ll == 0)  // linear in m1, m2, DL
    {
    m1 = params[0]*TSUN;    // Mass1
    m2 = params[1]*TSUN;    // Mass2
    distance = params[6]*1.0e9*PC_SI; // distance
    Mtot = (m1+m2);
    eta = m1*m2/((m1+m2)*(m1+m2));
    Mc = pow(m1*m2,3.0/5.0)/pow(m1+m2,1.0/5.0);
    }
    else if(ll == 1)  // log in m1, m2, DL
    {
     m1 = exp(params[0])*TSUN;    // Mass1
     m2 = exp(params[1])*TSUN;    // Mass2
     distance = exp(params[6])*1.0e9*PC_SI; // distance
     Mtot = (m1+m2);
     eta = m1*m2/((m1+m2)*(m1+m2));
     Mc = pow(m1*m2,3.0/5.0)/pow(m1+m2,1.0/5.0);
    }
    else // log in Mc, Mt, DL
    {
    distance = exp(params[6])*1.0e9*PC_SI; // distance
    Mc = exp(params[0])*TSUN;
    Mtot = exp(params[1])*TSUN;
    eta = pow((Mc/Mtot), (5.0/3.0));
     if(eta > 0.25)
     {
        dm = 0.0;
     }
     else
     {
        dm = sqrt(1.0-4.0*eta);
     }
    m1 = Mtot*(1.0+dm)/2.0;
    m2 = Mtot*(1.0-dm)/2.0;
    }
    
    //printf("%e %e %e\n", m1, m2, Mtot);
    
    m1_SI = m1*MSUN_SI/TSUN;
    m2_SI = m2*MSUN_SI/TSUN;
    
    chi1 = params[2];  // Spin1
    chi2 = params[3];  // Spin2
    tc = params[5];    // merger time

    AmpPhaseFDWaveform *ap = NULL;
    double fRef_in=PDfref;
    double *AF, *TF;
    int ret, flag1, flag2;
    
    AF = (double*)malloc(sizeof(double)* (N/2));
    TF = (double*)malloc(sizeof(double)* (N/2));
    
    StartStop(ll, params, Tobs, Tzero, Tend, &fstart, &fstop, &fr);
    
    if(fstop > fny) fstop = fny;
    
    //printf("%.15e %e %e\n", fstart, fstop, fr);
    
    nfmin = (int)(fstart*Tobs);
    if(nfmin < 1) nfmin = 1;
    nfmax = (int)(fstop*Tobs);
    if(nfmax <= nfmin) nfmax = nfmin+1;
    if(nfmax > N/2) nfmax = N/2;
    nf = nfmax-nfmin;
    
    if(nf < 1)
    {
        nf = 1;
        if(nfmax < 2)
        {
            nfmin = 1;
            nfmax = 2;
        }
        if(nfmax == N/2)
        {
            nfmin = N/2-1;
        }
    }
    
    fmin = (double)(nfmin)/Tobs;
    fmax = (double)(nfmax)/Tobs;
    
    deltaF = 1.0/Tobs;
    
    //printf("%e %e %d %ld\n", fmin, fmax, nfmin, nf);
    
    
    ta = (double*)malloc(sizeof(double)* (nf));
    xia = (double*)malloc(sizeof(double)* (nf));
    FpAR = (double*)malloc(sizeof(double)* (nf));
    FcAR = (double*)malloc(sizeof(double)* (nf));
    FpER = (double*)malloc(sizeof(double)* (nf));
    FcER = (double*)malloc(sizeof(double)* (nf));
    FpAI = (double*)malloc(sizeof(double)* (nf));
    FcAI = (double*)malloc(sizeof(double)* (nf));
    FpEI = (double*)malloc(sizeof(double)* (nf));
    FcEI = (double*)malloc(sizeof(double)* (nf));
    FF = (double*)malloc(sizeof(double)* (nf));
    
    RealVector *freq;
    freq = CreateRealVector(nf);
    for (i=0; i< nf; i++) freq->data[i] = fmin+(double)(i)*deltaF;
    for (i=0; i< nf; i++) FF[i] = freq->data[i];
    
    if(distance != distance) printf("distance nan %e\n", params[6]);
    if(distance <= 0.0) printf("<= 0 %e\n", params[6]);
    
    ret = IMRPhenomDGenerateh22FDAmpPhase(
                                          &ap,      /**< [out] FD waveform */
                                          freq, /**< Input: frequencies (Hz) on which to evaluate h22 FD - will be copied in the output AmpPhaseFDWaveform. Frequencies exceeding max freq covered by PhenomD will be given 0 amplitude and phase. */
                                          0.0,                  /**< Orbital phase at fRef (rad) */
                                          fRef_in,               /**< reference frequency (Hz) */
                                          m1_SI,                 /**< Mass of companion 1 (kg) */
                                          m2_SI,                 /**< Mass of companion 2 (kg) */
                                          chi1,                  /**< Aligned-spin parameter of companion 1 */
                                          chi2,                  /**< Aligned-spin parameter of companion 2 */
                                          distance               /**< Distance of source (m) */
                                          );
    
    //end = clock();
    //cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    //printf("PhenomD call took %f seconds\n", cpu_time_used);
    
    // start = clock();
    
    
    // fill the Time - Frequency series
    
    flag = 0;
    told = ap->time[0]+tc;
    
    for (i=0; i< nf; i++)
    {
        t = ap->time[i]+tc;
        
        if(t < told && flag == 0)
        {
            flag = 1;
            ii = i-1;
            tx = told;
            fx = freq->data[ii];
        }
        
        TF[i] = t;
        if(flag == 1) TF[i] = tx + (freq->data[i]-fx)/fr*Mtot;
        told = t;
    }
    

    //   phic = params[4];  // merger phase
    //   cosi = params[10];  // inclination
    
    RAFstat(params, nf, TF, FF, xia, FpAR, FpAI, FcAR, FcAI, FpER, FpEI, FcER, FcEI);
    
    Aplus = 0.5;
    Across = 0.0;
    
    for (i=0; i< nf; i++)
    {
        AF[i] =  h22fac*ap->amp[i]/sqrtTobs;
        fonfs = freq->data[i]/fstar;
        AF[i] *= (4.0*fonfs*sin(fonfs));   // conversion to fractional frequency and leading order TDI transfer function
    }
    
     for (i=0; i< 4; i++)
      {
      for(n=0; n< N; n++)
      {
        AS[i][n] = 0.0;
        ES[i][n] = 0.0;
      }
     }
    
    /*   Main Loop   */
    
    
    nn = 0;
    for(n=nfmin; n< nfmax; n++)
    {
        // Barycenter time and frequency
        
        // The frequency and time arrays start at nfmin
        m = n-nfmin;
        
        if(m > -1 && m < nf)
        {
            t = TF[m];
            f = FF[m];
            xi = xia[m];
            
            kdotx = t-xi;
            
            // Tukey filter to match what is done in time domain
            x = 1.0;
            if(t < Tzero+t_tuke) x = 0.5*(1.0+cos(PI*((t-Tzero)/t_tuke-1.0)));
            if(t > Tend-t_tuke && t < Tend)
            {
            x = 0.5*(1.0-cos(PI*(t-Tend)/t_tuke));
            }
            if(t > Tend) x = 0.0;
            
            Amp = x*AF[m];
            Phase = ap->phase[m];
        
            HC = Amp*cos(2.0*PI*f*(Tend-tc+dt/2.0-kdotx)-Phase);
            HS = Amp*sin(2.0*PI*f*(Tend-tc+dt/2.0-kdotx)-Phase);
            
           // The cosps = 0, sinps = 1 case is the same as the cosps = 1, sinps = 0 case with fp(01) = fc(10) amd fc(01) = -fp(10)
            
                   // 1 0  HC -> HC, HS -> HS
                 
                   AS[0][n] = FpAR[m]*Aplus*HC - FpAI[m]*Aplus*HS;
                   AS[0][N-n] = FpAI[m]*Aplus*HC + FpAR[m]*Aplus*HS;
             
                   ES[0][n] = FpER[m]*Aplus*HC - FpEI[m]*Aplus*HS;
                   ES[0][N-n] = FpEI[m]*Aplus*HC + FpER[m]*Aplus*HS;
                 
                 // 0 1    HC -> -HC, HS -> -HS
                 
                   AS[1][n] = -FcAR[m]*Aplus*HC + FcAI[m]*Aplus*HS;
                   AS[1][N-n] = -FcAI[m]*Aplus*HC - FcAR[m]*Aplus*HS;
             
                   ES[1][n] = -FcER[m]*Aplus*HC + FcEI[m]*Aplus*HS;
                   ES[1][N-n] = -FcEI[m]*Aplus*HC - FcER[m]*Aplus*HS;
             

                 // 1 0   HC -> -HS, HS -> HC
            
                    AS[2][n] = -FpAR[m]*Aplus*HS - FpAI[m]*Aplus*HC;
                    AS[2][N-n] = -FpAI[m]*Aplus*HS + FpAR[m]*Aplus*HC;
                         
                    ES[2][n] = -FpER[m]*Aplus*HS - FpEI[m]*Aplus*HC;
                    ES[2][N-n] = -FpEI[m]*Aplus*HS + FpER[m]*Aplus*HC;
             
                 // 0 1    HC -> HS, HS -> -HC
            
                     AS[3][n] = FcAR[m]*Aplus*HS + FcAI[m]*Aplus*HC;
                     AS[3][N-n] = FcAI[m]*Aplus*HS - FcAR[m]*Aplus*HC;
                        
                     ES[3][n] = FcER[m]*Aplus*HS + FcEI[m]*Aplus*HC;
                     ES[3][N-n] = FcEI[m]*Aplus*HS - FcER[m]*Aplus*HC;
            
            
        }
        
        
    }
    
    /*   Deallocate Arrays   */
    
    free(AF);
    free(TF);
    
    free(ta);
    free(xia);
    free(FcAR);
    free(FpAR);
    free(FcER);
    free(FpER);
    free(FcAI);
    free(FpAI);
    free(FcEI);
    free(FpEI);
    free(FF);
    
    DestroyAmpPhaseFDWaveform(ap);
    DestroyRealVector(freq);
    
    return;
}



void RAantenna(double *params, int NF, double *TF, double *FF, double *xi, double *FpAR, double *FpAI, double *FcAR, double *FcAI,
               double *FpER, double *FpEI, double *FcER, double *FcEI)
{
    
    /*   Indicies   */
    int i,j, k, n, m, a, M;
    
    /*   Gravitational Wave basis vectors   */
    double *u,*v,*kv;
    
    /*   Polarization basis tensors   */
    double **eplus, **ecross;
    
    /*   Spacecraft position and separation vector   */
    double *x, *y, *z;
    double *r12, *r13, *r21, *r23, *r31, *r32;
    double *r10, *r20, *r30;
    double *vx, *vy, *vz;
    
    double q1, q2, q3, q4;
    
    /*   Dot products   */
    double kdotx;
    
    /*   Convenient quantities   */
    double **dplus, **dcross;
    
    /*   GW Source data   */
    double Mc, theta, phi, psi, D, iota, A, Aplus, Across, f0, fdot, phio;
    double costh, sinth, cosph, sinph, cosi, cosps, sinps;
    
    /*   Time and distance variables   */
    double t, xa, ya, za;
    
    /*   Miscellaneous  */
    double xm, fstep, power, om, mx;
    
    double delt;
    
    double fpx, fcx, fpy, fcy, fpz, fcz;
    
    double **TR, **TI, **kdr;
    
    double *kdg;
    
    double fr;
    
    /*   Allocating Arrays   */
    
    u = double_vector(4); v = double_vector(4); kv = double_vector(4);
    
    eplus  = double_matrix(4,4); ecross = double_matrix(4,4);
    
    dplus  = double_matrix(4,4); dcross = double_matrix(4,4);
    
    TR  = double_matrix(4,4); TI = double_matrix(4,4);
    
    kdr = double_matrix(4,4);
    
    kdg = double_vector(4);
    
    x = double_vector(4); y = double_vector(4); z = double_vector(4);
    
    r12 = double_vector(4); r21 = double_vector(4); r31 = double_vector(4);
    r13 = double_vector(4); r23 = double_vector(4); r32 = double_vector(4);
    r10 = double_vector(4); r20 = double_vector(4); r30 = double_vector(4);
    
    phi = params[8];   // EclipticLongitude
    psi = params[9];   // polarization
    //Calculate cos and sin of sky position, inclination, polarization
    costh = params[7];   sinth = sqrt(1.0-costh*costh);
    cosph = cos(phi);     sinph = sin(phi);
    cosps = cos(2.*psi);  sinps = sin(2.*psi);
    
    
    /*   Tensor basis  */
    v[1] =  -costh*cosph;
    v[2] =  -costh*sinph;
    v[3] = sinth;
    
    u[1] =  sinph;
    u[2] = -cosph;
    u[3] =  0.;
    
    
    kv[1] = -sinth*cosph;
    kv[2] = -sinth*sinph;
    kv[3] = -costh;
    
    
    
    for(i=1;i<=3;i++)
    {
        for(j=1;j<=3;j++)
        {
            eplus[i][j]  = u[i]*u[j] - v[i]*v[j];
            ecross[i][j] = u[i]*v[j] + v[i]*u[j];
        }
    }
    
    /*   Main Loop   */
    for(n=0; n< NF; n++)
    {
        // Barycenter time
        t = TF[n];
        fr = FF[n]/(2.0*fstar);
        
        spacecraft(t, x, y, z);
        
        // guiding center
        xa = (x[1]+x[2]+x[3])/3.0;
        ya = (y[1]+y[2]+y[3])/3.0;
        za = (z[1]+z[2]+z[3])/3.0;
        
        kdotx = (xa*kv[1]+ya*kv[2]+za*kv[3])/clight;
        
        // detector time and frequency
        xi[n]  = t - kdotx;
        
        //Unit separation vector from spacecraft i to j
        r12[1] = (x[2] - x[1])/Larm;   r13[1] = (x[3] - x[1])/Larm;   r23[1] = (x[3] - x[2])/Larm;
        r12[2] = (y[2] - y[1])/Larm;   r13[2] = (y[3] - y[1])/Larm;   r23[2] = (y[3] - y[2])/Larm;
        r12[3] = (z[2] - z[1])/Larm;   r13[3] = (z[3] - z[1])/Larm;   r23[3] = (z[3] - z[2])/Larm;
        
        // These are not unit vectors. Just pulling out the Larm scaling
        r10[1] = (xa-x[1])/Larm;   r10[2] = (ya-y[1])/Larm;  r10[3] = (za-z[1])/Larm;
        r20[1] = (xa-x[2])/Larm;   r20[2] = (ya-y[2])/Larm;  r20[3] = (za-z[2])/Larm;
        r30[1] = (xa-x[3])/Larm;   r30[2] = (ya-y[3])/Larm;  r30[3] = (za-z[3])/Larm;
        
        kdr[1][2] = 0.0;
        for(k=1; k<=3; k++) kdr[1][2] += kv[k]*r12[k];
        kdr[1][3] = 0.0;
        for(k=1; k<=3; k++) kdr[1][3] += kv[k]*r13[k];
        kdr[2][3] = 0.0;
        for(k=1; k<=3; k++) kdr[2][3] += kv[k]*r23[k];
        
        kdr[2][1] = -kdr[1][2];  kdr[3][1] = -kdr[1][3];  kdr[3][2] = -kdr[2][3];
        
        kdg[1] = 0.0;
        for(k=1; k<=3; k++) kdg[1] += kv[k]*r10[k];
        kdg[2] = 0.0;
        for(k=1; k<=3; k++) kdg[2] += kv[k]*r20[k];
        kdg[3] = 0.0;
        for(k=1; k<=3; k++) kdg[3] += kv[k]*r30[k];
        
        //Make use of symmetry
        for(i=1; i<=3; i++)
        {
            r21[i] = -r12[i];
            r31[i] = -r13[i];
            r32[i] = -r23[i];
        }
        
        
        for(i=1; i<=3; i++)
        {
            for(j=1; j<=3; j++)
            {
                q1 = fr*(1.0-kdr[i][j]);
                q2 = fr*(1.0+kdr[i][j]);
                q3 = -fr*(3.0+kdr[i][j]-2.0*kdg[i]);
                q4 = -fr*(1.0+kdr[i][j]-2.0*kdg[i]);
                q1 = (sin(q1)/q1);
                q2 = (sin(q2)/q2);
                TR[i][j] = 0.5*(q1*cos(q3)+q2*cos(q4));   // goes to 1 when f/fstat small
                TI[i][j] = 0.5*(q1*sin(q3)+q2*sin(q4));   // goes to 0 when f/fstat small
            }
        }
        

        
        dplus[1][2] = dplus[1][3] = dplus[2][1] = dplus[2][3] = dplus[3][1] = dplus[3][2] = 0.;
        dcross[1][2] = dcross[1][3] = dcross[2][1] = dcross[2][3] = dcross[3][1] = dcross[3][2] = 0.;
        //Convenient quantities d+ & dx
        for(i=1; i<=3; i++)
        {
            for(j=1; j<=3; j++)
            {
                dplus[1][2]  += r12[i]*r12[j]*eplus[i][j];   dcross[1][2] += r12[i]*r12[j]*ecross[i][j];
                dplus[2][3]  += r23[i]*r23[j]*eplus[i][j];   dcross[2][3] += r23[i]*r23[j]*ecross[i][j];
                dplus[1][3]  += r13[i]*r13[j]*eplus[i][j];   dcross[1][3] += r13[i]*r13[j]*ecross[i][j];
            }
        }
        
        dplus[2][1] = dplus[1][2];  dcross[2][1] = dcross[1][2];
        dplus[3][2] = dplus[2][3];  dcross[3][2] = dcross[2][3];
        dplus[3][1] = dplus[1][3];  dcross[3][1] = dcross[1][3];
        
        fpx = -0.5*( (dplus[1][2]*cosps+dcross[1][2]*sinps)*TR[1][2] - (dplus[1][3]*cosps+dcross[1][3]*sinps)*TR[1][3] );
        fcx = -0.5*( (-dplus[1][2]*sinps+dcross[1][2]*cosps)*TR[1][2] - (-dplus[1][3]*sinps + dcross[1][3]*cosps)*TR[1][3] );
        
        fpy = -0.5*( (dplus[2][3]*cosps+dcross[2][3]*sinps)*TR[2][3] - (dplus[2][1]*cosps+dcross[2][1]*sinps)*TR[2][1] );
        fcy = -0.5*( (-dplus[2][3]*sinps+dcross[2][3]*cosps)*TR[2][3] - (-dplus[2][1]*sinps + dcross[2][1]*cosps)*TR[2][1] );
        
        fpz = -0.5*( (dplus[3][1]*cosps+dcross[3][1]*sinps)*TR[3][1] - (dplus[3][2]*cosps+dcross[3][2]*sinps)*TR[3][2] );
        fcz = -0.5*( (-dplus[3][1]*sinps+dcross[3][1]*cosps)*TR[3][1] - (-dplus[3][2]*sinps + dcross[3][2]*cosps)*TR[3][2] );

        /* Original AET definition
        FpAR[n] = (2.0*fpx-fpy-fpz)/3.0;
        FcAR[n] = (2.0*fcx-fcy-fcz)/3.0;
        
        FpER[n] = (fpz-fpy)/sq3;
        FcER[n] = (fcz-fcy)/sq3;*/
        
        
        /* LDC's AET definition */
        FpAR[n] = (fpz-fpx)*0.707106781186547;
        FcAR[n] = (fcz-fcx)*0.707106781186547;
        
        FpER[n] = (fpx-2*fpy+fpz)*0.408248290463863;
        FcER[n] = (fcx-2*fcy+fcz)*0.408248290463863;
 
        fpx = -0.5*( (dplus[1][2]*cosps+dcross[1][2]*sinps)*TI[1][2] - (dplus[1][3]*cosps+dcross[1][3]*sinps)*TI[1][3] );
        fcx = -0.5*( (-dplus[1][2]*sinps+dcross[1][2]*cosps)*TI[1][2] - (-dplus[1][3]*sinps + dcross[1][3]*cosps)*TI[1][3] );
                                         
        fpy = -0.5*( (dplus[2][3]*cosps+dcross[2][3]*sinps)*TI[2][3] - (dplus[2][1]*cosps+dcross[2][1]*sinps)*TI[2][1] );
        fcy = -0.5*( (-dplus[2][3]*sinps+dcross[2][3]*cosps)*TI[2][3] - (-dplus[2][1]*sinps + dcross[2][1]*cosps)*TI[2][1] );
                                                               
        fpz = -0.5*( (dplus[3][1]*cosps+dcross[3][1]*sinps)*TI[3][1] - (dplus[3][2]*cosps+dcross[3][2]*sinps)*TI[3][2] );
        fcz = -0.5*( (-dplus[3][1]*sinps+dcross[3][1]*cosps)*TI[3][1] - (-dplus[3][2]*sinps + dcross[3][2]*cosps)*TI[3][2] );
            
        /* Original AET definition
        FpAI[n] = (2.0*fpx-fpy-fpz)/3.0;
        FcAI[n] = (2.0*fcx-fcy-fcz)/3.0;
                                                                                     
        FpEI[n] = (fpz-fpy)/sq3;
        FcEI[n] = (fcz-fcy)/sq3;*/
        
        
        /* LDC's AET definition */
        FpAI[n] = (fpz-fpx)*0.707106781186547;
        FcAI[n] = (fcz-fcx)*0.707106781186547;
        
        FpEI[n] = (fpx-2*fpy+fpz)*0.408248290463863;
        FcEI[n] = (fcx-2*fcy+fcz)*0.408248290463863;

        
    }
    
    free_double_vector(u); free_double_vector(v); free_double_vector(kv);
    
    free_double_matrix(eplus,4); free_double_matrix(ecross,4);
    
    free_double_matrix(dplus,4); free_double_matrix(dcross,4);
    
    free_double_matrix(TR,4); free_double_matrix(TI,4);
    
    free_double_matrix(kdr,4);
    
    free_double_vector(kdg);
    
    free_double_vector(x); free_double_vector(y); free_double_vector(z);
    
    free_double_vector(r12); free_double_vector(r21); free_double_vector(r31);
    free_double_vector(r13); free_double_vector(r23); free_double_vector(r32);
    free_double_vector(r10); free_double_vector(r20); free_double_vector(r30);
    
    return;
}

void RAFstat(double *params, int NF, double *TF, double *FF, double *xi, double *FpAR, double *FpAI, double *FcAR, double *FcAI, double *FpER, double *FpEI, double *FcER, double *FcEI)
{
    
    /*   Indicies   */
    int i,j, k, n, m, a, M;
    
    /*   Gravitational Wave basis vectors   */
    double *u,*v,*kv;
    
    /*   Polarization basis tensors   */
    double **eplus, **ecross;
    
    /*   Spacecraft position and separation vector   */
    double *x, *y, *z;
    double *r12, *r13, *r21, *r23, *r31, *r32;
    double *r10, *r20, *r30;
    double *vx, *vy, *vz;
    
    double q1, q2, q3, q4;
    
    /*   Dot products   */
    double kdotx;
    
    /*   Convenient quantities   */
    double **dplus, **dcross;
    
    /*   GW Source data   */
    double Mc, theta, phi, psi, D, iota, A, Aplus, Across, f0, fdot, phio;
    double costh, sinth, cosph, sinph, cosi, cosps, sinps;
    
    /*   Time and distance variables   */
    double t, xa, ya, za;
    
    /*   Miscellaneous  */
    double xm, fstep, power, om, mx;
    
    double delt;
    
    double fpx, fcx, fpy, fcy, fpz, fcz;
    
    double **TR, **TI, **kdr;
    
    double *kdg;
    
    double fr;
    
    /*   Allocating Arrays   */
    
    u = double_vector(4); v = double_vector(4); kv = double_vector(4);
    
    eplus  = double_matrix(4,4); ecross = double_matrix(4,4);
    
    dplus  = double_matrix(4,4); dcross = double_matrix(4,4);
    
    TR  = double_matrix(4,4); TI = double_matrix(4,4);
    
    kdr = double_matrix(4,4);
    
    kdg = double_vector(4);
    
    x = double_vector(4); y = double_vector(4); z = double_vector(4);
    
    r12 = double_vector(4); r21 = double_vector(4); r31 = double_vector(4);
    r13 = double_vector(4); r23 = double_vector(4); r32 = double_vector(4);
    r10 = double_vector(4); r20 = double_vector(4); r30 = double_vector(4);
    
    phi = params[8];   // EclipticLongitude
    psi = params[9];   // polarization
    //Calculate cos and sin of sky position, inclination, polarization
    costh = params[7];   sinth = sqrt(1.0-costh*costh);
    cosph = cos(phi);     sinph = sin(phi);
    
    
    /*   Tensor basis  */
    v[1] =  -costh*cosph;
    v[2] =  -costh*sinph;
    v[3] = sinth;
    
    u[1] =  sinph;
    u[2] = -cosph;
    u[3] =  0.;
    
    
    kv[1] = -sinth*cosph;
    kv[2] = -sinth*sinph;
    kv[3] = -costh;
    
    
    
    for(i=1;i<=3;i++)
    {
        for(j=1;j<=3;j++)
        {
            eplus[i][j]  = u[i]*u[j] - v[i]*v[j];
            ecross[i][j] = u[i]*v[j] + v[i]*u[j];
        }
    }
    
    /*   Main Loop   */
    for(n=0; n< NF; n++)
    {
        // Barycenter time
        t = TF[n];
        fr = FF[n]/(2.0*fstar);
        
        spacecraft(t, x, y, z);
        
        // guiding center
        xa = (x[1]+x[2]+x[3])/3.0;
        ya = (y[1]+y[2]+y[3])/3.0;
        za = (z[1]+z[2]+z[3])/3.0;
        
        kdotx = (xa*kv[1]+ya*kv[2]+za*kv[3])/clight;
        
        // detector time and frequency
        xi[n]  = t - kdotx;
        
        //Unit separation vector from spacecraft i to j
        r12[1] = (x[2] - x[1])/Larm;   r13[1] = (x[3] - x[1])/Larm;   r23[1] = (x[3] - x[2])/Larm;
        r12[2] = (y[2] - y[1])/Larm;   r13[2] = (y[3] - y[1])/Larm;   r23[2] = (y[3] - y[2])/Larm;
        r12[3] = (z[2] - z[1])/Larm;   r13[3] = (z[3] - z[1])/Larm;   r23[3] = (z[3] - z[2])/Larm;
        
        // These are not unit vectors. Just pulling out the Larm scaling
        r10[1] = (xa-x[1])/Larm;   r10[2] = (ya-y[1])/Larm;  r10[3] = (za-z[1])/Larm;
        r20[1] = (xa-x[2])/Larm;   r20[2] = (ya-y[2])/Larm;  r20[3] = (za-z[2])/Larm;
        r30[1] = (xa-x[3])/Larm;   r30[2] = (ya-y[3])/Larm;  r30[3] = (za-z[3])/Larm;
        
        kdr[1][2] = 0.0;
        for(k=1; k<=3; k++) kdr[1][2] += kv[k]*r12[k];
        kdr[1][3] = 0.0;
        for(k=1; k<=3; k++) kdr[1][3] += kv[k]*r13[k];
        kdr[2][3] = 0.0;
        for(k=1; k<=3; k++) kdr[2][3] += kv[k]*r23[k];
        
        kdr[2][1] = -kdr[1][2];  kdr[3][1] = -kdr[1][3];  kdr[3][2] = -kdr[2][3];
        
        kdg[1] = 0.0;
        for(k=1; k<=3; k++) kdg[1] += kv[k]*r10[k];
        kdg[2] = 0.0;
        for(k=1; k<=3; k++) kdg[2] += kv[k]*r20[k];
        kdg[3] = 0.0;
        for(k=1; k<=3; k++) kdg[3] += kv[k]*r30[k];
        
        //Make use of symmetry
        for(i=1; i<=3; i++)
        {
            r21[i] = -r12[i];
            r31[i] = -r13[i];
            r32[i] = -r23[i];
        }
        
        
        for(i=1; i<=3; i++)
        {
            for(j=1; j<=3; j++)
            {
                q1 = fr*(1.0-kdr[i][j]);
                q2 = fr*(1.0+kdr[i][j]);
                q3 = -fr*(3.0+kdr[i][j]-2.0*kdg[i]);
                q4 = -fr*(1.0+kdr[i][j]-2.0*kdg[i]);
                q1 = (sin(q1)/q1);
                q2 = (sin(q2)/q2);
                TR[i][j] = 0.5*(q1*cos(q3)+q2*cos(q4));   // goes to 1 when f/fstat small
                TI[i][j] = 0.5*(q1*sin(q3)+q2*sin(q4));   // goes to 0 when f/fstat small
            }
        }
        

        
        dplus[1][2] = dplus[1][3] = dplus[2][1] = dplus[2][3] = dplus[3][1] = dplus[3][2] = 0.;
        dcross[1][2] = dcross[1][3] = dcross[2][1] = dcross[2][3] = dcross[3][1] = dcross[3][2] = 0.;
        //Convenient quantities d+ & dx
        for(i=1; i<=3; i++)
        {
            for(j=1; j<=3; j++)
            {
                dplus[1][2]  += r12[i]*r12[j]*eplus[i][j];   dcross[1][2] += r12[i]*r12[j]*ecross[i][j];
                dplus[2][3]  += r23[i]*r23[j]*eplus[i][j];   dcross[2][3] += r23[i]*r23[j]*ecross[i][j];
                dplus[1][3]  += r13[i]*r13[j]*eplus[i][j];   dcross[1][3] += r13[i]*r13[j]*ecross[i][j];
            }
        }
        
        dplus[2][1] = dplus[1][2];  dcross[2][1] = dcross[1][2];
        dplus[3][2] = dplus[2][3];  dcross[3][2] = dcross[2][3];
        dplus[3][1] = dplus[1][3];  dcross[3][1] = dcross[1][3];
        
        //cosps = 1, sinps = 0
        
        // The cosps = 0, sinps = 1 case is the same as the cosps = 1, sinps = 0 case with fp(01) = fc(10) amd fc(01) = -fp(10)
        
        fpx = -0.5*( dplus[1][2]*TR[1][2] - dplus[1][3]*TR[1][3] );
        fcx = -0.5*( dcross[1][2]*TR[1][2] - dcross[1][3]*TR[1][3] );
        
        fpy = -0.5*(dplus[2][3]*TR[2][3] - dplus[2][1]*TR[2][1] );
        fcy = -0.5*( dcross[2][3]*TR[2][3] - dcross[2][1]*TR[2][1] );
        
        fpz = -0.5*(dplus[3][1]*TR[3][1] - dplus[3][2]*TR[3][2] );
        fcz = -0.5*( dcross[3][1]*TR[3][1] - dcross[3][2]*TR[3][2] );
        
        FpAR[n] = (2.0*fpx-fpy-fpz)/3.0;
        FcAR[n] = (2.0*fcx-fcy-fcz)/3.0;
               
        FpER[n] = (fpz-fpy)/sq3;
        FcER[n] = (fcz-fcy)/sq3;
        
        fpx = -0.5*( dplus[1][2]*TI[1][2] - dplus[1][3]*TI[1][3] );
        fcx = -0.5*( dcross[1][2]*TI[1][2] - dcross[1][3]*TI[1][3] );
        
        fpy = -0.5*(dplus[2][3]*TI[2][3] - dplus[2][1]*TI[2][1] );
        fcy = -0.5*( dcross[2][3]*TI[2][3] - dcross[2][1]*TI[2][1] );
        
        fpz = -0.5*(dplus[3][1]*TI[3][1] - dplus[3][2]*TI[3][2] );
        fcz = -0.5*( dcross[3][1]*TI[3][1] - dcross[3][2]*TI[3][2] );
                                                                                     
        FpAI[n] = (2.0*fpx-fpy-fpz)/3.0;
        FcAI[n] = (2.0*fcx-fcy-fcz)/3.0;
                                                                                     
        FpEI[n] = (fpz-fpy)/sq3;
        FcEI[n] = (fcz-fcy)/sq3;
        
    }
    
    free_double_vector(u); free_double_vector(v); free_double_vector(kv);
    
    free_double_matrix(eplus,4); free_double_matrix(ecross,4);
    
    free_double_matrix(dplus,4); free_double_matrix(dcross,4);
    
    free_double_matrix(TR,4); free_double_matrix(TI,4);
    
    free_double_matrix(kdr,4);
    
    free_double_vector(kdg);
    
    free_double_vector(x); free_double_vector(y); free_double_vector(z);
    
    free_double_vector(r12); free_double_vector(r21); free_double_vector(r31);
    free_double_vector(r13); free_double_vector(r23); free_double_vector(r32);
    free_double_vector(r10); free_double_vector(r20); free_double_vector(r30);
    
    return;
}

void spacecraft(double t, double *x, double *y, double *z)
{

  double alpha;
  double beta1, beta2, beta3;
  double sa, sb, ca, cb;
 
  alpha = 2.*pi*fm*t + kappa0;

  beta1 = 0. + lambda0;
  beta2 = 2.*pi/3. + lambda0;
  beta3 = 4.*pi/3. + lambda0;

  sa = sin(alpha);
  ca = cos(alpha);


  sb = sin(beta1);
  cb = cos(beta1);
  x[1] = AU*ca + AU*ec*(sa*ca*sb - (1. + sa*sa)*cb);
  y[1] = AU*sa + AU*ec*(sa*ca*cb - (1. + ca*ca)*sb);
  z[1] = -sq3*AU*ec*(ca*cb + sa*sb);

 
  sb = sin(beta2);
  cb = cos(beta2);
  x[2] = AU*ca + AU*ec*(sa*ca*sb - (1. + sa*sa)*cb);
  y[2] = AU*sa + AU*ec*(sa*ca*cb - (1. + ca*ca)*sb);
  z[2] = -sq3*AU*ec*(ca*cb + sa*sb);

  sb = sin(beta3);
  cb = cos(beta3);
  x[3] = AU*ca + AU*ec*(sa*ca*sb - (1. + sa*sa)*cb);
  y[3] = AU*sa + AU*ec*(sa*ca*cb - (1. + ca*ca)*sb);
  z[3] = -sq3*AU*ec*(ca*cb + sa*sb);
  
}

void quadratures(double *ct, double *st, double *data1, double *data2, double *Sn, int n)
{
    int nb2, i, l, k, j;
    int imax, imin;
    double cx, sx;
    
    nb2 = n/2;
    
    cx = 0.0;
    sx = 0.0;
    
    for (i=1; i < nb2; i++)
    {
        l=i;
        k=n-i;
        cx += (data1[l]*data2[l] + data1[k]*data2[k])/Sn[i];
        sx += (data1[k]*data2[l] - data1[l]*data2[k])/Sn[i];
    }
    
    *ct = 4.0*cx;
    *st = 4.0*sx;
    
    
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
        
        corr[l]	= (data1[l]*data2[l] + data1[k]*data2[k])/Sn[i];
        corr[k]	= (data1[k]*data2[l] - data1[l]*data2[k])/Sn[i];
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


void getfreq(double *fnew, double *tf, double t, double fguess, double m1_SI, double m2_SI, double chi1, double chi2, double tc)
{
    AmpPhaseFDWaveform *ap = NULL;
    double ep, u, v, tnew, x;
    double delT, delF, dtdf, fonfs;
    int ret;
    double M_sec;
    RealVector *freq;
    
    M_sec = (m1_SI+m2_SI) * MTSUN_SI/MSUN_SI;
    
    ep = 1.0e-6/M_sec;
    v = (4.0*PI*ep);
    u = (2.0*PI*ep*ep);
    
    freq = CreateRealVector((3));
    
    freq->data[0] = fguess-ep;
    freq->data[1] = fguess;
    freq->data[2] = fguess+ep;
    
    if(freq->data[0] < 0.0)
    {
        freq->data[0] = fguess;
        freq->data[1] = fguess+ep;
        freq->data[2] = fguess+2.0*ep;
    }
    
    ret = IMRPhenomDGenerateh22FDAmpPhase(&ap, freq, 0.0, PDfref, m1_SI, m2_SI, chi1, chi2, PC_SI);
    
    tnew = (ap->phase[2]-ap->phase[0])/v +tc;
    
    dtdf = (ap->phase[2]+ap->phase[0]-2.0*ap->phase[1])/u;
    
    delT = t-tnew;
    
    delF = delT/dtdf;
    
    *fnew = fguess + delF;
    
    *tf = tnew;
    
    DestroyAmpPhaseFDWaveform(ap);
    DestroyRealVector(freq);
    
}



void bwlf(double *in, double *out, int fwrv, int M, int n, double s, double f)
{
    /* Butterworth bandpass filter
     n = filter order 2,4,6,8,...
     s = sampling frequency
     f = half power frequency
     */
    
    if(n % 2){ printf("Order must be 2,4,6,8,...\n"); return;}
    
    int i, j;
    n = n/2;
    double a = tan(PI*f/s);
    double a2 = a*a;
    double r;
    double *A = (double *)malloc(n*sizeof(double));
    double *d1 = (double *)malloc(n*sizeof(double));
    double *d2 = (double *)malloc(n*sizeof(double));
    double *w0 = (double *)calloc(n, sizeof(double));
    double *w1 = (double *)calloc(n, sizeof(double));
    double *w2 = (double *)calloc(n, sizeof(double));
    double x;
    
    for(i=0; i<n; ++i)
    {
        r = sin(PI*(2.0*i+1.0)/(4.0*n));
        s = a2 + 2.0*a*r + 1.0;
        A[i] = a2/s;
        d1[i] = 2.0*(1-a2)/s;
        d2[i] = -(a2 - 2.0*a*r + 1.0)/s;
        w0[i] = 0.0;
        w1[i] = 0.0;
        w2[i] = 0.0;
    }
    
    for(j=0; j< M; ++j)
    {
        if(fwrv == 1) x = in[j];
        if(fwrv == -1) x = in[M-j-1];
        for(i=0; i<n; ++i)
        {
            w0[i] = d1[i]*w1[i] + d2[i]*w2[i] + x;
            x = A[i]*(w0[i] + 2.0*w1[i] + w2[i]);
            w2[i] = w1[i];
            w1[i] = w0[i];
        }
        if(fwrv == 1) out[j] = x;
        if(fwrv == -1) out[M-j-1] = x;
    }
    
    free(A);
    free(d1);
    free(d2);
    free(w0);
    free(w1);
    free(w2);
    
    return;
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

double *double_vector(int N)
{
    return malloc( (N+1) * sizeof(double) );
}

void free_double_vector(double *v)
{
    free(v);
}

int *int_vector(int N)
{
    return malloc( (N+1) * sizeof(int) );
}

void free_int_vector(int *v)
{
    free(v);
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

double ***double_tensor(int N, int M, int P)
{
    int i,j;
    
    double ***t = malloc( (N+1) * sizeof(double **));
    for(i=0; i<N+1; i++)
    {
        t[i] = malloc( (M+1) * sizeof(double *));
        for(j=0; j<M+1; j++)
        {
            t[i][j] = malloc( (P+1) * sizeof(double));
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


