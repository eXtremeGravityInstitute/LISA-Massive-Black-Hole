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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

#include <time.h>

#include "Declarations.h"

#ifndef _OPENMP
#define omp ignore
#endif

// OSX
// clang -Xpreprocessor -fopenmp -lomp -w -o  PTMCMC PTMCMC.c Utils.c Response.c IMRPhenomD_internals.c IMRPhenomD.c -lgsl -lgslcblas  -lm

// Linux
// gcc -std=gnu99 -fopenmp -w -o PTMCMC PTMCMC.c IMRPhenomD_internals.c IMRPhenomD.c Utils.c Response.c -lgsl -lgslcblas  -lm




//##############################################
//MT modifications

gsl_rng **rvec;
//##############################################

int main(int argc,char **argv)
{

  double f, fdot, theta, phi, A, iota, psi, phase;
  char Gfile[50];
  char filename[50];
  double *params, *premove, *pnew;
  double *AS, *ES;
  double *AQ, *EQ;
  double AR, AI, ER, EI;
  double fonfs, Sn, Sm, Acut;
  double Aar, Aai;
  double x, y, z;
  long M, q;
  long i, j, k, cnt, mult, id, NS;
  double Mc, fstart, fstop, fr;
  double SNR, logL;
  double HH, HD, HDQ, DD, Match, MX, ts, ps;
  double HHA, HDA, DDA, HHE, HDE, DDE, MA, ME;
    
  struct Het *het   = malloc(sizeof(struct Het));
  struct Data *dat  = malloc(sizeof(struct Data));
    
  double *AC, *EC, *TC;
    
    int seg, rep;
    
    
    clock_t start, end;
    double cpu_time_used;
    
    const gsl_rng_type * P;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    
    P = gsl_rng_default;
    r = gsl_rng_alloc (P);
    
  FILE* in;
  FILE* out;
    
    //##############################################
    //open MP modifications
    omp_set_num_threads(NC);
    rvec = (gsl_rng **)malloc(sizeof(gsl_rng *) * (NC+1));
    for(i = 0 ; i<= NC; i++){
        rvec[i] = gsl_rng_alloc(P);
        gsl_rng_set(rvec[i] , i);
    }
    //##############################################
    
    
    params = (double*)malloc(sizeof(double)* (NP));
    premove = (double*)malloc(sizeof(double)* (NP));
    
    
    if(LDC == 1)  // LDC = 1 tells the code to analyze LDC data. Otherwise it does a noise-free run on the parameters provided
    {
    
    // setting the segment # to -1 causes the code to use the full data set
    
    if(argc<3)
    {
        printf("./PTMCMC segment# source#\n");
        printf("segment numbers run from 0 to 12\n");
        printf("source numbers start at 0\n");
        return 0;
    }
    
    seg = atoi(argv[1]);
    rep = atoi(argv[2]);
    
    if(seg > 0)
    {
        dat->Tobs = Tsegment;
    }
    else
    {
        dat->Tobs = 16.0*Tsegment;
    }
    
    dat->sqrtTobs = sqrt(dat->Tobs);
        
    dat->dt = cadence;
    
    dat->Nch = 2;  // only analyze A, E
    dat->N = (int)(dat->Tobs/dat->dt);
    
    dat->SN = double_matrix(dat->Nch,dat->N/2);
    dat->SM = double_matrix(dat->Nch,dat->N/2);
    dat->data = double_matrix(dat->Nch,dat->N);
    
    // count the sources
    in = fopen("search_sources.dat","r");
    NS = -1;
    while ((j = fgetc(in)) != EOF)
    {
     fscanf(in,"%lf%lf", &x, &x);
     for(i=0; i< NP; i++) fscanf(in,"%lf", &x);
    NS++;
    }
    
    rewind(in);

    // extract the source to be analyzed
    for(k=0; k < NS; k++)
    {
        if(k == rep)
        {
            fscanf(in,"%lf%lf", &x, &x);
            for(i=0; i< NP; i++) fscanf(in,"%lf", &params[i]);
            //printf("%e\n", params[5]);
        }
        else // wind the file
        {
        fscanf(in,"%lf%lf", &x, &x);
        for(i=0; i< NP; i++) fscanf(in,"%lf", &x);
        }
    }
    fclose(in);
    
    // read in the data and PSDs
    if(seg > -1)
    {
    
    dat->Tstart = (double)(seg)*dat->Tobs;
    dat->Tend = (double)(seg+1)*dat->Tobs;
    
    // read in the previously estimated smooth and full PSDs
    for(id=0; id < dat->Nch; id++)
    {
      sprintf(filename, "specfit_%d_%d.dat", id, seg);
      in = fopen(filename,"r");
      for(i=0; i< dat->N/2; i++)
      {
        fscanf(in,"%lf%lf%lf%lf\n", &f, &dat->SM[id][i], &x, &dat->SN[id][i]);
      }
      fclose(in);
    }
 
    // Read in FFTed LDC data
    sprintf(filename, "AET_seg%d_f.dat", seg);
    in = fopen(filename,"r");
    for(i=0; i< dat->N; i++)
    {
        fscanf(in,"%lf%lf%lf%lf\n", &f, &dat->data[0][i], &dat->data[1][i], &x);
    }
    fclose(in);
        
    }
    else  // using full data
    {
        dat->Tstart = 0.0;
        dat->Tend = dat->Tobs;
        
        // Read in FFTed LDC data
        in = fopen("AET_f.dat","r");
        for(i=0; i< dat->N; i++)
        {
            fscanf(in,"%lf%lf%lf%lf\n", &f, &dat->data[0][i], &dat->data[1][i], &x);
        }
        fclose(in);
        
        // read in the previously estimated smooth and full PSDs
        for(id=0; id < dat->Nch; id++)
        {
          sprintf(filename, "specav_%d.dat", id);
          in = fopen(filename,"r");
          for(i=0; i< dat->N/2; i++)
          {
            fscanf(in,"%lf%lf%lf\n", &f, &dat->SM[id][i], &dat->SN[id][i]);
          }
          fclose(in);
        }

    }
    
    printf("%.0f %.0f %d\n", dat->Tobs, dat->Tstart, dat->N);
    
    
    // form up residual by subtracting other BH signals
    AS = double_vector(dat->N);
    ES = double_vector(dat->N);
    in = fopen("search_sources.dat","r");
   
    for(k=0; k < NS; k++)
    {
        fscanf(in,"%lf%lf", &x, &x);
        for(i=0; i< NP; i++) fscanf(in,"%lf", &premove[i]);
        
        if(k == rep)
        {
            if(premove[5] < dat->Tstart || premove[5] > dat->Tend) printf("WARNING: source does not merge during the chosen time interval\n");
        }
        
        if(k != rep)
          {
          // only subtract sources that have not merged
            if(premove[5] > dat->Tstart)
            {
            map_params(2, premove);
            ResponseFreq(dat, 2, premove, AS, ES);
             for(i=0; i< dat->N; i++)
             {
              dat->data[0][i] -= AS[i];
              dat->data[1][i] -= ES[i];
             }
            }
          }
    }
    fclose(in);
    free_double_vector(AS);
    free_double_vector(ES);
        
        // change to better intrinsic parameterization
        map_params(2, params);
        
    }
    else // Do a noise-free run on the parameters provided (if looping over a catalog, have the parameters copied into source.dat
    {
        
        if(argc<3)
           {
               printf("./PTMCMC Tobs cadence\n");
               // typical cadence is around 5 seconds.
               return 0;
           }
          
           dat->Tobs = atof(argv[1]);
           dat->sqrtTobs = sqrt(dat->Tobs);
           dat->Tstart = 0.0;
           dat->Tend = dat->Tobs;
           dat->dt = atof(argv[2]);
           dat->Nch = 2;  // only analyze A, E
           dat->N = (int)(dat->Tobs/dat->dt);
           dat->SN = double_matrix(dat->Nch,dat->N/2);
           dat->SM = double_matrix(dat->Nch,dat->N/2);
           dat->data = double_matrix(dat->Nch,dat->N);
        
          in = fopen("source.dat","r");
          for(i=0; i< NP; i++) fscanf(in,"%lf", &params[i]);
          fclose(in);
        
          // change to better intrinsic parameterization
          map_params(2, params);
        
      
        
           for(i=1; i< dat->N/2; i++)
           {
               f = (double)(i)/dat->Tobs;
               instrument_noise(f, &dat->SM[0][i]);
               dat->SN[0][i] = dat->SM[0][i];
               dat->SM[1][i] = dat->SM[0][i];
               dat->SN[1][i] = dat->SM[0][i];
           }
        
            dat->SM[0][0] = dat->SM[0][1];
            dat->SN[0][0] = dat->SM[0][1];
            dat->SM[1][0] = dat->SM[0][1];
            dat->SN[1][0] = dat->SM[0][1];
        
        out = fopen("spec.dat","w");
        for(i=1; i< dat->N/2; i++)
         {
            f = (double)(i)/dat->Tobs;
            fprintf(out, "%e %e\n", f, dat->SM[0][i]);
         }
        fclose(out);
        
        AS = double_vector(dat->N);
        ES = double_vector(dat->N);
        
        // generate the reference signal (i.e. the "data")
        ResponseFreq(dat, 2, params, AS, ES);
        for(i=0; i< dat->N; i++)
        {
         dat->data[0][i] = AS[i];
         dat->data[1][i] = ES[i];
        }
          
    }

    
    double **paramx;
    paramx = double_matrix(NC,NP);
    
    // here we copy the best fit source into all the chains
    // in the global fit the previous state for each chain
    // will instead be passed in
    for (i=0; i< NC; i++)
    {
        for (j=0; j< NP; j++) paramx[i][j] = params[j];
    }
    
    // set up the whos-who reference. In the global fit this will be passed
    // over from the previous state of the BH model
    int *who;
    who = int_vector(NC);
    for (i=0; i< NC; i++) who[i] = i;
    
    // perform the MCMC
    MCMC(dat, het, 2, who, paramx);
    
    //###############################################
    //MT modification
   	for(i =0 ;i<= NC; i++){
        gsl_rng_free(rvec[i]);
    }
    free(rvec);
    //###############################################
    
    free(params);
 
    return 1;
    
}

void MCMC(struct Data *dat, struct Het *het, int ll, int *who, double **paramx)
{
    double *AAmp, *EAmp, *APhase, *EPhase;
    double kxm, *logLx, logLy, SNR, x, y, z, Lref, fc, fe;
    double *nhx;
    double lm1, lm2, chi1, chi2;
    double ctheta, phi;
    double m1, m2, DL;
    int i, j, k, km, n, q, mc, NF, N;
    int NFmax = 100000;
    double *AF, *PF, *FF, *TF, *SAE, *FS, *ASD;
    double *FN, *ASN;
    int NH=1000;
    int MP;
    int M = 2000000;
    double *heat;
    double **paramy;
    double *pnew;
    double *pref, *scale;
    double ***Fisher;
    double ***history;
    double *min, *max;
    double **ejump, ***evec, **diag;
    double **ejumpI, ***evecI;
    double alpha, beta, tau, f, df;
    double phic, tc, px, tm;
    double **sx, **sy;
    int *m;
    int *scount, *sacc, hold;
    int mcount, macc, typ;
    int **av, **cv;
    double Mc, Mtot, eta, dm;
    double a, b, c;
    
    double logLmax;
    double *pmax;
    
    double *ATR, *ATI, *ETR, *ETI;
    double *AS, *ES;
    double *ASR, *ESR;
    
    double **Cov, **Chl;
    double *zv;
    
    double ***iChl, ***eChl;
    
    double thetaL, phiL;
    
    clock_t start, end;
    double cpu_time_used;
    
    FILE *in;
    FILE *chain;
    FILE *out;
    FILE *levels;
    FILE *swaps;
    FILE *lcheck;
    FILE *lcheckh;
    
    const gsl_rng_type * P;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    
    P = gsl_rng_default;
    r = gsl_rng_alloc (P);
    
    MP = 6;
    
    N = dat->N;
    
    eChl = double_tensor(NC,NE,NE);
    iChl = double_tensor(NC,NI,NI);

    pmax = (double*)malloc(sizeof(double)* (NP));
    
    pnew = (double*)malloc(sizeof(double)* (NP));
    pref = (double*)malloc(sizeof(double)* (NP));
    scale = (double*)malloc(sizeof(double)* (NP));
    
    for (i=0; i< NP; i++) pref[i] = paramx[who[0]][i];

    FS = (double*)malloc(sizeof(double)* (NFmax));

    SetUp(dat, ll, pref, NFmax, &NF, FS);
    
    AF = (double*)malloc(sizeof(double)* (NF));
    PF = (double*)malloc(sizeof(double)* (NF));
    TF = (double*)malloc(sizeof(double)* (NF));
    
    Intrinsic(ll, pref, dat->Tstart, NF, FS, TF, PF, AF);
    
    AAmp = (double*)malloc(sizeof(double)* (NF));
    EAmp = (double*)malloc(sizeof(double)* (NF));
    APhase = (double*)malloc(sizeof(double)* (NF));
    EPhase = (double*)malloc(sizeof(double)* (NF));

    sacc = int_vector(NC);
    scount = int_vector(NC);
    heat = double_vector(NC);
    logLx = double_vector(NC);
    paramy = double_matrix(NC,NP);
    history = double_tensor(NC,NH,NP);
    Fisher = double_tensor(NC,NP,NP);
    ejump = double_matrix(NC,NP);
    evec = double_tensor(NC,NP,NP);
    
    // scale the PSDs
    sx = double_matrix(NC,dat->Nch);
    sy = double_matrix(NC,dat->Nch);

    Extrinsic(pref, dat->Tstart, dat->Tend, NF, FS, TF, PF, AF, AAmp, EAmp, APhase, EPhase, &kxm);
    
    FF = (double*)malloc(sizeof(double)* (NF+1));
    for(n=0; n< NF; n++) FF[n] = FS[n];
    
    free(FS);
    free(TF);
    free(PF);
    free(AF);
    
   
    // prior boundaries
    max = (double*)malloc(sizeof(double)* (NP));
    min = (double*)malloc(sizeof(double)* (NP));
    
    // prior boundaries
    max = (double*)malloc(sizeof(double)* (NP));
    min = (double*)malloc(sizeof(double)* (NP));
    
    // ll = 0  [0] Mass1  [1] Mass2
    // ll = 1  [0] ln Mass1  [1] ln Mass2
    // ll = 2  [0] ln Mc  [1] ln Mtot
    // [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln distance
    // [7] cos EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] cos inclination
    
    if(ll == 0)
    {
    max[0] = 5.0e8;
    max[1] = 5.0e8;
    min[0] = 1.0e3;
    min[1] = 1.0e3;
    }
    
    if(ll == 1)
    {
    max[0] = log(5.0e8);
    max[1] = log(5.0e8);
    min[0] = log(1.0e3);
    min[1] = log(1.0e3);
    }
    
    
    if(ll == 2)
    {
    max[0] = log(0.44*5.0e8);
    max[1] = log(5.0e8);
    min[0] = log(1.0e2);
    min[1] = log(1.0e3);
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

    
    // use one of the cold chains to produce the reference waveform
    het_space(dat, het, ll, paramx[who[0]], min, max);
    heterodyne(dat, het, ll, paramx[who[0]]);
    
    for (i=0; i< NC; i++)
       {
         for (j=0; j< dat->Nch; j++) sx[i][j] = 1.0;
       }
   

    for (i=0; i< NC; i++)
    {
        for (k=0; k< NH; k++)
        {
            for (j=0; j< NP; j++) history[i][k][j] = paramx[i][j];
        }
    }


    /*
     start = clock();
     for(i=0; i<100; i++) x = log_likelihood_het(dat, het, ll, paramx[0], sx[0]);
     end = clock();
     cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
     printf("heterodyned likelihood calculation took %f seconds\n", cpu_time_used/100.0);
    
      start = clock();
      for(i=0; i<100; i++) x = Likelihood(dat, ll, paramx[0]);
      end = clock();
      cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
      printf("fast likelihood calculation took %f seconds\n", cpu_time_used/100.0);
    
      start = clock();
      for(i=0; i<100; i++) x = Likelihood_Slow(dat, ll, paramx[0]);
      end = clock();
      cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
      printf("full likelihood calculation took %f seconds\n", cpu_time_used/100.0);
    */
    
    

    for (i=0; i< NC; i++)
    {
    logLx[i] = log_likelihood_het(dat, het, ll, paramx[i], sx[i]);
    }
    
    logLmax = logLx[who[0]];
    for (j=0; j< NP; j++) pmax[j] = paramx[who[0]][j];
     
     // run NCC cold chains
     for (i=0; i< NCC; i++) heat[i] = 1.0;
     
     
     SNR = het->SNR;
     x = pow((SNR/5.0),1.0/(double)(NC-NCC));
     if(x > 1.3) x = 1.3;
     for (i=NCC; i< NC; i++) heat[i] = heat[i-1]*x;
     printf("SNR %f increment %f max heat %f SNReff = %f\n", SNR, x, heat[NC-1], SNR/heat[NC-1]);

    #pragma omp parallel for
    for (i=0; i< NC; i++)
    {
    FisherHet(dat, het, ll, paramx[i], Fisher[i]);
    FisherEvec(Fisher[i], ejump[i], evec[i], NP);
    efix(dat, het, 1, ll, paramx[i], min, max, ejump[i], evec[i], 1.0);
    }

    av = int_matrix(5,NC);
    cv = int_matrix(5,NC);
    
    for(j = 0; j < 5; j++)
    {
        for(k=0; k < NC; k++)
        {
            av[j][k] = 0;
            cv[j][k] = 1;
        }
    }
    
    m = int_vector(NC);
    
    for(k=0; k < NC; k++)
    {
        m[k] = 0;
        sacc[k] = 0;
        scount[k] = 1;
    }
    
    mcount = 1;
    macc = 0;
    
    swaps = fopen("swaps.dat","w");
    levels = fopen("levels.dat","w");
    chain = fopen("chain.dat","w");
    lcheck = fopen("lcheck.dat","w");
    lcheckh = fopen("lcheck_hot.dat","w");
    
    if(nflag == 1)
    {
      a = 0.6;
      b = 0.3;
      c = 0.1;
    }
    else
    {
        a = 0.5;
        b = 0.2;
        c = 0.0;
    }
    
    for(mc = 1; mc <= M; mc++)
    {
        /*
        if(mc%1000==0)
        {
        freehet(het);
        het_space(dat, het, ll, pmax, min, max);
        heterodyne(dat, het, ll, pmax);
        }
        */
        
        if(mc%10000==0  && lhold == 0)
        {
            // update the Fisher matrices
            #pragma omp parallel for
            for(k = 0; k < NC; k++)
            {
                FisherHet(dat, het, ll, paramx[k], Fisher[k]);
                FisherEvec(Fisher[k], ejump[k], evec[k], NP);
                efix(dat, het, 1, ll, paramx[k], min, max, ejump[k], evec[k], 1.0);
            }
        }
        
        alpha = gsl_rng_uniform(r);
        
        if((NC > 1) && (alpha < 0.2))  // decide if we are doing a MCMC update of all the chains or a PT swap
        {
            
            // chain swap
           
            alpha = (double)(NC-1)*gsl_rng_uniform(r);
            j = (int)(alpha);
            beta = exp((logLx[who[j]]-logLx[who[j+1]])/heat[j+1] - (logLx[who[j]]-logLx[who[j+1]])/heat[j]);
            alpha = gsl_rng_uniform(r);
            if(beta > alpha)
            {
                hold = who[j];
                who[j] = who[j+1];
                who[j+1] = hold;
                sacc[j]++;
            }
            
             scount[j]++;
            
        }
        else      // MCMC update
        {
            
            mcount++;
           
            for(j = 0; j < NC; j++)
            {
                for(i = 0; i < NP; i++) paramy[j][i] = paramx[j][i];
                if(nflag == 1) for(i = 0; i < dat->Nch; i++) sy[j][i] = sx[j][i];
            }

        // all chains do the same type of update since some (especially type 2) are much slower than the others. Saves them waiting on others to finish
            alpha = gsl_rng_uniform(r);
            if(alpha > a)
            {
            typ = 0;
            }
           else if(alpha > b)
            {
            typ = 1;
            }
          else if(alpha > c)
            {
            typ = 2;
            }
           else
           {
             typ = 3;
           }
            
        #pragma omp parallel for
        for(k=0; k < NC; k++)
        {
          update(dat, het, typ, k, ll, logLx, paramx, paramy, sx, sy, min, max, who, heat, history, NH, ejump, evec, cv, av, rvec[k]);
         }
            
            for(k=0; k < NCC; k++)  // update maxL solution
            {
              q = who[k];
             if(logLx[q] > logLmax)
             {
              logLmax = logLx[q];
              for (j=0; j< NP; j++) pmax[j] = paramx[q][j];
             }
            }
            
            // add to the history file
            if(mc%10 == 0)
            {
                for(k=0; k < NC; k++)
                {
                    q = who[k];
                    i = m[k]%NH;
                    // the history file is kept for each temperature
                    for(j=0; j<NP; j++) history[k][i][j] = paramx[q][j];
                    m[k]++;
                }
            }
            
        }
        
            if(mc > 1000 && mc%20 == 0)
            {
                /*
                for(k=0; k < NC; k++)
                {
                    q = who[k];
                    x = log_likelihood_het(dat, het, ll, paramx[q]);
                    y = Likelihood_check(dat, het, ll, paramx[q]);
                    if(k < 4)
                    {
                        fprintf(lcheck,"%.12e %.12e\n", y, y-x);
                    }
                    else
                    {
                        fprintf(lcheckh,"%.12e %.12e\n", y, y-x);
                    }
                }*/
                
                for(k = 0; k < NCC; k++)
                {
                    q = who[k];
                    
                    if(ll == 0)
                    {
                    m1 = paramx[q][0];
                    m2 = paramx[q][1];
                    }
                    
                    if(ll == 1)
                    {
                    m1 = exp(paramx[q][0]);
                    m2 = exp(paramx[q][1]);
                    }
                    
                    if(ll == 2)
                    {
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
                    fprintf(chain,"%d %f %f %f %f\n", q, thetaL, phiL, sx[q][0], sx[q][1]);
                }
                
                
                fprintf(swaps, "%d ", mc);
                for(k = 0; k < NC-1; k++)
                {
                    fprintf(swaps, "%f ", (double)(sacc[k])/(double)(scount[k]));
                }
                fprintf(swaps, "\n");
                
                for(k = 0; k < NC; k++)
                {
                    fprintf(levels, "%.12e ", logLx[who[k]]);
                }
                fprintf(levels, "\n");
                
                
            }
            
            if(mc%100 == 0)
            {
                

                q = who[0];
                
                if(ll == 0)
                {
                m1 = paramx[q][0];
                m2 = paramx[q][1];
                }
                
                if(ll == 1)
                {
                m1 = exp(paramx[q][0]);
                m2 = exp(paramx[q][1]);
                }
                
                if(ll == 2)
                {
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
                
                
               if(nflag == 1)
               {
                    // only want to record the reduced log likelihood
                    x = logLx[q];
                   
                    for (i=0; i< dat->Nch; i++)
                    {
                      x += 0.5*het->DD[i]/sx[q][i];
                      x += (double)(het->MM-het->MN)*log(sx[q][i]);
                    }
                   
                   printf("%d %e %e %e %f %f %f %f %f %f %f\n", mc, x, m1, m2, sx[q][0], sx[q][1],
                       (double)(sacc[3])/(double)(scount[3]),
                       (double)(av[0][q])/(double)(cv[0][q]),
                       (double)(av[1][q])/(double)(cv[1][q]),
                       (double)(av[2][q])/(double)(cv[2][q]),
                       (double)(av[3][q])/(double)(cv[3][q]));
               }
                else
                {
                    printf("%d %e %e %e %f %f %f %f\n", mc, logLx[q], m1, m2,
                    (double)(sacc[3])/(double)(scount[3]),
                    (double)(av[0][q])/(double)(cv[0][q]),
                    (double)(av[1][q])/(double)(cv[1][q]),
                    (double)(av[2][q])/(double)(cv[2][q]));
                }
            
            }

            
        }
        
    fclose(chain);
    fclose(levels);
    fclose(swaps);
    fclose(lcheck);
    fclose(lcheckh);

    free(sx);
    free(sy);
    free_int_matrix(av,5);
    free_int_matrix(cv,5);
    free_int_vector(m);
    free_int_vector(sacc);
    free_int_vector(scount);
    free_double_vector(heat);
    free_double_vector(logLx);
    free_double_matrix(paramy,NC);
    free_double_tensor(history,NC,NH);
    free_double_tensor(Fisher,NC,NP);
    free_double_matrix(ejump,NC);
    free_double_tensor(evec,NC,NP);

    free(AAmp);
    free(EAmp);
    free(APhase);
    free(EPhase);
    free(FF);
    
}

void update(struct Data *dat, struct Het *het, int typ, int k, int ll, double *logLx, double **paramx, double **paramy, double **sx, double **sy, double *min, double *max, int *who, double *heat, double ***history, int NH, double **ejump, double ***evec, int **cv, int **av, gsl_rng *r)
{
    int q, i, j;
    int fstat;
    double qx, qy, a, b, c, aa, bb, tx, nhy;
    double alpha, beta, DL, pDx, pDy, H;
    double logLy, eta, leta, pMcMx, pMcMy, lglx;
    double w, v;
    double lm1, lm2, chi1, chi2;
    int flag, id1, id2;
    int Ndraw, Nden;
    long itheta, iphi, iDL, iiota, ipsi, iphic, itc, label;
    double dtheta, dphi, dDL, diota, dpsi, dphic, dtc;
    double cth, phi;
    double size, x, y, z;
    double Lf;
    double *jump, *zv;
    double *u;
    
    // set to 1 to use Fstat likelihood
    fstat = 0;
    
    u = double_vector(5);
    
    jump = double_vector(NI);
    zv = double_vector(NI);
    
    q = who[k];
    
    qx = qy = 0.0;    // log proposal densities
    
    if(typ == 0) // fisher jump
    {
        
        // pick an eigendirection to jump in
        beta = gsl_rng_uniform(r);
        i = (int)(beta*(NP));
        // draw the jump size
        beta = sqrt(heat[k])*ejump[q][i]*gsl_ran_gaussian(r,1.0);
        for(j = 0; j < NP; j++) paramy[q][j] = paramx[q][j]+beta*evec[q][i][j];
        
        tx = -1.0;
    }
    else if(typ == 1)// differential evolution
    {
        
        // the history file is kept for each temperature
        
        de_jump(paramx[q], paramy[q], history[k], NH, NP, r);
        
        tx = Tmerger(paramx[q],paramx[q][5]);
        
    }
    else if(typ == 2)// big sky
    {
             beta = gsl_rng_uniform(r);
             if(beta > 0.7)
             {
                 cth = -1.0+2.0*gsl_rng_uniform(r);
                 phi = TPI*gsl_rng_uniform(r);
             }
             else if (beta > 0.3)
             {
                 cth = paramx[q][7] + gsl_ran_gaussian(r,0.1);
                 phi = paramx[q][8] + gsl_ran_gaussian(r,0.1);
             }
             else
             {
                 cth = paramx[q][7] + gsl_ran_gaussian(r,0.05);
                 phi = paramx[q][8] + gsl_ran_gaussian(r,0.05);
             }
             
             if(phi < 0.0) phi += TPI;
             if(phi > TPI) phi -= TPI;
        
                paramy[q][7] = cth;
                paramy[q][8] = phi;

        if(fabs(cth) < 1.0)
          {
        
        int *pmap;
        double **Fish;
        double *params;
        double **Svec;
        double *Sval;
              
        Svec = double_matrix(4,4);
        Sval = double_vector(4);
        
        pmap = (int*)malloc(sizeof(int)* (NP));
        Fish = double_matrix(4,4);
        params = (double*)malloc(sizeof(double)* (NP));
        
        pmap[0] = pmap[1] = pmap[2] = pmap[3] = -1;
        pmap[5] = pmap[7] = pmap[8] = -1;
        pmap[4] = 0;
        pmap[6] = 1;
        pmap[9] = 2;
        pmap[10] = 3;
        
        
        for(j = 0; j < NP; j++) params[j] = paramx[q][j];
        tx = -1.0;  // time offset will be at correct value, no need to maximize
             
        lglx = Fstat_het(dat, het, ll, params, sx[q], tx);
        
        // The signal is invariant under these shifts. Find out which one gives
        // the smallest parameter deltas
        j = 0;
        u[0] = fabs(params[4]-paramx[q][4]) + fabs(params[9]-paramx[q][9]);
        u[1] = fabs(params[4]+PI/2-paramx[q][4]) + fabs(params[9]+PI/2-paramx[q][9]);
        u[2] = fabs(params[4]+PI/2-paramx[q][4]) + fabs(params[9]-PI/2-paramx[q][9]);
        u[3] = fabs(params[4]-PI/2-paramx[q][4]) + fabs(params[9]+PI/2-paramx[q][9]);
        u[4] = fabs(params[4]-PI/2-paramx[q][4]) + fabs(params[9]-PI/2-paramx[q][9]);
        
        j = 0;
        y = u[0];
        for (i=1; i< 5; i++)
        {
           if(u[i] < y)
            {
                y = u[i];
                j = i;
            }
        }
        
        if(j == 1)
        {
            params[4] += PI/2;
            params[9] += PI/2;
        }
        
        if(j == 2)
        {
            params[4] += PI/2;
            params[9] -= PI/2;
        }
        
        if(j == 3)
        {
            params[4] -= PI/2;
            params[9] += PI/2;
        }
        
        if(j == 4)
        {
            params[4] -= PI/2;
            params[9] -= PI/2;
        }
        
        //FisherSub(dat, ll, pmap, params, Fish);
        FisherSubHet(dat, het, ll, pmap, params, Fish);
        
        qx = 0.0;
        for (i=0; i< NP; i++)
        {
            if(pmap[i] > -1)
            {
              for (j=0; j< NP; j++)
               {
                if(pmap[j] > -1) qx -= 0.5*Fish[pmap[i]][pmap[j]]*(paramx[q][i]-params[i])*(paramx[q][j]-params[j]);
               }
            }
        }
        
         //y = logLx[q]-lglx;
        //printf("%f %f\n", y, qx);
        
        x = det(Fish,4)/(heat[k]*heat[k]*heat[k]*heat[k]);
        
        qx /= heat[k];
        
       // printf("%f %f\n", qx, 0.5*log(x));
        
        qx += 0.5*log(x);
        
        
        for(j = 0; j < NP; j++) params[j] = paramx[q][j];
        params[7] = cth;
        params[8] = phi;
        
        tx = Tmerger(paramx[q],paramx[q][5]);
             
        Lf = Fstat_het(dat, het, ll, params, sx[q], tx);
        
        //FisherSub(dat, ll, pmap, params, Fish);
        FisherSubHet(dat, het, ll, pmap, params, Fish);
        
        FisherEvec(Fish, Sval, Svec, 4);
        
        for(j = 0; j < NP; j++) paramy[q][j] = params[j];
        
        // pick an eigendirection to jump in
        beta = gsl_rng_uniform(r);
        i = (int)(beta*4);
        // draw the jump size
        beta = sqrt(heat[k])*Sval[i]*gsl_ran_gaussian(r,1.0);
        for(j = 0; j < NP; j++)
        {
          if(pmap[j] > -1) paramy[q][j] = params[j]+beta*Svec[i][pmap[j]];
        }
        
        qy = 0.0;
        for (i=0; i< NP; i++)
        {
            if(pmap[i] > -1)
            {
                for (j=0; j< NP; j++)
                {
                    if(pmap[j] > -1) qy -= 0.5*Fish[pmap[i]][pmap[j]]*(paramy[q][i]-params[i])*(paramy[q][j]-params[j]);
                }
            }
        }
        
        //printf("%f %f\n", v-w, qy);
        
        x = det(Fish,4)/(heat[k]*heat[k]*heat[k]*heat[k]);
        
        qy /= heat[k];
        
        // printf("%f %f\n", qy, 0.5*log(x));
        
        qy += 0.5*log(x);
        
       // printf("%f %f\n", qx, qy);
        
       // printf("%f %f %f %f %f ", qx, qy, lglx, logLx[q], Lf);
        
        
        // Need to include a Jacobian that accounts for the deteministic time mapping.
        // The Tmerger function is used to re-set tc so that the merger occurs at the
        // same time in the detector. Depending on the sky location, the time range dt
        // mapped to by a unit dcostheta dphi around the reference location will be different.
        
        // Test were unclear on wether this kind of term was needed
        
        // tvol is the time volume surrounding theta, phi
        //qy += log(tvol(paramy[q]));
        //qx += log(tvol(paramx[q]));
        
        
        // To cover the full range, apply a pi shift to 2psi and 2phic
        beta = gsl_rng_uniform(r);
        if(beta > 0.5)
        {
            paramy[q][4] += PI/2.0;
            paramy[q][9] += PI/2.0;
        }
              
        
        free_double_matrix(Fish,4);
        free(pmap);
        free(params);
        free_double_matrix(Svec,4);
        free_double_vector(Sval);
              
          }
        
         
         // Note that testing this proposal using logL=const requires is not very informative since the
         // Ftstat likelihood mostly gets skipped since the intrinsic parameters are not close to the true.
        
    }
    else  // update PSD scaling
    {
        beta = gsl_rng_uniform(r);
        if(beta > 0.7)
        {
            x = 0.1;
        }
        else if (beta > 0.3)
        {
            x = 0.01;
        }
        else
        {
            x = 0.001;
        }
        
        for(i = 0; i < dat->Nch; i++) sy[q][i] = sx[q][i] + gsl_ran_gaussian(r,x);
                                                  
    }
    
    
    // [0] ln Mass1  [1] ln Mass2  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln distance
    // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination

    if(ll == 0 || ll == 1)
    {
    if(paramy[q][1] > paramy[q][0])  // catch if m2 > m1 and flip
    {
        lm1 = paramy[q][1];
        chi1 = paramy[q][3];
        lm2 = paramy[q][0];
        chi2 = paramy[q][2];
        paramy[q][0] = lm1;
        paramy[q][1] = lm2;
        paramy[q][2] = chi1;
        paramy[q][3] = chi2;
    }
    }
    
    
    
    cv[typ][k]++;
    
    // [0] ln Mass1  [1] ln Mass2  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln distance
    // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
    

    // re-map angular parameters to their proper range
    if(paramy[q][4] > PI)   paramy[q][4] -= PI;
    if(paramy[q][4] < 0.0)  paramy[q][4] += PI;
    if(paramy[q][8] > TPI)  paramy[q][8] -= TPI;
    if(paramy[q][8] < 0.0)  paramy[q][8] += TPI;
    if(paramy[q][9] > PI)   paramy[q][9] -= PI;
    if(paramy[q][9] < 0.0)  paramy[q][9] += PI;

    
    // check proposed values are in prior range
    flag = 0;
    for(i = 0; i < NP; i++)
    {
        if(paramy[q][i] > max[i] || paramy[q][i] < min[i]) flag = 1;
    }
    
    if(nflag == 1)
    {
        for(i = 0; i < dat->Nch; i++) if(sy[q][i] < 0.25) flag = 1;
        for(i = 0; i < dat->Nch; i++) if(sy[q][i] > 4.0) flag = 1;
    }
     
    
    if(ll == 2 && flag == 0)
    {
     // eta cannot exceed 0.25
     leta = (5.0/3.0)*(paramy[q][0]-paramy[q][1]);
     if(leta > log(0.25)) flag = 1;
    
      // Jacobian that makes the prior flat in m1, m2.
      if(flag == 0)
      {
        eta = exp(leta);
        pMcMy = 2.0*paramy[q][1]+leta-0.5*log(1.0-4.0*eta);
        
        leta = (5.0/3.0)*(paramx[q][0]-paramx[q][1]);
        eta = exp(leta);
        pMcMx = 2.0*paramx[q][1]+leta-0.5*log(1.0-4.0*eta);
      }
        
    }
    
    if(flag == 0)
    {
        // Jacobian that makes the prior flat in m1, m2.
        // Jacobian is m1*m2, but we combine the probablities as logs
        if(ll == 0)
        {
        pMcMy = log(paramy[q][0]*paramy[q][1]);
        pMcMx = log(paramx[q][0]*paramx[q][1]);
        }

        if(ll == 1)
        {
        pMcMy = paramy[q][0]+paramy[q][1];
        pMcMx = paramx[q][0]+paramx[q][1];
        }
        
        logLy = 0.0;
        nhy = 0.0;
        
        if(lhold == 0)
        {
            //logLy = Likelihood(dat, ll, paramy[q]);
            logLy = log_likelihood_het(dat, het, ll, paramy[q], sy[q]);
        }
        
    
    // variable in MCMC is x=logD, so p(x) = dD/dx p(D) = D p(D) = D^3
    
    //pDx = 3.0*paramx[q][6];   // uniform in volume prior
    //pDy = 3.0*paramy[q][6];   // uniform in volume prior
    
        if(ll == 0)
        {
        pDx = log(paramx[q][6]);   // uniform in distance prior
        pDy = log(paramy[q][6]);   // uniform in distance prior
        }
        else
        {
        pDx = paramx[q][6];   // uniform in distance prior
        pDy = paramy[q][6];   // uniform in distance prior
        }
    
    
    H = (logLy-logLx[q])/heat[k] + pMcMy + pDy - qy - pDx - pMcMx + qx;
        
   // if(typ == 2) printf("%e %e\n", qx, qy);
    
    
    alpha = log(gsl_rng_uniform(r));

    
    if(H > alpha)
    {
        // copy over new state if accepted
        logLx[q] = logLy;
        for(i = 0; i < NP; i++) paramx[q][i] = paramy[q][i];
        if(nflag == 1) for(i = 0; i < dat->Nch; i++) sx[q][i] = sy[q][i];
        av[typ][k]++;
    }
        
    }
    
    free_double_vector(u);
    free_double_vector(jump);
    free_double_vector(zv);
    
}




