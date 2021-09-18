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
#include <omp.h>

#include "mbh.h"

#ifndef _OPENMP
#define omp ignore
#endif

// OSX
// clang -Xpreprocessor -fopenmp -lomp -w -o  PTMCMC PTMCMC.c Utils.c Response.c IMRPhenomD_internals.c IMRPhenomD.c -lgsl -lgslcblas  -lm

// Linux
// gcc -std=gnu99 -fopenmp -w -o PTMCMC PTMCMC.c IMRPhenomD_internals.c IMRPhenomD.c Utils.c Response.c -lgsl -lgslcblas  -lm

/* total number of chains */
#define NC 24

void MCMC(struct MBH_Data *dat, struct Het *het, int ll, int *who, double **params);

int main(int argc,char **argv)
{
    
    double f;
    char filename[50];
    double *params, *premove;
    double *AS, *ES;
    double x;
    long i, j, k, id, NS;
    
    struct Het *het = malloc(sizeof(struct Het));
    struct MBH_Data *dat = malloc(sizeof(struct MBH_Data));
    
    int seg, rep;
    
    const gsl_rng_type * P;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    
    P = gsl_rng_default;
    r = gsl_rng_alloc (P);
    
    FILE* in;
    FILE* out;
    
    params = (double*)malloc(sizeof(double)* (NParams));
    premove = (double*)malloc(sizeof(double)* (NParams));
    
    
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
            for(i=0; i< NParams; i++) fscanf(in,"%lf", &x);
            NS++;
        }
        
        rewind(in);
        
        // extract the source to be analyzed
        for(k=0; k < NS; k++)
        {
            if(k == rep)
            {
                fscanf(in,"%lf%lf", &x, &x);
                for(i=0; i< NParams; i++) fscanf(in,"%lf", &params[i]);
                //printf("%e\n", params[5]);
            }
            else // wind the file
            {
                fscanf(in,"%lf%lf", &x, &x);
                for(i=0; i< NParams; i++) fscanf(in,"%lf", &x);
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
                sprintf(filename, "specfit_%ld_%d.dat", id, seg);
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
                sprintf(filename, "specav_%ld.dat", id);
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
            for(i=0; i< NParams; i++) fscanf(in,"%lf", &premove[i]);
            
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
        for(i=0; i< NParams; i++) fscanf(in,"%lf", &params[i]);
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
    paramx = double_matrix(NC,NParams);
    
    // here we copy the best fit source into all the chains
    // in the global fit the previous state for each chain
    // will instead be passed in
    for (i=0; i< NC; i++)
    {
        for (j=0; j< NParams; j++) paramx[i][j] = params[j];
    }
    
    // set up the whos-who reference. In the global fit this will be passed
    // over from the previous state of the BH model
    int *who;
    who = int_vector(NC);
    for (i=0; i< NC; i++) who[i] = i;
    
    // perform the MCMC
    MCMC(dat, het, 2, who, paramx);
    
    
    free(params);
    
    return 1;
    
}

void MCMC(struct MBH_Data *dat, struct Het *het, int ll, int *who, double **paramx)
{
    double *logLx, SNR, x;
    double m1, m2;
    int i, j, k, q, mc, N;
    int NH=1000;
    int MP;
    int M = 2000000;
    double *heat;
    double **paramy;
    double ***Fisher;
    double ***history;
    double *min, *max;
    double **ejump, ***evec;
    double alpha,beta;
    double **sx, **sy;
    int *m;
    int *scount, *sacc, hold;
    int mcount, macc, typ;
    int **av, **cv;
    double Mc, Mtot, eta, dm;
    double a, b, c;
    
    double logLmax;
    double *pmax;
    
    FILE *chain;
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
    
    pmax = (double*)malloc(sizeof(double)* (NParams));
    
    sacc = int_vector(NC);
    scount = int_vector(NC);
    heat = double_vector(NC);
    logLx = double_vector(NC);
    paramy = double_matrix(NC,NParams);
    history = double_tensor(NC,NH,NParams);
    Fisher = double_tensor(NC,NParams,NParams);
    ejump = double_matrix(NC,NParams);
    evec = double_tensor(NC,NParams,NParams);
    
    // scale the PSDs
    sx = double_matrix(NC,dat->Nch);
    sy = double_matrix(NC,dat->Nch);
    
    // prior boundaries
    max = (double*)malloc(sizeof(double)* (NParams));
    min = (double*)malloc(sizeof(double)* (NParams));
    set_mbh_priors(dat, ll, min, max);
    
    // use one of the cold chains to produce the reference waveform
    het_space(dat, het, ll, paramx[who[0]], min, max);
    heterodyne(dat, het, ll, paramx[who[0]]);
    
    //initialize noise model
    for (i=0; i< NC; i++)
    {
        for (j=0; j< dat->Nch; j++) sx[i][j] = 1.0;
    }
    
    //initialize parameter chain buffer
    for (i=0; i< NC; i++)
    {
        for (k=0; k< NH; k++)
        {
            for (j=0; j< NParams; j++) history[i][k][j] = paramx[i][j];
        }
    }
    
    //initialize likelihood for each chain
    for (i=0; i< NC; i++)
    {
        logLx[i] = log_likelihood_het(dat, het, ll, paramx[i], sx[i]);
    }
    
    //store max likelihood & parameters
    logLmax = logLx[who[0]];
    for (j=0; j< NParams; j++) pmax[j] = paramx[who[0]][j];
    
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
        FisherEvec(Fisher[i], ejump[i], evec[i], NParams);
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
    
    //##############################################
    //open MP modifications
    gsl_rng **rvec;
    omp_set_num_threads(NC);
    rvec = (gsl_rng **)malloc(sizeof(gsl_rng *) * (NC+1));
    for(i = 0 ; i<= NC; i++){
        rvec[i] = gsl_rng_alloc(P);
        gsl_rng_set(rvec[i] , i);
    }
    //##############################################
    
    
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
                FisherEvec(Fisher[k], ejump[k], evec[k], NParams);
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
                for(i = 0; i < NParams; i++) paramy[j][i] = paramx[j][i];
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
                    for (j=0; j< NParams; j++) pmax[j] = paramx[q][j];
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
                    for(j=0; j<NParams; j++) history[k][i][j] = paramx[q][j];
                    m[k]++;
                }
            }
            
        }
        
        if(mc > 1000 && mc%20 == 0)
        {
             
            print_mbh_chain_file(dat, het, who, paramx, logLx, sx, ll, mc, chain);
            
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
            
            //get_component_masses(paramx[q], ll, &m1, &m2);
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
    free_double_tensor(Fisher,NC,NParams);
    free_double_matrix(ejump,NC);
    free_double_tensor(evec,NC,NParams);
    
    //###############################################
    //MT modification
    for(i =0 ;i<= NC; i++){
        gsl_rng_free(rvec[i]);
    }
    free(rvec);
    //###############################################
    
    
}
