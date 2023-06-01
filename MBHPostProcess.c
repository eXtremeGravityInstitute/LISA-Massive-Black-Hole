/*******************************************************************************************
 
 Copyright (c) 2021 Tyson Littenberg
 
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

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include <gsl/gsl_rng.h>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#include <time.h>
#include "mbh.h"

#define FMIN 1e-5 //!<minimum frequency of reconstructions (Hz)
#define FMAX 0.02 //!<maximum frequency of reconstructions (Hz)
#define NTDI 2    //!<number of TDI channels 2={A,E}
#define NH 1000   //!<number of reconstructions used for statistics
#define NC 18     //!<number of columns in chain.dat

#define MAXSTRINGSIZE 512 //!<maximum number of characters in strings

int main(int argc,char *argv[])
{
    if(argc!=4)
    {
        printf("Usage: mbh_post /path/to/chain/file T /path/to/output/file\n");
        return 0;
    }
    /*
     parse command line
     argv[1] = /path/to/chain/file
     argv[2] = Tobs (s)
     argv[3] = /path/to/output/file
    */
    char filename[MAXSTRINGSIZE];
    sprintf(filename,argv[1]);
    FILE *chainfile = fopen(filename,"r");
    double Tobs = (double)atof(argv[2]);
    int NMAX = (int)(Tobs*(FMAX-FMIN));
    
    /*
     parse chain files
     */
    
    //first count samples in the chain file
    char* line;
    char lineBuffer[MAXSTRINGSIZE];
    int MBHSTEPS = 0;
    while((line = fgets(lineBuffer, MAXSTRINGSIZE, chainfile)) != NULL) MBHSTEPS++;
    MBHSTEPS--;
    rewind(chainfile);
    
    //check that there are enough samples
    //Note: Half of samples discarded as burn in.  Need at least NH left over
    if(MBHSTEPS < 2*NH)
    {
        fprintf(stderr,"ERROR: Not enough chain samples for post processing\n");
        fprintf(stderr,"  You need at least %i samples in the chain\n", NH*2);
        fprintf(stderr,"  Get more samples or change definition of NH\n");
        return 1;
    }

    //store parameter samples
    double **params = double_matrix(MBHSTEPS,NParams);
    char *column=NULL;
    double contents[NC];
    for(int i=0; i<MBHSTEPS; i++)
    {
        line=fgets(lineBuffer,MAXSTRINGSIZE,chainfile);
        
        column = strtok(line," ");
        for(int n=0; n<NC; n++)
        {
            sscanf(column, "%lg", &contents[n]);
            column=strtok(NULL," ");
        }
        
        
        for(int n=0; n<NParams; n++) params[i][n] = contents[n+2];
        
        //reparameterize for waveform generator
        map_params(2, params[i]);
    }
    fclose(chainfile);
    
    /*
     generate waveforms
     */
    
    //allocate memory for waveform generator
    struct MBH_Data *dat  = malloc(sizeof(struct MBH_Data));
    dat->Tobs = Tobs;
    dat->sqrtTobs = sqrt(dat->Tobs);
    dat->dt = 1./(2.*FMAX);
    dat->Nch = NTDI;
    dat->N = (int)(dat->Tobs/dat->dt);
    dat->data = double_matrix(dat->Nch,dat->N);
    dat->Tstart = 0.0;
    dat->Tend = dat->Tstart + dat->Tobs;

    //allocate frequency, phase, and amplitude arrays
    double *f = malloc(NMAX*sizeof(double));
    double *A_amp = malloc(NMAX*sizeof(double));
    double *E_amp = malloc(NMAX*sizeof(double));
    double *A_phi = malloc(NMAX*sizeof(double));
    double *E_phi = malloc(NMAX*sizeof(double));
    
    for(int i=0; i<NMAX; i++) f[i] = FMIN + (double)(i)/Tobs;

    //allocate memory for holding reconstructions
    double ***hrec = double_tensor(NMAX, NTDI, NH);
    for(int i=0; i<NMAX; i++)
        for(int j=0; j<NTDI; j++)
            for(int k=0; k<NH; k++)
                hrec[i][j][k] = 0.0;
    
    //downsample chain and fill hrec
    int j=0;
    fprintf(stdout,"Computing waveform reconstructions:\n");
    for(int i=MBHSTEPS/2; i<MBHSTEPS; i+=MBHSTEPS/2/NH)
    {
        fprintf(stdout,"   waveform %.4i\r",j);
        fflush(stdout);
        fullphaseamp(dat, 2, NMAX, params[i], f, A_amp, E_amp, A_phi, E_phi);
        
        for(int n=0; n<NMAX; n++)
        {
            
            double A_re = A_amp[n]*cos(A_phi[n]);
            double A_im = A_amp[n]*sin(A_phi[n]);
            double E_re = E_amp[n]*cos(E_phi[n]);
            double E_im = E_amp[n]*sin(E_phi[n]);
            
            hrec[n][0][j] = A_re*A_re+A_im*A_im;
            hrec[n][1][j] = E_re*E_re+E_im*E_im;
        }
        j++;
        if(j>=NH) break;
    }
    

    
    /*
     do statistics and write files
     */
    sprintf(filename,"%s",argv[3]);
    FILE *recfile=fopen(filename,"w");

    double A_med,A_lo_50,A_hi_50,A_lo_90,A_hi_90;
    double E_med,E_lo_50,E_hi_50,E_lo_90,E_hi_90;

    for(int n=0; n<NMAX; n++)
        for(int m=0; m<NTDI; m++)
            gsl_sort(hrec[n][m],1,NH);

    for(int i=0; i<NMAX; i++)
    {
        double f = (double)(i)/dat->Tobs;

        A_med   = gsl_stats_median_from_sorted_data   (hrec[i][0], 1, 100);
        A_lo_50 = gsl_stats_quantile_from_sorted_data (hrec[i][0], 1, 100, 0.25);
        A_hi_50 = gsl_stats_quantile_from_sorted_data (hrec[i][0], 1, 100, 0.75);
        A_lo_90 = gsl_stats_quantile_from_sorted_data (hrec[i][0], 1, 100, 0.05);
        A_hi_90 = gsl_stats_quantile_from_sorted_data (hrec[i][0], 1, 100, 0.95);

        E_med   = gsl_stats_median_from_sorted_data   (hrec[i][1], 1, 100);
        E_lo_50 = gsl_stats_quantile_from_sorted_data (hrec[i][1], 1, 100, 0.25);
        E_hi_50 = gsl_stats_quantile_from_sorted_data (hrec[i][1], 1, 100, 0.75);
        E_lo_90 = gsl_stats_quantile_from_sorted_data (hrec[i][1], 1, 100, 0.05);
        E_hi_90 = gsl_stats_quantile_from_sorted_data (hrec[i][1], 1, 100, 0.95);

        
        fprintf(recfile,"%.12g ",f);
        fprintf(recfile,"%.12lg ",A_med);
        fprintf(recfile,"%.12lg ",A_lo_50);
        fprintf(recfile,"%.12lg ",A_hi_50);
        fprintf(recfile,"%.12lg ",A_lo_90);
        fprintf(recfile,"%.12lg ",A_hi_90);
        fprintf(recfile,"%.12lg ",E_med);
        fprintf(recfile,"%.12lg ",E_lo_50);
        fprintf(recfile,"%.12lg ",E_hi_50);
        fprintf(recfile,"%.12lg ",E_lo_90);
        fprintf(recfile,"%.12lg ",E_hi_90);
        fprintf(recfile,"\n");

    }
    fclose(recfile);

    return 0;
}
