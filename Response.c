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

void het_space(struct Data *dat, struct Het *het, int ll, double *params, double *min, double *max)
{
    int i, id, j, k, ii, jj, kk, MN, MM, jmin, flag, djmax, N, M, J, NR;
    int jmid, boost;
    int pflag;
    double **Ar, ***Ap, **Pr, ***Pp;
    double ***Cs, ***Sn, ***DH;
    double *dfa;
    double *fgrid;
    double fstop;
    char filename[1024];
    double k0, k1, k2;
    double norm0, norm1;
    double x, y, z;
    double tol, tl, tola, tla;
    double cmid, smid, dmid;
    FILE *out;
    double *ratio;
    int *stst;
    double **fisher, **cov;
    double **evec, *eval;
    double **dparams, *px;
    double alpha0, z0, zs, alphanew;
    double leta, eta;
    double alphax, dz, LDLmin, Lmax;
    double kxm, dfmax;
    
    N = dat->N;
    
    pflag = 1;  // print diagnostics if pflag = 1
    
    int Nch = dat->Nch;
    
    het->Nch = Nch;
    
    int NF;
    int NFmax = 100000;
    
    double *AF, *PF, *FF, *TF;
       
    FF = (double*)malloc(sizeof(double)* (NFmax));
       
    SetUp(dat, ll, params, NFmax, &NF, FF);
    
    x = SNRstart(dat, ll, params);
    
    if(x > FF[0]) FF[0] = x;
    
   // printf("%f %f\n", FF[0], FF[NF-1]);
    
    MN = (int)(ceil(FF[0]*dat->Tobs));
    MM = (int)(floor(FF[NF-1]*dat->Tobs));
    

    J = 2;
    tol = 0.00001;

    dparams = double_matrix(NP,NP);
    fisher = double_matrix(NP,NP);
    evec = double_matrix(NP,NP);
    eval = double_vector(NP);
    px = double_vector(NP);
    
    dfa = double_vector(N/2);
    
    FisherFast(dat, 2, params, fisher);
    FisherEvec(fisher, eval, evec, NP);
    efix(dat, het, 0, 2, params, min, max, eval, evec, 50.0);
    

    for (i = 0; i < NP; ++i)
    {
        for (j = 0; j < NP; ++j) px[j] = params[j] + eval[i]*evec[i][j];
        
         if(ll == 2)
          {
           leta = (5.0/3.0)*(px[0]-px[1]);
           eta = exp(leta);
           if(eta > 0.25)
           {
           for (j = 0; j < NP; ++j) px[j] = params[j] - alpha0*eval[i]*evec[i][j];
           }
          }
          
          // re-map angular parameters to their proper range
          
          x = px[4]/PI;
          px[4] = (x-floor(x))*PI;
          if(x < 0.0) px[4] += 1.0;
          x = px[8]/TPI;
          px[8] = (x-floor(x))*TPI;
          if(x < 0.0) px[8] += 1.0;
          x = px[9]/PI;
          px[9] = (x-floor(x))*PI;
          if(x < 0.0) px[9] += 1.0;
          
          while(px[4] > PI)  px[4] -= PI;
          while(px[4] < 0.0)  px[4] += PI;
          while(px[8] > TPI)  px[8] -= TPI;
          while(px[8] < 0.0)  px[8] += TPI;
          while(px[9] > PI)   px[9] -= PI;
          while(px[9] < 0.0)  px[9] += PI;
        
          for (j = 0; j < NP; ++j) if(px[j] > max[j]) px[j] = max[j];
          for (j = 0; j < NP; ++j) if(px[j] < min[j]) px[j] = min[j];
            
          if(ll == 2)
            {
            leta = (5.0/3.0)*(px[0]-px[1]);
            eta = exp(leta);
            if(eta > 0.25) px[0] = px[1] + 3.0/5.0*log(0.2499);
            if(eta < etamin) px[0] = px[1] + 3.0/5.0*log(etamin);
            }

        
        for (j = 0; j < NP; ++j) dparams[i][j] = px[j];
                   
        }
    
    
   
    
    fgrid = double_vector(2*NF);
    
    
    Ar = double_matrix(Nch,NF);
    Ap = double_tensor(NP,Nch,NF);
    Pr = double_matrix(Nch,NF);
    Pp = double_tensor(NP,Nch,NF);
    Cs = double_tensor(NP,Nch,NF);
    Sn = double_tensor(NP,Nch,NF);
    DH = double_tensor(NP,Nch,NF);
    
    TF = (double*)malloc(sizeof(double)* (NF));
    PF = (double*)malloc(sizeof(double)* (NF));
    AF = (double*)malloc(sizeof(double)* (NF));
    
    double phic, tc;
    
    // reference amplitude and phase
    fullphaseamp(dat, ll, NF, params, FF, Ar[0], Ar[1], Pr[0], Pr[1]);
    
    // perturbed amplitude and phase
    for(i = 0; i < NP; i++) fullphaseamp(dat, ll, NF, dparams[i], FF, Ap[i][0], Ap[i][1], Pp[i][0], Pp[i][1]);
    
    
    // phase difference vs frequency and ampltide ratio vs frequency
    for(i = 0; i < NP; i++)
    {
        for(id = 0; id < Nch; id++)
        {
         for(j = 0; j < NF; j++)
          {
              Cs[i][id][j] = 0.0;
              Sn[i][id][j] = 0.0;
              DH[i][id][j] = 0.0;
              
              if(Ar[id][j] > 0.0)
              {
              Cs[i][id][j] = 1.0-Ap[i][id][j]/Ar[id][j]*cos(Pr[id][j]-Pp[i][id][j]);
              Sn[i][id][j] = Ap[i][id][j]/Ar[id][j]*sin(Pr[id][j]-Pp[i][id][j]);
              DH[i][id][j] = Cs[i][id][j]*Cs[i][id][j]+Sn[i][id][j]*Sn[i][id][j];
              }
              
          }
        }
    }
    
    if(pflag == 1)
    {
    out = fopen("fspace.dat","w");
        
      for(j = 0; j < NF; j++)
        {
        fprintf(out, "%e ", FF[j]);
        for(i = 0; i < NP; i++) fprintf(out, "%e %e %e ", Cs[i][0][j], Sn[i][0][j], DH[i][0][j]);
        fprintf(out, "\n");
        }
        
        fclose(out);
        
    }
    
   
    
    dfmax = 1.0e-3;
 
    
    M = 0;
    fgrid[M] = FF[0];

        jmin = 0;
    
        for(j = 1; j < NF; j++)
        {
            
            flag = 0;
            
           // tl = tol*ratio[j];  // scale the tolerance by SNR contribution
           // if(tl < tol) tl = tol;
            
            tl = tol;
            
          // printf("%f %e\n", (double)(j)/Tobs, tl);
            
            // linear interpolation at the mid-point
            jmid = (j+jmin)/2;
            
           // printf("%d %d\n", j, jmid);
           
            for(id = 0; id < Nch; id++)
              {
                  for(i = 0; i < NP; i++)
                  {
                      cmid = (Cs[i][id][jmin]*(FF[j]-FF[jmid])+ Cs[i][id][j]*(FF[jmid]-FF[jmin]))/(FF[j]-FF[jmin]);
                      smid = (Sn[i][id][jmin]*(FF[j]-FF[jmid])+ Sn[i][id][j]*(FF[jmid]-FF[jmin]))/(FF[j]-FF[jmin]);
                      dmid = (DH[i][id][jmin]*(FF[j]-FF[jmid])+ DH[i][id][j]*(FF[jmid]-FF[jmin]))/(FF[j]-FF[jmin]);
                      
                      //printf("%d %d %d %e %e\n", j, id, i, Cs[i][id][jmid], cmid);
                      
                      if(fabs(Cs[i][id][jmid]-cmid) > tl || fabs(Sn[i][id][jmid]-smid) > tl || fabs(DH[i][id][jmid]-dmid) > tl)
                      {
                          flag = 1;
                          //printf("%d %d %d %f\n", j, id, i, tl);
                      }
                  }
              }
                
            
            // max frequency spacing
            if(FF[j]-FF[jmin] > dfmax) flag = 1;
     
            
            if(flag == 1) // re-start
            {
                jmin = j;
                M++;
                fgrid[M] = FF[j];
            }
            
        }
        
        
    
    if(pflag == 1)
       {
       out = fopen("df.dat","w");
       for (i = 0; i < M; ++i)
       {
           fprintf(out,"%e %e\n", fgrid[i], fgrid[i+1]-fgrid[i]);
       }
       fclose(out);
       }

    
    for (i = 0; i < N/2; ++i) dfa[i] = 1.0/dat->Tobs;
    
    for(ii = 0;  ii < M-1; ii++)
    {
        x = fgrid[ii+1]-fgrid[ii];
        
        kk = (int)(dat->Tobs*fgrid[ii]);
        jj = (int)(dat->Tobs*fgrid[ii+1]);
        
        for (i = kk; i < jj; ++i)
        {
           dfa[i] = x;
        }
        

     }

    
    for (i = jj; i < N/2; ++i) dfa[i] = dfmax;
    
    
    int *fgflat;
    
    fgflat = int_vector(2*M);
    
     j = MN;
     ii = 0;
     fgflat[ii] = j;
    
     do
     {
     jmin = (int)(rint(dat->Tobs*dfa[j]));
     //printf("%d %d %d %e\n", j, j+jmin*J, jmin, dfa[j]);
     do
     {
     flag = 0;
     k = j+jmin*J;
     for (i = j; i <= k; ++i)
     {
         if(i < N/2)
         {
          jj = (int)(rint(dat->Tobs*dfa[i]));
           if(jj < jmin)
           {
             jmin = jj;
             flag = 1;
           }
         }
         else
         {
             flag = 0;
         }
     }
     //printf("%d %d %d\n", j, j+jmin*J, jmin);
     }while(flag == 1);
    
     j = k;
      
     for (i = 0; i < J; ++i)
     {
         ii++;
         fgflat[ii] = fgflat[ii-1]+jmin;
        // printf("%d %d %d %d\n", ii, fgflat[ii], jmin, MM);
     }
         
         
     }while(j < MM);
    
    // fix the last section so it ends at MM
    jmin = (MM-fgflat[ii-J])/J;
    
    if(jmin >= 1)
    {

    ii -= J;
    
   // printf("\n%d %d %d\n", jmin, fgflat[ii-1], fgflat[ii-1]+J*jmin);
    
    for (i = 0; i < J; ++i)
    {
        ii++;
        fgflat[ii] = fgflat[ii-1]+jmin;
    }
        
    }
    
    
   // printf("%d %d %d\n", fgflat[ii], fgflat[ii-J], MM);
    
    if(pflag == 1)
    {
     out = fopen("df_flat.dat","w");
     for (i = 0; i < ii; ++i)
     {
         fprintf(out,"%e %e\n", (double)(fgflat[i])/dat->Tobs, (double)(fgflat[i+1]-fgflat[i])/dat->Tobs);
     }
     fclose(out);
    }
    
    
       M = ii+1;
       NR = ii/J;
       het->J = J;
       het->NR = NR;
       het->M = M;
       
       het->fgrid = int_vector(M);
       het->freq = double_vector(M);
       for (i = 0; i < M; ++i) het->fgrid[i] = fgflat[i];
       for (i = 0; i < M; ++i) het->freq[i] = (double)(fgflat[i])/dat->Tobs;
    
        het->MN = fgflat[0];
        het->MM = fgflat[M-1]+1;

    
    free(FF);
    free(TF);
    free(PF);
    free(AF);

    free_int_vector(fgflat);
    free_double_vector(dfa);
    //free_double_vector(ratio);
    //free_int_vector(stst);
    free_double_matrix(dparams,NP);
    free_double_matrix(fisher,NP);
    free_double_matrix(evec,NP);
    free_double_vector(eval);
    free_double_vector(px);
    free_double_vector(fgrid);
    free_double_matrix(Ar,Nch);
    free_double_tensor(Ap,NP,Nch);
    free_double_matrix(Pr,Nch);
    free_double_tensor(Pp,NP,Nch);
    free_double_tensor(Cs,NP,Nch);
    free_double_tensor(Sn,NP,Nch);
    free_double_tensor(DH,NP,Nch);
    
}

void fullphaseamp(struct Data *dat, int ll, int K, double *params, double *freq, double *Aamp, *Eamp, *Aphase, *Ephase)
{
    // full frequency template, amplitude and phase
    
    int i;
    double tc, phic, kxm;
    double *TF, *PF, *AF;
    FILE *out;
    
    TF = (double*)malloc(sizeof(double)* (K));
    PF = (double*)malloc(sizeof(double)* (K));
    AF = (double*)malloc(sizeof(double)* (K));
    
    // reference amplitude and phase
    Intrinsic(ll, params, dat->Tobs, K, freq, TF, PF, AF);
    
    // Have to add in the merger time and phase since they are not in the intrinsic or extrinsic subroutines
    phic = params[4];  // merger phase
    tc = params[5];    // merger time
    for(i=0; i< K; i++)
    {
        PF[i] -= TPI*freq[i]*(dat->Tobs-tc + dat->dt/2.0)-2.0*phic;
    }
    
    Extrinsic(params, dat->Tstart, dat->Tend, K, freq, TF, PF, AF, Aamp, Eamp, Aphase, Ephase, &kxm);
    
    free(TF);
    free(PF);
    free(AF);
    
}

void freehet(struct Het *het)
{
    free_double_tensor(het->IP, het->NR,1+het->J);
    free_double_tensor(het->SL,het->NR,het->Nch);
    free_double_matrix(het->amp,het->Nch);
    free_double_matrix(het->phase,het->Nch);
    free_double_matrix(het->rc,het->Nch);
    free_double_matrix(het->rs,het->Nch);
    free_double_matrix(het->dc,het->Nch);
    free_double_matrix(het->ds,het->Nch);
    free_double_matrix(het->aa,het->Nch);
    free_double_tensor(het->lc,het->NR,het->Nch);
    free_double_tensor(het->ls,het->NR,het->Nch);
    free_double_tensor(het->ldc,het->NR,het->Nch);
    free_double_tensor(het->lds,het->NR,het->Nch);
    free_double_vector(het->pref);
    free_int_vector(het->fgrid);
    free_double_vector(het->freq);
    free_double_vector(het->logLR);
    free_double_vector(het->DD);
}

void heterodyne(struct Data *dat, struct Het *het, int ll, double *params)
{
    int i, id, j, k, n;
    double f, fstop, x, y, z;
    double **hb, **rb, **ampb, **phaseb;
    double cp, sp;
    double Tobs;
    double **RC, **RS, **AA;
    double **DC, **DS;
    double logL, HH, HR;
    double kxm;
    double phic, tc, dt;
    int M, N, MM, MN, NR, ds, J, U;
    int Nch = het->Nch;
    FILE *out;
    
    M = het->M;
    MM = het->MM;
    MN = het->MN;
    J = het->J;
    NR = het->NR;
    
    N = dat->N;
    
    Tobs = dat->Tobs;
  
    het->IP = double_tensor(NR,J+1,J+1);
    het->SL = double_tensor(NR,Nch,J+1);
    het->amp = double_matrix(Nch,M);
    het->phase = double_matrix(Nch,M);
    het->rc = double_matrix(Nch,NR);
    het->rs = double_matrix(Nch,NR);
    het->dc = double_matrix(Nch,NR);
    het->ds = double_matrix(Nch,NR);
    het->aa = double_matrix(Nch,NR);
    het->lc = double_tensor(NR,Nch,J+1);
    het->ls = double_tensor(NR,Nch,J+1);
    het->ldc = double_tensor(NR,Nch,J+1);
    het->lds = double_tensor(NR,Nch,J+1);
    het->pref = double_vector(NP);  // keep a copy of the source parameters used for the heterodyne
    het->logLR = double_vector(Nch); // reference likelihood in each channel
    het->DD = double_vector(Nch); // data inner product (used when fitting noise level)
    
    for(j = 0; j < NP; j++) het->pref[j] = params[j];
    
    hb = double_matrix(Nch,N);
    rb = double_matrix(Nch,N);
    RS = double_matrix(Nch,N);
    RC = double_matrix(Nch,N);
    DS = double_matrix(Nch,N);
    DC = double_matrix(Nch,N);
    AA = double_matrix(Nch,N);
    ampb = double_matrix(Nch,MM);
    phaseb = double_matrix(Nch,MM);
    
    // full frequency template, amplitude and phase
    
    double *FF;
    
    FF = double_vector(MM);
    FF[0] = 1.0/dat->Tobs;
    for (i = 1; i < MM; ++i) FF[i] = (double)(i)/dat->Tobs;
    
    fullphaseamp(dat, ll, MM, params, FF, ampb[0], ampb[1], phaseb[0], phaseb[1]);
    
    /*
    out = fopen("wtf.dat","w");
    for (i = 1; i < MM; ++i) fprintf(out,"%e %e %e %e %e\n", FF[i], ampb[0][i], ampb[1][i], phaseb[0][i], phaseb[1][i]);
    fclose(out);
    */
    
        for(id = 0; id < Nch; id++)
           {
           for (i = 0; i < N; ++i)
           {
               hb[id][i] = 0.0;
           }
           }
    
      for(id = 0; id < Nch; id++)
        {
        for (i = 1; i < MM; ++i)
        {
            hb[id][i] = ampb[id][i]*cos(phaseb[id][i]);
            hb[id][N-i] = ampb[id][i]*sin(phaseb[id][i]);
        }
        }
    
    // store the downsampled reference amplitude and phase
    //fullphaseamp(dat, ll, M, params, het->freq, het->amp[0], het->amp[1], het->phase[0], het->phase[1]);
    
    
       for(id = 0; id < Nch; id++)
           {
           for (i = 0; i < M; ++i)
           {
               n = het->fgrid[i];
               //printf("%d %d %e %e %e\n", i, n, FF[n], het->freq[i], FF[n]-het->freq[i]);
               //printf("%d %e %e %e\n", i, het->amp[id][i], ampb[id][n], het->amp[id][i]-ampb[id][n]);
               het->amp[id][i] = ampb[id][n];
               het->phase[id][i] = phaseb[id][n];
           }
           }

    /*
    out = fopen("wtf2.dat","w");
    for (i = 0; i < M; ++i) fprintf(out,"%e %e %e\n", het->freq[i], het->amp[0][i], het->amp[1][i]);
    fclose(out);
    */
    
     // form up resdiuals
      for(id = 0; id < Nch; id++)
       {
        for(j = 0; j < N; j++) rb[id][j] = dat->data[id][j] - hb[id][j];
       }

    
    logL = 0.0;
    x = 0.0;
    for(id = 0; id < Nch; id++)
    {
        HH = fourier_nwip2(hb[id], hb[id], dat->SN[id], MN, MM, N);
        HR = fourier_nwip2(rb[id], hb[id], dat->SN[id], MN, MM, N);
        het->DD[id] = fourier_nwip2(dat->data[id], dat->data[id], dat->SN[id], MN, MM, N);
        het->logLR[id] = HR+0.5*HH;
        logL += het->logLR[id];
        x += HH;
        //printf("HetRef %d %f %f\n", id, HH, HR);
    }
    
    het->SNR = sqrt(x);
    
    printf("ref logL = %e\n", logL);
    
    for(id = 0; id < Nch; id++)
       {
         for(j = 0; j < N; j++)
         {
            RC[id][j] = 0.0;
            RS[id][j] = 0.0;
            DC[id][j] = 0.0;
            DS[id][j] = 0.0;
            AA[id][j] = 0.0;
         }
       }

    
    for(id = 0; id < Nch; id++)
    {
      for(j = 1; j < MM; j++)
      {
        cp = cos(phaseb[id][j]);
        sp = sin(phaseb[id][j]);
        x = ampb[id][j]/dat->SN[id][j];
        RC[id][j] = 2.0*(rb[id][j]*cp+rb[id][N-j]*sp)*x;
        RS[id][j] = 2.0*(-rb[id][j]*sp+rb[id][N-j]*cp)*x;
        AA[id][j] = ampb[id][j]*ampb[id][j]/dat->SN[id][j];
        DC[id][j] = 2.0*(dat->data[id][j]*cp+dat->data[id][N-j]*sp)*x;
        DS[id][j] = 2.0*(-dat->data[id][j]*sp+dat->data[id][N-j]*cp)*x;
      }
    }
    

    
    // store boundary values for correcting overcount
    for(id = 0; id < Nch; id++)
    {
        for (i = 0; i < NR; ++i)
          {
            k = het->fgrid[i*J];
            het->rc[id][i] = RC[id][k];
            het->rs[id][i] = RS[id][k];
            het->dc[id][i] = DC[id][k];
            het->ds[id][i] = DS[id][k];
            het->aa[id][i] = AA[id][k];
        }
    }
    

    
    double **P, **PP;
    
    PP = double_matrix(J+1,J+1); // downsampled Legendres
    
    for (i = 0; i < NR; ++i)
    {
      for(id = 0; id < Nch; id++)
        {
         for(j = 0; j <= J; j++)
         {
             het->SL[i][id][j] = 0.0;
             het->lc[i][id][j] = 0.0;
             het->ls[i][id][j] = 0.0;
             het->ldc[i][id][j] = 0.0;
             het->lds[i][id][j] = 0.0;
         }
        }
    }

    
    for (i = 0; i < NR; ++i)
    {
        //printf("%d %d %d\n", (i+1)*J, M, het->fgrid[(i+1)*J]);
        
        if((i+1)*J < M)
        {
            
         U = het->fgrid[(i+1)*J]-het->fgrid[i*J];
         
         P = double_matrix(J+1,U+1);
        
         legendre_maker(J, U, P);
        
        // downsampled Legendre
        for (j = 0; j <= J; ++j)
          {
            for (k = 0; k <= J; ++k)
              {
                  PP[j][k] = P[j][het->fgrid[i*J+k]-het->fgrid[i*J]];
              }
          }
        
        // store inverse of downsampled Legendre
        Inverse(PP, het->IP[i], J+1);
            
            /*
            printf("%d\n", i);
            for (j = 0; j <= J; ++j)
            {
              for (k = 0; k <= J; ++k)
                {
                    printf("%f ", het->IP[i][j][k]);
                }
                printf("\n");
            }
            printf("\n");
             */
        
        for(id = 0; id < Nch; id++)
        {
          for(n = 0; n <= U; n++)
            {
                k = het->fgrid[i*J]+n;
               for(j = 0; j <= J; j++)
                {
                    het->SL[i][id][j] += AA[id][k]*P[j][n];
                    het->lc[i][id][j] += RC[id][k]*P[j][n];
                    het->ls[i][id][j] += RS[id][k]*P[j][n];
                    het->ldc[i][id][j] += DC[id][k]*P[j][n];
                    het->lds[i][id][j] += DS[id][k]*P[j][n];
                }
            }
         }
        
        
        free_double_matrix(P,J+1);
            
        }
    }
    
    //printf("done with setup\n");
    
    free_double_matrix(ampb,Nch);
    free_double_matrix(phaseb,Nch);
    free_double_vector(FF);
    free_double_matrix(PP,J+1);
    free_double_matrix(RC,Nch);
    free_double_matrix(RS,Nch);
    free_double_matrix(DC,Nch);
    free_double_matrix(DS,Nch);
    free_double_matrix(AA,Nch);
    free_double_matrix(hb,Nch);
    free_double_matrix(rb,Nch);
    
        
}


double log_likelihood_het(struct Data *dat, struct Het *het, int ll, double *params, double *sx)
{
    double **amp, **phase;
    double **hc, **hs;
    double HH, HR, logL;
    double x, y;
    double L0,L1;
    double **cc, **ss;
    double *uu, *vv;
    int i, j, k, id, M, J, NR;
    int Nch = het->Nch;
    
    FILE *out;
    
    logL = 0.0;
    
    if(lhold == 0)
    {
    
    M = het->M;
    J = het->J;
    NR = het->NR;

    amp = double_matrix(Nch,M);
    phase = double_matrix(Nch,M);
    cc = double_matrix(Nch,M);
    ss = double_matrix(Nch,M);
    hc = double_matrix(Nch,M);
    hs = double_matrix(Nch,M);
    uu = double_vector(J+1);
    vv = double_vector(J+1);
  
        
    fullphaseamp(dat, ll, M, params, het->freq, amp[0], amp[1], phase[0], phase[1]);
    
    //out = fopen("check.dat","w");
    for(id = 0; id < Nch; id++)
    {
        for(j = 0; j < M; j++)
        {
            cc[id][j] = 0.0;
            hc[id][j] = 0.0;
            hs[id][j] = 0.0;
            if(het->amp[id][j] > 0.0)
            {
            cc[id][j] = 2.0*(het->amp[id][j]-amp[id][j]*cos(het->phase[id][j]-phase[id][j]));
            ss[id][j] = 2.0*amp[id][j]*sin(het->phase[id][j]-phase[id][j]);
            hc[id][j] = cc[id][j]/het->amp[id][j];
            hs[id][j] = ss[id][j]/het->amp[id][j];
            cc[id][j] = hc[id][j]*hc[id][j]+hs[id][j]*hs[id][j];
            }
            //if(id == 0) fprintf(out,"%d %.12e %e %f %f %f\n", j, het->freq[j], het->amp[id][j], hc[id][j], hs[id][j], cc[id][j]);
        }
    }
       // fclose(out);
        
        logL = 0.0;
        
        for(id = 0; id < Nch; id++)
        {
          HH = 0.0;
          HR = 0.0;
       
            for(i = 0; i < NR; i++)
            {
                
                for(j = 0; j <= J; j++)
                {
                  vv[j] = cc[id][i*J+j];
                  uu[j] = 0.0;
                }
                
                // get Legendre coefficicients for slow term
                for(j = 0; j <= J; j++)
                    {
                     for(k = 0; k <= J; k++)
                        {
                            uu[j] += het->IP[i][j][k]*vv[k];
                        }
                    }
                
                for(j = 0; j <= J; j++) HH += uu[j]*het->SL[i][id][j];
                if(i < NR-2) HH -= het->aa[id][i+1]*vv[J];  // correction for overcount
                
                
            if(LDC == 1)  // don't need these if noise free
              {
                    
                for(j = 0; j <= J; j++)
                {
                  vv[j] = hc[id][i*J+j];
                  uu[j] = 0.0;
                }
                
                for(j = 0; j <= J; j++)
                    {
                     for(k = 0; k <= J; k++)
                        {
                            uu[j] += het->IP[i][j][k]*vv[k];
                        }
                    }
                
                for(j = 0; j <= J; j++) HR += uu[j]*het->lc[i][id][j];
                
                if(i < NR-2) HR -= het->rc[id][i+1]*vv[J];  // correction for overcount
                
                
                    for(j = 0; j <= J; j++)
                           {
                             vv[j] = hs[id][i*J+j];
                             uu[j] = 0.0;
                           }
                           
                           for(j = 0; j <= J; j++)
                               {
                                for(k = 0; k <= J; k++)
                                   {
                                       uu[j] += het->IP[i][j][k]*vv[k];
                                   }
                               }
                           
                           for(j = 0; j <= J; j++) HR += uu[j]*het->ls[i][id][j];
                           
                           if(i < NR-2) HR -= het->rs[id][i+1]*vv[J];  // correction for overcount
                  
              }
                
               // printf("%d %e %e %e %e\n", i, HH, HR, x, y);
        }
    
         if(nflag == 1)
         {
             logL += (het->logLR[id]-(HR+0.5*HH)-0.5*het->DD[id])/sx[id];
         }
         else
         {
           logL += het->logLR[id]-(HR+0.5*HH);
         }
        
       // printf("Het %d %f %f %f %f\n", id, HH, HR, x, y);
    }
       
      if(nflag == 1)
      {
      for(id = 0; id < Nch; id++)
       {
           logL -= (double)(het->MM-het->MN)*log(sx[id]);
       }
      }
    
    free_double_vector(uu);
    free_double_vector(vv);
    free_double_matrix(cc,Nch);
    free_double_matrix(ss,Nch);
    free_double_matrix(hc,Nch);
    free_double_matrix(hs,Nch);
    free_double_matrix(amp,Nch);
    free_double_matrix(phase,Nch);
        
    }
    
    return(logL);
    
}

double Fstat_het(struct Data *dat, struct Het *het, int ll, double *params, double *sx, double tm)
{
    double ***amp, ***phase;
    double **Fparams;
    double HH, HR, logL;
    double cosi, scale, psi, phic, Ac, Ap, A;
    double x, y, u, v;
    double QQ;
    double L0,L1,tx;
    double ***cc, ***ss;
    double *uu, *vv;
    int i, j, k, id, M, J, NR;
    int ii, jj;
    int Nch = het->Nch;
    
    FILE *out;
    
    logL = 0.0;
    
    if(lhold == 0)
    {
    
    M = het->M;
    J = het->J;
    NR = het->NR;

    Fparams = double_matrix(4,NP);
    amp = double_tensor(4,Nch,M);
    phase = double_tensor(4,Nch,M);
        
    // time maximize
    if(tm > 0.0)
    {
     tx = Tmerger(params,params[5]);
     params[5] += (tm-tx);
    }
        
    for(i = 0; i < 4; i++)
      {
          for(j = 0; j < NP; j++) Fparams[i][j] = params[j];
          Fparams[i][10] = 0.0;
      }
        

     
    Fparams[0][4] = 0.0;
    Fparams[0][9] = 0.0;
    Fparams[1][4] = 0.0;
    Fparams[1][9] = PI/4.0;
    Fparams[2][4] = PI/4.0;
    Fparams[2][9] = 0.0;
    Fparams[3][4] = PI/4.0;
    Fparams[3][9] = PI/4.0;
        
    for(i = 0; i < 4; i++) fullphaseamp(dat, ll, M, Fparams[i], het->freq, amp[i][0], amp[i][1], phase[i][0], phase[i][1]);
        
    // Fstat convention
     for(i = 0; i < M; i++)
     {
         for(id = 0; id < Nch; id++)
              {
               amp[1][id][i] *= -1.0;
               amp[2][id][i] *= -1.0;
              }
     }
        
    cc = double_tensor(4,Nch,M);
    ss = double_tensor(4,Nch,M);
   
    uu = double_vector(J+1);
    vv = double_vector(J+1);
  
    
    for(i = 0; i < 4; i++)
     {
      for(id = 0; id < Nch; id++)
       {
        for(j = 0; j < M; j++)
        {
            cc[i][id][j] = 0.0;
            ss[i][id][j] = 0.0;
            if(het->amp[id][j] > 0.0)
            {
            cc[i][id][j] = 2.0*amp[i][id][j]*cos(het->phase[id][j]-phase[i][id][j])/het->amp[id][j];
            ss[i][id][j] = 2.0*amp[i][id][j]*sin(het->phase[id][j]-phase[i][id][j])/het->amp[id][j];
            }
        }
      }
     }
     
    double **MM, **MI;
    double *NV, *aV;
        
    NV = (double*)malloc(sizeof(double)* (4));
    aV = (double*)malloc(sizeof(double)* (4));
    MM = double_matrix(4,4);
    MI = double_matrix(4,4);
        
        
        
    for(ii = 0; ii < 4; ii++)
        {
            NV[ii] = 0.0;
        for(jj = 0; jj < 4; jj++)
         {
             MM[ii][jj]  = 0.0;
         }
        }
        
    for(id = 0; id < Nch; id++)
     {
      for(ii = 0; ii < 4; ii++)
        {
         QQ = 0.0;
          for(i = 0; i < NR; i++)
            {
                for(j = 0; j <= J; j++)
                {
                  vv[j] = cc[ii][id][i*J+j];
                  uu[j] = 0.0;
                }
                
                for(j = 0; j <= J; j++)
                    {
                     for(k = 0; k <= J; k++)
                        {
                            uu[j] += het->IP[i][j][k]*vv[k];
                        }
                    }
                
                for(j = 0; j <= J; j++) QQ += uu[j]*het->ldc[i][id][j];
                
                if(i < NR-2) QQ -= het->dc[id][i+1]*vv[J];  // correction for overcount
                
                
                    for(j = 0; j <= J; j++)
                           {
                             vv[j] = ss[ii][id][i*J+j];
                             uu[j] = 0.0;
                           }
                           
                           for(j = 0; j <= J; j++)
                               {
                                for(k = 0; k <= J; k++)
                                   {
                                       uu[j] += het->IP[i][j][k]*vv[k];
                                   }
                               }
                           
                           for(j = 0; j <= J; j++) QQ += uu[j]*het->lds[i][id][j];
                           
                           if(i < NR-2) QQ -= het->ds[id][i+1]*vv[J];  // correction for overcount
                
           } // end i
            if(nflag == 1)
            {
                NV[ii] += QQ/sx[id];
            }
            else
            {
                NV[ii] += QQ;
            }
         }  // end ii
      }   // end id
        
      for(id = 0; id < Nch; id++)
          {
          for(ii = 0; ii < 4; ii++)
             {
             for(jj = ii; jj < 4; jj++)
                {
                   QQ = 0.0;
                   for(i = 0; i < NR; i++)
                     {
                       
                       for(j = 0; j <= J; j++)
                        {
                         vv[j] = cc[ii][id][i*J+j]*cc[jj][id][i*J+j]+ss[ii][id][i*J+j]*ss[jj][id][i*J+j];
                         uu[j] = 0.0;
                        }
                                   
                        for(j = 0; j <= J; j++)
                          {
                         for(k = 0; k <= J; k++)
                            {
                            uu[j] += het->IP[i][j][k]*vv[k];
                            }
                         }
                        for(j = 0; j <= J; j++) QQ += uu[j]*het->SL[i][id][j];;
                        if(i < NR-2) QQ -= het->aa[id][i+1]*vv[J];  // correction for overcount
                       }
                    
                       if(nflag == 1)
                        {
                        MM[ii][jj] += QQ/sx[id];
                        }
                       else
                       {
                        MM[ii][jj] += QQ;
                       }
                }
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
                aV[i] += MI[i][j]*NV[j];
                logL += 0.5*MI[i][j]*NV[i]*NV[j];
            }
        }
        
        if(nflag == 1)
        {
        for(id = 0; id < Nch; id++)
         {
             logL -= 0.5*het->DD[id]/sx[id];
             logL -= (double)(het->MM-het->MN)*log(sx[id]);
         }
        }
        
        //printf("%e %f\n", logL, sqrt(2.0*logL));
        
        //printf("%e %e %e %e\n", aV[0], aV[1], aV[2], aV[3]);
        
        x = (aV[0]+aV[3]);
        x = x*x;
        y = (aV[1]-aV[2]);
        y = y*y;
        u = (aV[0]-aV[3]);
        u = u*u;
        v = (aV[1]+aV[2]);
        v = v*v;
    
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
        
        scale = 4.0/A;
        
     // printf("%f %f %f %f %f %f %f %f %f\n", scale, x, y, phic, psi, cosi, params[4], params[9], params[10]);
        
        
        params[10] = cosi;
        params[6] += log(scale);
        params[9] = psi;
        params[4] = phic;
         
    
    free(NV);
    free(aV);
    free_double_matrix(MM,4);
    free_double_matrix(MI,4);
    free_double_vector(uu);
    free_double_vector(vv);
    free_double_tensor(cc,4,Nch);
    free_double_tensor(ss,4,Nch);
    free_double_tensor(amp,4,Nch);
    free_double_tensor(phase,4,Nch);
    free_double_matrix(Fparams,4);
         
    }
    
    return(logL);
    
}

// Uses the full data (slow)
void FstatFull(struct Data *dat, int ll, double *params, double *pnew)
{
    
    /*   Indicies   */
    int i,j, k, n, m, imax;
    
    /*   Miscellaneous  */
    double xm, fstep, power, om, mx, tx;
    
    double *pfilt;
    
    double **filtA, **filtE;
    
    double **MM, **MI;
    
    double *NV, *aV;
    
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
    
    
    clock_t start, end;
    double cpu_time_used;
    
    pfilt = (double*)malloc(sizeof(double)* (NP));
    
    filtA = double_matrix(4,dat->N);
    filtE = double_matrix(4,dat->N);
    
    NV = (double*)malloc(sizeof(double)* (4));
    aV = (double*)malloc(sizeof(double)* (4));
    MM = double_matrix(4,4);
    MI = double_matrix(4,4);
    
    for (i=0; i< NP; i++) pfilt[i] = params[i];
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
       // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
    
    pfilt[10] = 0.0;
    
    pfilt[4] = 0.0;
    pfilt[9] = 0.0;
    ResponseFreq(dat, ll, pfilt, filtA[0], filtE[0]);
    
    
    pfilt[4] = 0.0;
    pfilt[9] = PI/4.0;
    ResponseFreq(dat, ll, pfilt, filtA[1], filtE[1]);
    
    
    pfilt[4] = PI/4.0;
    pfilt[9] = 0.0;
    ResponseFreq(dat, ll, pfilt, filtA[2], filtE[2]);
    
    
    pfilt[4] = PI/4.0;
    pfilt[9] = PI/4.0;
    ResponseFreq(dat, ll, pfilt, filtA[3], filtE[3]);
    
    
    for (i=0; i< dat->N; i++)
    {
        filtA[1][i] *= -1.0;
        filtE[1][i] *= -1.0;
        filtA[2][i] *= -1.0;
        filtE[2][i] *= -1.0;
    }
    
   //printf("\n");
        for (i=0; i< 4; i++)
        {
            NV[i] = 2.0*(fourier_nwip(filtA[i], dat->data[0], dat->SN[0], dat->N) + fourier_nwip(filtE[i], dat->data[1], dat->SN[1], dat->N));
            for (j=i; j< 4; j++)
            {
                MM[i][j] =  4.0*(fourier_nwip(filtA[i], filtA[j], dat->SN[0], dat->N) + fourier_nwip(filtE[i], filtE[j], dat->SN[1], dat->N));
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
        
        Inverse(MM, MI, 4);
        
        logL = 0.0;
        
        for (i=0; i< 4; i++)
        {
            aV[i] = 0.0;
            for (j=0; j< 4; j++)
            {
                aV[i] += MI[i][j]*NV[j];
                logL += 0.5*MI[i][j]*NV[i]*NV[j];
            }
        }
        
        printf("Fstat %e %f\n", logL, sqrt(2.0*logL));
        
        //printf("%e %e %e %e\n", aV[0], aV[1], aV[2], aV[3]);
        
        x = (aV[0]+aV[3]);
        x = x*x;
        y = (aV[1]-aV[2]);
        y = y*y;
        u = (aV[0]-aV[3]);
        u = u*u;
        v = (aV[1]+aV[2]);
        v = v*v;
    
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
        
       printf("%f %f %f %f %f %f %f %f %f\n", scale, x, y, phic, psi, cosi, params[4], params[9], params[10]);
    
       // maxed on data without noise get
       // scale = 0.996465 phic = 2.837624 psi =  2.930106 cosi = 0.339382
        // maxed on data with noise get
        // scale = 0.992009 phic = 2.835180 psi = 2.927718 cosi = 0.337941
    
        // reference values
        // 1 2.848705 2.937147 0.339386
    
      for (i=0; i< NP; i++) pnew[i] = params[i];
    
        pnew[10] = cosi;
        pnew[6] += log(scale);
        pnew[9] = psi;
        pnew[4] = phic;

    /*   Deallocate Arrays   */
    
     free(NV);
     free(aV);
     free_double_matrix(MM,4);
     free_double_matrix(MI,4);
     free_double_matrix(filtA,4);
     free_double_matrix(filtE,4);
    
    
    
    return;
}





double SNRFast(struct Data *dat, int ll, double *params)
{
    double *AAmp, *EAmp, *APhase, *EPhase;
    double af, fr, df, DT, fac, deltaF, f, fmax, fmin, x, y, t, t10;
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double Amp, Phase, tf, fguess, kxm;
    double m1_SI, m2_SI, distance, tc, phic;
    int i, j, NF, NW, flag;
    double A, P;
    double SNR;
    double px, fnew;
    double *Aint, *Eint;
    double *SNRSQA, *SNRSQE;

    FILE *out;
    
    int NFmax = 100000;
 
    double *AF, *PF, *FF, *TF;
    
    //clock_t start, end;
    //double cpu_time_used;
    
    FF = (double*)malloc(sizeof(double)* (NFmax));
    
    SetUp(dat, ll, params, NFmax, &NF, FF);
                         
    AF = (double*)malloc(sizeof(double)* (NF));
    PF = (double*)malloc(sizeof(double)* (NF));
    TF = (double*)malloc(sizeof(double)* (NF));
    
    Intrinsic(ll, params, dat->Tobs, NF, FF, TF, PF, AF);
    
    AAmp = (double*)malloc(sizeof(double)* (NF));
    EAmp = (double*)malloc(sizeof(double)* (NF));
    APhase = (double*)malloc(sizeof(double)* (NF));
    EPhase = (double*)malloc(sizeof(double)* (NF));
                         
    Extrinsic(params, dat->Tstart, dat->Tend, NF, FF, TF, PF, AF, AAmp, EAmp, APhase, EPhase, &kxm);
    
    Aint = (double*)malloc(sizeof(double)* (NF));
    Eint = (double*)malloc(sizeof(double)* (NF));
    
    out = fopen("intg.dat","w");
    for (i=0; i< NF; i++)
    {
        j = (int)(FF[i]*dat->Tobs);
        //printf("%d %d\n", i, j);
        Aint[i] = 4.0*AAmp[i]*AAmp[i]/dat->SM[0][j];
        Eint[i] = 4.0*EAmp[i]*EAmp[i]/dat->SM[1][j];
        fprintf(out,"%e %e %e %e\n", TF[i], FF[i], Aint[i], Eint[i]);
    }
    fclose(out);
     
    
    gsl_interp_accel *AAacc = gsl_interp_accel_alloc();
    gsl_spline *AAspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(AAspline, FF, AAmp, NF);
    
    gsl_interp_accel *AEacc = gsl_interp_accel_alloc();
    gsl_spline *AEspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(AEspline, FF, EAmp, NF);
    
    gsl_interp_accel *IAacc = gsl_interp_accel_alloc();
    gsl_spline *IAspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(IAspline, FF, Aint, NF);
    
    gsl_interp_accel *IEacc = gsl_interp_accel_alloc();
    gsl_spline *IEspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(IEspline, FF, Eint, NF);
    
    SNR = gsl_spline_eval_integ(IAspline, FF[0], FF[NF-1], IAacc)+gsl_spline_eval_integ(IEspline, FF[0], FF[NF-1], IEacc);
    SNR *= dat->Tobs;
    SNR = sqrt(SNR);
    printf("SNR = %e SNRsq = %e\n", SNR, SNR*SNR);
    
    /*
     
     // This section is used to get SNR versus time
    
    SNRSQA = (double*)malloc(sizeof(double)* (NF));
    SNRSQE = (double*)malloc(sizeof(double)* (NF));
    
    SNRSQA[0] = 0.0;
    SNRSQE[0] = 0.0;
    for (i=1; i< NF; i++)
    {
     SNRSQA[i] = SNRSQA[i-1]+gsl_spline_eval_integ(IAspline, FF[i-1], FF[i], IAacc);
     SNRSQE[i] = SNRSQE[i-1]+gsl_spline_eval_integ(IEspline, FF[i-1], FF[i], IEacc);
    }
    
    double SnrA, SnrE, SAT, SET, tmin, tmax;
    
    SnrA = SNRSQA[NF-1];
    SnrE = SNRSQE[NF-1];
    
    SNR = sqrt(SnrA+SnrE);
    
    printf("SNR = %e SNRsq = %e\n", SNR, SNR*SNR);
    
    //printf("%f %f %e\n", SnrA, SnrE, SNRSQA[NF-1]+SNRSQE[NF-1]);

    gsl_interp_accel *SAacc = gsl_interp_accel_alloc();
    gsl_spline *SAspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(SAspline, TF, SNRSQA, NF);
    
    gsl_interp_accel *SEacc = gsl_interp_accel_alloc();
    gsl_spline *SEspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(SEspline, TF, SNRSQE, NF);
    
    
    NW = (int)(dat->Tobs/1000.0);
    tmin = TF[0];
    tmax = TF[NF-1];
    out = fopen("SNRofT.dat","w");
    
    flag = 0;
    
    for (i=1; i< NW; i++)
    {
        t = tmin+(double)(i)*1000.0;
        SAT = SnrA;  // make sure we get to the final value
        SET = SnrE;
        if(t < tmax)
        {
            SAT = gsl_spline_eval (SAspline, t, SAacc);
            SET = gsl_spline_eval (SEspline, t, SEacc);
        }
        
        x = sqrt(SAT+SET);
        
        if(x > 10.0 && flag == 0)
        {
            flag = 1;
            y = t;
        }
        
        fprintf(out,"%e %e %e %e\n", t, x, sqrt(SAT), sqrt(SET));
    }
    fclose(out);
    
    t10 = y;
    
    if(flag == 1)
    {
    // refine estimate of when SNR = 10 is reached
    
    flag = 0;
    
    for (i=-1000; i< 1000; i++)
    {
        t = y+(double)(i);
        
        SAT = SnrA;  // make sure we get to the final value
        SET = SnrE;
        if(t < tmax)
        {
            SAT = gsl_spline_eval (SAspline, t, SAacc);
            SET = gsl_spline_eval (SEspline, t, SEacc);
        }
        
        x = sqrt(SAT+SET);
        
        if(x > 10.0 && flag == 0)
        {
            flag = 1;
            t10 = t;
        }
        
    }
 
    
    printf("SNR=10 at t = %e, tc-t = %e\n", t10, params[5]-t10);
     
    }
    
    
    free(SNRSQA);
    free(SNRSQE);
     
     */
     
     
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
    
    gsl_spline_free(IAspline);
    gsl_spline_free(IEspline);
    gsl_spline_free(AAspline);
    gsl_spline_free(AEspline);
    gsl_interp_accel_free(IAacc);
    gsl_interp_accel_free(IEacc);
    gsl_interp_accel_free(AAacc);
    gsl_interp_accel_free(AEacc);
    
    return(SNR);
    
}

void Extrinsic(double *params, double Tstart, double Tend, int NF, double *FF, double *TF, double *PF, double *AF, double *AAmp, double *EAmp, double *APhase, double *EPhase, double *kxm)
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
    
    double t, f, kdotx, Amp, Phase, RR, II;
    
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
         if(t > Tend || t < Tstart)
         {
             x = 0.0;
         }
         else
         {
          if(t < Tstart+t_tuke) x = 0.5*(1.0+cos(PI*((t-Tstart)/t_tuke-1.0)));
          if(t > Tend-t_tuke && t < Tend) x = 0.5*(1.0-cos(PI*(t-Tend)/t_tuke));
         }
        
        
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
        
        //printf("%e %e %e\n", f, AAmp[n], EAmp[n]);
        
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


void SetUp(struct Data *dat, int ll, double *params, int NFmax, int *NFS, double *FF)
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

    
    tc = params[5];    // merger time
    
    StartStop(ll, params, dat->Tstart, dat->Tend, dat->dt, &fmin, &fmax, &fr);
    
   // printf("%f %f %f\n", fmin, fr, fmax);

    
    // this can happen when tc is really small and the masses are large
    if(fmax < fmin)
    {
        fmin = 0.5*fmax;
    }
    

    dfmin = 1.0/dat->Tobs;
    //dfmax = fmax/100.0;
   // DT = 5.0e4;
    
    dfmax = fmax/1000.0;
    DT = 1.0e4;
    
    fac = DT*pow(8.0*PI, 8.0/3.0)*3.0/40.0*pow(Mc,5.0/3.0);
    
    f = fmin;
    NF = 1;
    do
    {
        df = fac*pow(f,11.0/3.0);
       // printf("%e %e %e %e\n", f, df, dfmin, dfmax);
        if(df < dfmin) df = dfmin;
        if(df > dfmax) df = dfmax;
        f += df;
        NF++;
         //printf("%d %e\n", NF, f);
    }while(f < fmax);
    
    //printf("%d %e %e\n", NF, fmin, fmax);
    
     // Need to catch is NF > NFmax
    
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
        df = 1.0/(20.0*dat->dt);
        FF[0] = dfmin;
        for (i=1; i< NF; i++) FF[i] = FF[i-1]+df;
    }
        
    //printf("NF = %d\n", NF);
    
    //for (i=0; i< NF; i++) printf("%d %e\n", i, FF[i]);

    *NFS = NF;
    
}




// This uses the smooth component of the power spectra
void FisherFast(struct Data *dat, int ll, double *params, double **Fisher)
{
    double *AAmp, *EAmp, *APhase, *EPhase, *TF, *FF, *PF, *AF;
    double *PFref, *AFref, *TFref;
    double af, fr, df, DT, fac, deltaF, f, fmin, x, t;
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double Amp, Phase, tf, fguess;
    double m1_SI, m2_SI, distance, tc, phic;
    double *paramsP, *paramsM;
    int i, j, k, ii, NF, NW, lx;
    double A, P;
    double px, fnew, kxm;
    double *Aint, *Eint, *SADS, *SEDS;
    double *AphaseP, *AphaseM, *AampP, *AampM;
    double *EphaseP, *EphaseM, *EampP, *EampM;
    double **DphaseA, **DampA;
    double **DphaseE, **DampE;
    double *epsilon;
    double Tstart, Tend, Tobs;
    
    int NFmax = 100000;
    
    int kmin, kmax, flag;
    
    Tstart = dat->Tstart;
    Tend = dat->Tend;
    Tobs = dat->Tobs;
    
    epsilon = (double*)malloc(sizeof(double)* (NP));
    paramsP = (double*)malloc(sizeof(double)* (NP));
    paramsM = (double*)malloc(sizeof(double)* (NP));
    
    FF = (double*)malloc(sizeof(double)* (NFmax));
    
    SetUp(dat, ll, params, NFmax, &NF, FF);
    
    AF = (double*)malloc(sizeof(double)* (NF));
    PF = (double*)malloc(sizeof(double)* (NF));
    TF = (double*)malloc(sizeof(double)* (NF));
    
    TFref = (double*)malloc(sizeof(double)* (NF));
    PFref = (double*)malloc(sizeof(double)* (NF));
    AFref = (double*)malloc(sizeof(double)* (NF));
    
    Intrinsic(ll, params, Tobs, NF, FF, TFref, PFref, AFref);
    
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
    // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
    
    // will take log derivatives for m1, m2, DL
    
      epsilon[0] = 1.0e-6;
      epsilon[1] = 1.0e-6;
      epsilon[2] = 1.0e-4;
      epsilon[3] = 1.0e-4;
      epsilon[4] = 1.0e-5;
      epsilon[5] = 1.0e-3;
      epsilon[6] = 1.0e-5;
      epsilon[7] = 1.0e-5;
      epsilon[8] = 1.0e-5;
      epsilon[9] = 1.0e-5;
      epsilon[10] = 1.0e-5;
    
  
    
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
    SADS = (double*)malloc(sizeof(double)* (NF));
    SEDS = (double*)malloc(sizeof(double)* (NF));
    
   
    
    // reference amplitude and phase
    Extrinsic(params, Tstart, Tend, NF, FF, TFref, PFref, AFref, AAmp, EAmp, APhase, EPhase, &kxm);
    
    
    // need to resample SAS and SES to get SADS and SEDS
    
    for (k=0; k< NF; k++)
     {
         i = (int)(FF[k]*Tobs);
         //printf("%d %d\n", k, i);
         SADS[k] = dat->SM[0][i];
         SEDS[k] = dat->SM[1][i];
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

            Intrinsic(ll,paramsP, Tobs, NF, FF, TF, PF, AF);  // intrinsic
            Extrinsic(paramsP, Tstart, Tend, NF, FF, TF, PF, AF, AampP, EampP, AphaseP, EphaseP, &kxm); // extrinsic
                
            Intrinsic(ll,paramsM, Tobs, NF, FF, TF, PF, AF); // intrinsic
            Extrinsic(paramsM, Tstart, Tend, NF, FF, TF, PF, AF, AampM, EampM, AphaseM, EphaseM, &kxm); // extrinsic

            
            for (k=0; k< NF; k++)
            {
                DphaseA[i][k] = (AphaseP[k]-AphaseM[k]);
                DphaseE[i][k] = (EphaseP[k]-EphaseM[k]);
                DampA[i][k] = (AampP[k]-AampM[k]);
                DampE[i][k] = (EampP[k]-EampM[k]);
            }
        
    }

    

      // Put in coallesence phase manually. Note this is orbital phase, so times two in GW phase
        i = 4;
        x = -4.0*epsilon[i];
        for (k=0; k< NF; k++)
        {
            DphaseA[i][k] = x;
            DphaseE[i][k] = x;
            DampA[i][k] = 0.0;
            DampE[i][k] = 0.0;
        }
    
    
       // tc is not included in the phase subroutine
        i = 5;
        for (k=0; k< NF; k++)
        {
            x = -4.0*PI*FF[k]*epsilon[i];
            DphaseA[i][k] = x;
            DphaseE[i][k] = x;
            DampA[i][k] = 0.0;
            DampE[i][k] = 0.0;
        }
        
         // ln DL derivative
           i = 6;
           for (k=0; k< NF; k++)
           {
               DphaseA[i][k] = 0;
               DphaseE[i][k] = 0;
               DampA[i][k] = -2.0*epsilon[i]*AAmp[k];
               DampE[i][k] = -2.0*epsilon[i]*EAmp[k];
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

        Extrinsic(paramsP, Tstart, Tend, NF, FF, TFref, PFref, AFref, AampP, EampP, AphaseP, EphaseP, &kxm); // extrinsic
        Extrinsic(paramsM, Tstart, Tend, NF, FF, TFref, PFref, AFref, AampM, EampM, AphaseM, EphaseM, &kxm); // extrinsic
        
        for (k=0; k< NF; k++)
        {
            DphaseA[i][k] = (AphaseP[k]-AphaseM[k]);
            DphaseE[i][k] = (EphaseP[k]-EphaseM[k]);
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
  
            for (k=0; k< NF; k++)
            {
             Aint[k] = 4.0*(DampA[i][k]*DampA[j][k]+AAmp[k]*AAmp[k]*DphaseA[i][k]*DphaseA[j][k])/SADS[k];
             Eint[k] = 4.0*(DampE[i][k]*DampE[j][k]+EAmp[k]*EAmp[k]*DphaseE[i][k]*DphaseE[j][k])/SEDS[k];
            }
            
            gsl_spline_init(IAspline, FF, Aint, NF);
            gsl_spline_init(IEspline, FF, Eint, NF);
            
            Fisher[i][j] = Tobs*(gsl_spline_eval_integ(IAspline, FF[0], FF[NF-1], IAacc)+gsl_spline_eval_integ(IEspline, FF[0], FF[NF-1], IEacc));
            Fisher[i][j] /= (4.0*epsilon[i]*epsilon[j]);
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
    free(SADS);
    free(SEDS);
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

void FisherHet(struct Data *dat, struct Het *het, int ll, double *params, double **Fisher)
{
    double af, fr, df, DT, fac, deltaF, f, fmin, x, t;
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double Amp, Phase, tf, fguess;
    double m1_SI, m2_SI, distance, tc, phic;
    double *paramsP, *paramsM;
    int i, j, k, ii, jj, kk, NF, NW, lx;
    double A, P;
    double px, fnew, kxm;
    double *epsilon;
    double Tstart, Tend, Tobs;

    int kmin, kmax, flag, Nch;
    
    int M, J, NR, id;
    double HH;
    double **cc, *uu, *vv;
    double **amp, **phase;
    double ***ampP, ***phaseP;
    double ***ampM, ***phaseM;
    
    Nch = het->Nch;
    M = het->M;
    J = het->J;
    NR = het->NR;
    cc = double_matrix(Nch,M);
    uu = double_vector(J+1);
    vv = double_vector(J+1);
    
    Tstart = dat->Tstart;
    Tend = dat->Tend;
    Tobs = dat->Tobs;
    
    epsilon = (double*)malloc(sizeof(double)* (NP));
    paramsP = (double*)malloc(sizeof(double)* (NP));
    paramsM = (double*)malloc(sizeof(double)* (NP));
    
    amp = double_matrix(Nch,M);
    phase = double_matrix(Nch,M);
    ampP = double_tensor(NP,Nch,M);
    phaseP = double_tensor(NP,Nch,M);
    ampM = double_tensor(NP,Nch,M);
    phaseM = double_tensor(NP,Nch,M);
    
    // Reference phase and amplitude
    fullphaseamp(dat, ll, M, params, het->freq, amp[0], amp[1], phase[0], phase[1]);
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
    // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
    
    // will take log derivatives for m1, m2, DL
    
      epsilon[0] = 1.0e-6;
      epsilon[1] = 1.0e-6;
      epsilon[2] = 1.0e-4;
      epsilon[3] = 1.0e-4;
      epsilon[4] = 1.0e-5;
      epsilon[5] = 10.0;
      epsilon[6] = 1.0e-4;
      epsilon[7] = 1.0e-5;
      epsilon[8] = 1.0e-5;
      epsilon[9] = 1.0e-5;
      epsilon[10] = 1.0e-5;
    
    for (i=0; i< NP; i++)
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
        
           fullphaseamp(dat, ll, M, paramsP, het->freq, ampP[i][0], ampP[i][1], phaseP[i][0], phaseP[i][1]);
           fullphaseamp(dat, ll, M, paramsM, het->freq, ampM[i][0], ampM[i][1], phaseM[i][0], phaseM[i][1]);
           
          // store the central difference in the Plus arrays
          for (id = 0 ; id < Nch ; id++)  // loop over channels
          {
          for (k = 0 ; k < M ; k++)
          {
              ampP[i][id][k] -= ampM[i][id][k];
              phaseP[i][id][k] -= phaseM[i][id][k];
              ampP[i][id][k] /= (2.0*epsilon[i]);
              phaseP[i][id][k] /= (2.0*epsilon[i]);
          }
          }
        
    }


    for (i = 0; i < NP; i++)
    {
        for (j = 0; j < NP; j++) Fisher[i][j] = 0.0;
    }
    
        for (i = 0 ; i < NP ; i++)
        {
            for (j = i ; j < NP ; j++)
            {
                
                for (id = 0 ; id < Nch ; id++)  // loop over detectors
                {
                
                HH = 0.0;
                
                // these are the slow terms in the Fisher matrix sum
                for(jj = 0; jj < M; jj++)
                {
                 cc[id][jj] = 0.0;
                 if(het->amp[id][jj] > 0.0)
                  {
                   cc[id][jj] = 4.0*(ampP[i][id][jj]*ampP[j][id][jj]
                                  +amp[id][jj]*amp[id][jj]*phaseP[i][id][jj]*phaseP[j][id][jj]);
                  cc[id][jj] /= (het->amp[id][jj]*het->amp[id][jj]);
                  }
                }
         
            for(ii = 0; ii < NR; ii++)
             {
                                                     
                for(jj = 0; jj <= J; jj++)
                {
                  vv[jj] = cc[id][ii*J+jj];
                  uu[jj] = 0.0;
                }
                                                     
                // get Legendre coefficients for slow term
               for(jj = 0; jj <= J; jj++)
                {
                for(kk = 0; kk <= J; kk++)
                 {
                  uu[jj] += het->IP[ii][jj][kk]*vv[kk];
                 }
                }
                                                     
                for(jj = 0; jj <= J; jj++) HH += uu[jj]*het->SL[ii][id][jj];
                                                     
                if(ii < NR-2) HH -= het->aa[id][ii+1]*vv[J];  // correction for overcount
                                      
             }
                Fisher[i][j] += HH;
            }
        }
        
    }


    for (i=0; i< NP; i++)
    {
        for (j=i+1; j< NP; j++)
        {
            Fisher[j][i] =  Fisher[i][j];
        }
    }
    

    free_double_matrix(cc,Nch);
    free_double_vector(uu);
    free_double_vector(vv);
 
    free(epsilon);
    free(paramsP);
    free(paramsM);
    
    free_double_matrix(amp,Nch);
    free_double_matrix(phase,Nch);
    free_double_tensor(ampP,NP,Nch);
    free_double_tensor(phaseP,NP,Nch);
    free_double_tensor(ampM,NP,Nch);
    free_double_tensor(phaseM,NP,Nch);
    
    
    
}

void FisherSubHet(struct Data *dat, struct Het *het, int ll, int *pmap, double *params, double **Fisher)
{
    double af, fr, df, DT, fac, deltaF, f, fmin, x, t;
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double Amp, Phase, tf, fguess;
    double m1_SI, m2_SI, distance, tc, phic;
    double *paramsP, *paramsM;
    int i, j, k, ii, jj, kk, NF, NW, lx;
    double A, P;
    double px, fnew, kxm;
    double *epsilon;
    double Tstart, Tend, Tobs;

    int kmin, kmax, flag, Nch;
    
    int M, J, NR, id;
    double HH;
    double **cc, *uu, *vv;
    double **amp, **phase;
    double ***ampP, ***phaseP;
    double ***ampM, ***phaseM;
    
    Nch = het->Nch;
    M = het->M;
    J = het->J;
    NR = het->NR;
    cc = double_matrix(Nch,M);
    uu = double_vector(J+1);
    vv = double_vector(J+1);
    
    Tstart = dat->Tstart;
    Tend = dat->Tend;
    Tobs = dat->Tobs;
    
    epsilon = (double*)malloc(sizeof(double)* (NP));
    paramsP = (double*)malloc(sizeof(double)* (NP));
    paramsM = (double*)malloc(sizeof(double)* (NP));
    
    amp = double_matrix(Nch,M);
    phase = double_matrix(Nch,M);
    ampP = double_tensor(NP,Nch,M);
    phaseP = double_tensor(NP,Nch,M);
    ampM = double_tensor(NP,Nch,M);
    phaseM = double_tensor(NP,Nch,M);
    
    // Reference phase and amplitude
    fullphaseamp(dat, ll, M, params, het->freq, amp[0], amp[1], phase[0], phase[1]);
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
    // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
    
    // will take log derivatives for m1, m2, DL
    
      epsilon[0] = 1.0e-6;
      epsilon[1] = 1.0e-6;
      epsilon[2] = 1.0e-4;
      epsilon[3] = 1.0e-4;
      epsilon[4] = 1.0e-5;
      epsilon[5] = 10.0;
      epsilon[6] = 1.0e-4;
      epsilon[7] = 1.0e-5;
      epsilon[8] = 1.0e-5;
      epsilon[9] = 1.0e-5;
      epsilon[10] = 1.0e-5;
    
    for (i=0; i< NP; i++)
    {
        if(pmap[i] > -1)
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
        
           fullphaseamp(dat, ll, M, paramsP, het->freq, ampP[i][0], ampP[i][1], phaseP[i][0], phaseP[i][1]);
           fullphaseamp(dat, ll, M, paramsM, het->freq, ampM[i][0], ampM[i][1], phaseM[i][0], phaseM[i][1]);
           
          // store the central difference in the Plus arrays
          for (id = 0 ; id < Nch ; id++)  // loop over channels
          {
          for (k = 0 ; k < M ; k++)
          {
              ampP[i][id][k] -= ampM[i][id][k];
              phaseP[i][id][k] -= phaseM[i][id][k];
              ampP[i][id][k] /= (2.0*epsilon[i]);
              phaseP[i][id][k] /= (2.0*epsilon[i]);
          }
          }
             
         }
        
    }


    for (i = 0; i < NP; i++)
    {
        if(pmap[i] > -1)
            {
        for (j = 0; j < NP; j++) if(pmap[j] > -1) Fisher[pmap[i]][pmap[j]] = 0.0;
            }
    }
    
        for (i = 0 ; i < NP ; i++)
        {
            if(pmap[i] > -1)
            {
            for (j = i ; j < NP ; j++)
            {
                if(pmap[j] > -1)
                {
                    
                for (id = 0 ; id < Nch ; id++)  // loop over detectors
                {
                HH = 0.0;
                // these are the slow terms in the Fisher matrix sum
                   
                for(jj = 0; jj < M; jj++)
                {
                    cc[id][jj] = 0.0;
                    if(het->amp[id][jj] > 0.0)
                    {
                    cc[id][jj] = 4.0*(ampP[i][id][jj]*ampP[j][id][jj]
                                +amp[id][jj]*amp[id][jj]*phaseP[i][id][jj]*phaseP[j][id][jj]);
                    cc[id][jj] /= (het->amp[id][jj]*het->amp[id][jj]);
                    }
                }
         
            for(ii = 0; ii < NR; ii++)
             {
                                                     
                for(jj = 0; jj <= J; jj++)
                {
                  vv[jj] = cc[id][ii*J+jj];
                  uu[jj] = 0.0;
                }
                                                     
                // get Legendre coefficients for slow term
               for(jj = 0; jj <= J; jj++)
                {
                for(kk = 0; kk <= J; kk++)
                 {
                  uu[jj] += het->IP[ii][jj][kk]*vv[kk];
                 }
                }
                for(jj = 0; jj <= J; jj++) HH += uu[jj]*het->SL[ii][id][jj];
                if(ii < NR-2) HH -= het->aa[id][ii+1]*vv[J];  // correction for overcount
                } // loop over frequency segments
                Fisher[pmap[i]][pmap[j]] += HH;
                } // id loop
                    
            } // fi pmap[j] > -1
          } // j loop
        }  // fi pmap[i] > -1
      }  // i loop

    for (i=0; i< NP; i++)
       {
           if(pmap[i] > -1)
           {
               
           for (j=i+1; j< NP; j++)
           {
               if(pmap[j] > -1)
               {
               Fisher[pmap[j]][pmap[i]] =  Fisher[pmap[i]][pmap[j]];
               }
           }
               
           }
       }

    free_double_matrix(cc,Nch);
    free_double_vector(uu);
    free_double_vector(vv);
 
    free(epsilon);
    free(paramsP);
    free(paramsM);
    
    free_double_matrix(amp,Nch);
    free_double_matrix(phase,Nch);
    free_double_tensor(ampP,NP,Nch);
    free_double_tensor(phaseP,NP,Nch);
    free_double_tensor(ampM,NP,Nch);
    free_double_tensor(phaseM,NP,Nch);
    
    
    
}


void FisherSub(struct Data *dat, int ll, int *pmap, double *params, double **Fisher)
{
    double *AAmp, *EAmp, *APhase, *EPhase, *TF, *FF, *PF, *AF;
    double *PFref, *AFref, *TFref;
    double af, fr, df, DT, fac, deltaF, f, fmin, x, t;
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double Amp, Phase, tf, fguess;
    double m1_SI, m2_SI, distance, tc, phic;
    double *paramsP, *paramsM;
    int i, j, k, ii, NF, NW;
    double A, P;
    double px, fnew, kxm;
    double *Aint, *Eint, *SADS, *SEDS;
    double *AphaseP, *AphaseM, *AampP, *AampM;
    double *EphaseP, *EphaseM, *EampP, *EampM;
    double **DphaseA, **DampA;
    double **DphaseE, **DampE;
    double *epsilon;
    double Tstart, Tend, Tobs;
    
    int NFmax = 100000;
    
    int kmin, kmax, flag;
    
    Tstart = dat->Tstart;
    Tend = dat->Tend;
    Tobs = dat->Tobs;

    epsilon = (double*)malloc(sizeof(double)* (NP));
    paramsP = (double*)malloc(sizeof(double)* (NP));
    paramsM = (double*)malloc(sizeof(double)* (NP));
    
    FF = (double*)malloc(sizeof(double)* (NFmax));
    
    SetUp(dat, ll, params, NFmax, &NF, FF);
    
    AF = (double*)malloc(sizeof(double)* (NF));
    PF = (double*)malloc(sizeof(double)* (NF));
    TF = (double*)malloc(sizeof(double)* (NF));
    
    TFref = (double*)malloc(sizeof(double)* (NF));
    PFref = (double*)malloc(sizeof(double)* (NF));
    AFref = (double*)malloc(sizeof(double)* (NF));
    
    Intrinsic(ll, params, Tobs, NF, FF, TFref, PFref, AFref);
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
    // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
    
    // will take log derivatives for m1, m2, DL
    
    epsilon[0] = 1.0e-6;
    epsilon[1] = 1.0e-6;
    epsilon[2] = 1.0e-4;
    epsilon[3] = 1.0e-4;
    epsilon[4] = 1.0e-5;
    epsilon[5] = 1.0e-3;
    epsilon[6] = 1.0e-4;
    epsilon[7] = 1.0e-5;
    epsilon[8] = 1.0e-5;
    epsilon[9] = 1.0e-5;
    epsilon[10] = 1.0e-5;
    
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
    SADS = (double*)malloc(sizeof(double)* (NF));
    SEDS = (double*)malloc(sizeof(double)* (NF));
   
    
    // reference amplitude and phase
    Extrinsic(params, Tstart, Tend, NF, FF, TFref, PFref, AFref, AAmp, EAmp, APhase, EPhase, &kxm);
    
    // need to resample SAS and SES to get SADS and SEDS
       
       for (k=0; k< NF; k++)
        {
            i = (int)(FF[k]*Tobs);
            //printf("%d %d\n", k, i);
            SADS[k] = dat->SM[0][i];
            SEDS[k] = dat->SM[1][i];
        }
    
    // masses and spins
    for (i=0; i< 4; i++)
    {
        if(pmap[i] > -1)
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
        
        
        
        Intrinsic(ll,paramsP, Tobs, NF, FF, TF, PF, AF);  // intrinsic
        Extrinsic(paramsP, Tstart, Tend, NF, FF, TF, PF, AF, AampP, EampP, AphaseP, EphaseP, &kxm); // extrinsic
        
        Intrinsic(ll,paramsM, Tobs, NF, FF, TF, PF, AF); // intrinsic
        Extrinsic(paramsM, Tstart, Tend, NF, FF, TF, PF, AF, AampM, EampM, AphaseM, EphaseM, &kxm); // extrinsic
        
        
        for (k=0; k< NF; k++)
        {
            DphaseA[i][k] = (AphaseP[k]-AphaseM[k]);
            DphaseE[i][k] = (EphaseP[k]-EphaseM[k]);
            DampA[i][k] = (AampP[k]-AampM[k]);
            DampE[i][k] = (EampP[k]-EampM[k]);
        }
            
       }
        
    }
    
    
    
    for (i=4; i< 7; i++)
    {
        if(pmap[i] > -1)
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
    }
    
    for (i=7; i< NP; i++)
    {
        if(pmap[i] > -1)
        {
            
        for (j=0; j< NP; j++)
        {
            paramsP[j] = params[j];
            paramsM[j] = params[j];
        }
        
        paramsP[i] += epsilon[i];
        paramsM[i] -= epsilon[i];
        
        Extrinsic(paramsP, Tstart, Tend, NF, FF, TFref, PFref, AFref, AampP, EampP, AphaseP, EphaseP, &kxm); // extrinsic
        Extrinsic(paramsM, Tstart, Tend, NF, FF, TFref, PFref, AFref, AampM, EampM, AphaseM, EphaseM, &kxm); // extrinsic
        
        for (k=0; k< NF; k++)
        {
            DphaseA[i][k] = (AphaseP[k]-AphaseM[k]);
            DphaseE[i][k] = (EphaseP[k]-EphaseM[k]);
            DampA[i][k] = (AampP[k]-AampM[k]);
            DampE[i][k] = (EampP[k]-EampM[k]);
        }
            
        }
        
    }
    
    
    gsl_interp_accel *IAacc = gsl_interp_accel_alloc();
    gsl_spline *IAspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    
    
    gsl_interp_accel *IEacc = gsl_interp_accel_alloc();
    gsl_spline *IEspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    
    
    
    for (i=0; i< NP; i++)
    {
        if(pmap[i] > -1)
        {
            
        for (j=i; j< NP; j++)
        {
            if(pmap[j] > -1)
            {
                
            for (k=0; k< NF; k++)
            {
                Aint[k] = 4.0*(DampA[i][k]*DampA[j][k]+AAmp[k]*AAmp[k]*DphaseA[i][k]*DphaseA[j][k])/SADS[k];
                Eint[k] = 4.0*(DampE[i][k]*DampE[j][k]+EAmp[k]*EAmp[k]*DphaseE[i][k]*DphaseE[j][k])/SEDS[k];
            }
            
            gsl_spline_init(IAspline, FF, Aint, NF);
            gsl_spline_init(IEspline, FF, Eint, NF);
            
            Fisher[pmap[i]][pmap[j]] = Tobs*(gsl_spline_eval_integ(IAspline, FF[0], FF[NF-1], IAacc)+gsl_spline_eval_integ(IEspline, FF[0], FF[NF-1], IEacc));
            Fisher[pmap[i]][pmap[j]] /= (4.0*epsilon[i]*epsilon[j]);
                
            }
        }
            
        }
        
    }
    
    for (i=0; i< NP; i++)
    {
        if(pmap[i] > -1)
        {
            
        for (j=i+1; j< NP; j++)
        {
            if(pmap[j] > -1)
            {
            Fisher[pmap[j]][pmap[i]] =  Fisher[pmap[i]][pmap[j]];
            }
        }
            
        }
    }
    
    free_double_matrix(DphaseA,NP);
    free_double_matrix(DampA,NP);
    free_double_matrix(DphaseE,NP);
    free_double_matrix(DampE,NP);
    
    free(epsilon);
    free(SADS);
    free(SEDS);
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

void FisherDirect(struct Data *dat, int ll, double *params, double **Fisher)
{
    double af, fr, df, DT, fac, deltaF, f, fmin, x, t;
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double Amp, Phase, tf, fguess;
    double m1_SI, m2_SI, distance, tc, phic;
    double *paramsP, *paramsM;
    int i, j, k, ii, NF, NW, lx, N;
    double A, P;
    double px, fnew, kxm;
    double *epsilon;
    double **AH, **AP, **AM;
    double **EH, **EP, **EM;
    
    int kmin, kmax, flag;
    
    N = dat->N;

    
    epsilon = (double*)malloc(sizeof(double)* (NP));
    paramsP = (double*)malloc(sizeof(double)* (NP));
    paramsM = (double*)malloc(sizeof(double)* (NP));
    
    AH = double_matrix(NP,N);
    AP = double_matrix(NP,N);
    AM = double_matrix(NP,N);
    EH = double_matrix(NP,N);
    EP = double_matrix(NP,N);
    EM = double_matrix(NP,N);
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
    // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
    
    // will take log derivatives for m1, m2, DL
    
    epsilon[0] = 1.0e-6;
    epsilon[1] = 1.0e-6;
    epsilon[2] = 1.0e-4;
    epsilon[3] = 1.0e-4;
    epsilon[4] = 1.0e-4;
    epsilon[5] = 1.0e-2;
    epsilon[6] = 1.0e-4;
    epsilon[7] = 1.0e-4;
    epsilon[8] = 1.0e-4;
    epsilon[9] = 1.0e-4;
    epsilon[10] = 1.0e-4;
    
    for (i=0; i< NP; i++)
    {
        
        for (j=0; j< NP; j++)
        {
            paramsP[j] = params[j];
            paramsM[j] = params[j];
        }
        
        
        if(ll==0)
        {
            if(i==0 || i==1 || i==6) // log derivatives
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

        ResponseFreq(dat, ll, paramsP, AP[i], EP[i]);
        ResponseFreq(dat, ll, paramsM, AM[i], EM[i]);

        for (j=0; j< N; j++) AH[i][j] = (AP[i][j]-AM[i][j])/(2.0*epsilon[i]);
        for (j=0; j< N; j++) EH[i][j] = (EP[i][j]-EM[i][j])/(2.0*epsilon[i]);
        
    }

       free_double_matrix(AP,NP);
       free_double_matrix(AM,NP);
       free_double_matrix(EP,NP);
       free_double_matrix(EM,NP);
    
       free(epsilon);
       free(paramsP);
       free(paramsM);

   
    for (i=0; i< NP; i++)
    {
        for (j=i; j< NP; j++)
        {
  
            Fisher[i][j] =  fourier_nwip(AH[i], AH[j], dat->SN[0], N) + fourier_nwip(EH[i], EH[j], dat->SN[1], N);

        }
        
    }
    
    for (i=0; i< NP; i++)
    {
        for (j=i+1; j< NP; j++)
        {
            Fisher[j][i] =  Fisher[i][j];
        }
    }
    
    free_double_matrix(AH,NP);
    free_double_matrix(EH,NP);
    
    
}

double chisq_het(struct Data *dat, struct Het *het, int ll, double *params, double **ampR, double **phaseR)
{
    double **amp, **phase;
    double *uu, *vv, **cc;
    int Nch, M, J, NR;
    double csq;
    int i, id, ii, jj, kk;
    
    Nch = het->Nch;
    M = het->M;
    J = het->J;
    NR = het->NR;
    cc = double_matrix(Nch,M);
    uu = double_vector(J+1);
    vv = double_vector(J+1);
    amp = double_matrix(Nch,M);
    phase = double_matrix(Nch,M);
    
    //  phase and amplitude
    fullphaseamp(dat, ll, M, params, het->freq, amp[0], amp[1], phase[0], phase[1]);
    
    
    csq = 0.0;
    
           for (id = 0 ; id < Nch ; id++)  // loop over detectors
           {
           
           // these are the slow terms in the chi-squared sum
           for(jj = 0; jj < M; jj++)
           {
            cc[id][jj] = 0.0;
            if(het->amp[id][jj] > 0.0)
             {
              cc[id][jj] = 4.0*(ampR[id][jj]*ampR[id][jj]+amp[id][jj]*amp[id][jj]-2.0*amp[id][jj]*ampR[id][jj]*cos(phaseR[id][jj]-phase[id][jj]));
              cc[id][jj] /= (het->amp[id][jj]*het->amp[id][jj]);
             }
           }
    
       for(ii = 0; ii < NR; ii++)
        {
                                                
           for(jj = 0; jj <= J; jj++)
           {
             vv[jj] = cc[id][ii*J+jj];
             uu[jj] = 0.0;
           }
                                                
           // get Legendre coefficients for slow term
          for(jj = 0; jj <= J; jj++)
           {
           for(kk = 0; kk <= J; kk++)
            {
             uu[jj] += het->IP[ii][jj][kk]*vv[kk];
            }
           }
                                                
           for(jj = 0; jj <= J; jj++) csq += uu[jj]*het->SL[ii][id][jj];
                                                
           if(ii < NR-2) csq -= het->aa[id][ii+1]*vv[J];  // correction for overcount
                                 
        }
               
        }

    
     free_double_matrix(cc,Nch);
     free_double_vector(uu);
     free_double_vector(vv);
     free_double_matrix(amp,Nch);
     free_double_matrix(phase,Nch);
    
    return(csq);
    
}

double chisq(struct Data *dat, int ll, double *params, double *AR, double *ER)
{
    double *AS, *ES;
    double HH, HD;
    double csq;
    int i;
    
    AS = (double*)malloc(sizeof(double)* (dat->N));
    ES = (double*)malloc(sizeof(double)* (dat->N));
    
    ResponseFast(dat, ll, params, AS, ES);
    
    for (i=0; i< dat->N; i++)
       {
           AS[i] -= AR[i];
           ES[i] -= ER[i];
       }
    
    csq = (fourier_nwip(AS, AS, dat->SN[0], dat->N)+fourier_nwip(ES, ES, dat->SN[1], dat->N));
    
    
    free(AS);
    free(ES);
    
    return(csq);
    
}

void efix(struct Data *dat, struct Het *het, int hr, int ll, double *params, double *min, double *max, double *eval, double **evec, double zs)
{
    int i, j, k, flag;
    double alpha0, x, z, z0, alpha, alphanew;
    double dz, alphax;
    double leta, eta;
    double zmx, zmn;
    double dzmin, alpham, zm;
    double *px;
    double **ampR, **phaseR;
    double *AR, *ER;
    
    if(hr == 1)  // using heterodyne
    {
    ampR = double_matrix(het->Nch,het->M);
    phaseR = double_matrix(het->Nch,het->M);
    fullphaseamp(dat, ll, het->M, params, het->freq, ampR[0], ampR[1], phaseR[0], phaseR[1]);
    }
    else
    {
    AR = (double*)malloc(sizeof(double)* (dat->N));
    ER = (double*)malloc(sizeof(double)* (dat->N));
    ResponseFast(dat, ll, params, AR, ER);
    }
    
    // [0] ln(Mass1)  [1] ln(Mass2)  [2] Spin1 [3] Spin2 [4] phic [5] tc [6] ln(distance)
     // [7] EclipticCoLatitude, [8] EclipticLongitude  [9] polarization, [10] inclination
   
    zmx = zs*2.0;
    zmn = zs/2.0;
    
    px = double_vector(NP);

   for (i = 0; i < NP; ++i)
     {
         
        dzmin = 1.0e20;
        alpham = 1.0;
         
      alpha0 = 1.0;
         
        
      k = 0;
      do
      {
      for (j = 0; j < NP; ++j) px[j] = params[j] + alpha0*eval[i]*evec[i][j];
      
      if(ll == 2)
          {
        leta = (5.0/3.0)*(px[0]-px[1]);
        eta = exp(leta);
         if(eta > 0.25)
         {
         for (j = 0; j < NP; ++j) px[j] = params[j] - alpha0*eval[i]*evec[i][j];
         }
          }
          
          // re-map angular parameters to their proper range
          
          x = px[4]/PI;
          px[4] = (x-floor(x))*PI;
          if(x < 0.0) px[4] += 1.0;
          x = px[8]/TPI;
          px[8] = (x-floor(x))*TPI;
          if(x < 0.0) px[8] += 1.0;
          x = px[9]/PI;
          px[9] = (x-floor(x))*PI;
          if(x < 0.0) px[9] += 1.0;
          
          while(px[4] > PI)  px[4] -= PI;
          while(px[4] < 0.0)  px[4] += PI;
          while(px[8] > TPI)  px[8] -= TPI;
          while(px[8] < 0.0)  px[8] += TPI;
          while(px[9] > PI)   px[9] -= PI;
          while(px[9] < 0.0)  px[9] += PI;
          
         // for (j = 0; j < NP; ++j) printf("%d %e %e %e %e\n", j, px[j], params[j], min[j], max[j]);
          
          
          for (j = 0; j < NP; ++j) if(px[j] > max[j]) px[j] = max[j];
          for (j = 0; j < NP; ++j) if(px[j] < min[j]) px[j] = min[j];
      
          if(ll == 2)
          {
          leta = (5.0/3.0)*(px[0]-px[1]);
          eta = exp(leta);
          if(eta > 0.25) px[0] = px[1] + 3.0/5.0*log(0.2499);
          if(eta < etamin) px[0] = px[1] + 3.0/5.0*log(etamin);
          }
          
          if(hr == 1)  // using heterodyne
          {
          z0 = chisq_het(dat, het, ll, px, ampR, phaseR);
          }
          else
          {
          z0 = chisq(dat, ll, px, AR, ER);
          }
          
         //printf("B %f %f\n", alpha0, z0);
          
          dz = fabs(z0-zs);
          if(dz < dzmin)
          {
              dzmin = dz;
              alpham = alpha0;
              zm = z0;
          }
          
          if(z0 < zmn) alpha0 *= 2.0;
          if(z0 > zmx) alpha0 /= 1.9;
          
          k++;
          
      } while ( (z0 > zmx || z0 < zmn) && k < 15 && fabs(alpha0 < 1.0e4));
      
      
      alpha = alpha0*1.1;
          
        k = 0;
        do
        {
        for (j = 0; j < NP; ++j) px[j] = params[j] + alpha*eval[i]*evec[i][j];
         
        if(ll == 2)
        {
        leta = (5.0/3.0)*(px[0]-px[1]);
        eta = exp(leta);
        if(eta > 0.25)
        {
          for (j = 0; j < NP; ++j) px[j] = params[j] - alpha*eval[i]*evec[i][j];
        }
        }
            
        // re-map angular parameters to their proper range
            x = px[4]/PI;
            px[4] = (x-floor(x))*PI;
            if(x < 0.0) px[4] += 1.0;
            x = px[8]/TPI;
            px[8] = (x-floor(x))*TPI;
            if(x < 0.0) px[8] += 1.0;
            x = px[9]/PI;
            px[9] = (x-floor(x))*PI;
            if(x < 0.0) px[9] += 1.0;
            
        while(px[4] > PI)  px[4] -= PI;
        while(px[4] < 0.0)  px[4] += PI;
        while(px[8] > TPI)  px[8] -= TPI;
        while(px[8] < 0.0)  px[8] += TPI;
        while(px[9] > PI)   px[9] -= PI;
        while(px[9] < 0.0)  px[9] += PI;
            
        for (j = 0; j < NP; ++j) if(px[j] > max[j]) px[j] = max[j];
        for (j = 0; j < NP; ++j) if(px[j] < min[j]) px[j] = min[j];
            
            if(ll == 2)
            {
            leta = (5.0/3.0)*(px[0]-px[1]);
            eta = exp(leta);
            if(eta > 0.25) px[0] = px[1] + 3.0/5.0*log(0.2499);
            if(eta < etamin) px[0] = px[1] + 3.0/5.0*log(etamin);
            }
            
           
            
            if(hr == 1)  // using heterodyne
            {
            x = chisq_het(dat, het, ll, px, ampR, phaseR);
            }
            else
            {
              z = chisq(dat, ll, px, AR, ER);
            }
            
            //printf("R %f %f\n", alpha, z);
            
            if(alpha > alpha0)
            {
            alphanew = alpha0 +(zs-z0)/(z-z0)*(alpha-alpha0);
            }
            else
            {
             alphanew = alpha +(zs-z)/(z0-z)*(alpha0-alpha);
            }
            
            dz = fabs(z-zs);
            if(dz < dzmin)
            {
              dzmin = dz;
              alpham = alpha;
              zm = z;
            }
            
            z0 = z;
            alpha0 = alpha;
            alpha = alphanew;
 
             k++;
             
        } while (fabs(z-zs) > 0.2 && k < 10 && fabs(alpha) < 1.0e4);
         
          // printf("F %f %f\n", alpham, zm);
         
         // printf("\n");
         
         // kill any failed direction
         if(dzmin/zs > 1.0) alpham = 0.0;
         
         
         eval[i] *= alpham;
      
     }
    
    free_double_vector(px);
    
    if(hr == 1)  // using heterodyne
    {
    free_double_matrix(ampR,het->Nch);
    free_double_matrix(phaseR,het->Nch);
    }
    else
    {
     free(AR);
     free(ER);
    }
    
}

void instrument_noise(double f, double *SAE)
{
    //Power spectral density of the detector noise and transfer frequency
    double Sn, red, confusion_noise;
    double Sloc, fonfs;
    double f1, f2;
    double A1, A2, slope, LC;
    double Sps = 2.25e-22;
    double Sacc = 9.0e-30;
    
    fonfs = f/fstar;
    
    //LC = 16.0*fonfs*fonfs;
    
    // To match the LDC power spectra I have to divide by 10. No idea why...
    LC = 1.60*fonfs*fonfs;
    
    red = 16.0*((1.0e-4/f)*(1.0e-4/f));
    // red = 0.0;
    
    // Calculate the power spectral density of the detector noise at the given frequency
    
    *SAE = LC*16.0/3.0*pow(sin(fonfs),2.0)*( (2.0+cos(fonfs))*(Sps) + 2.0*(3.0+2.0*cos(fonfs)+cos(2.0*fonfs))*(Sacc/pow(2.0*PI*f,4.0)*(1.0+red)) ) / pow(2.0*Larm,2.0);
    
   // *SXYZ = LC*4.0*pow(sin(fonfs),2.0)*( 4.0*(Sps) + 8.0*(1.0+pow(cos(fonfs),2.0))*(Sacc/pow(2.0*PI*f,4.0)*(1.0+red)) ) / pow(2.0*Larm,2.0);
    
}

double Likelihood_check(struct Data *dat, struct Het *het, int ll, double *params)
{
    double *AS, *ES;
    double HH, HD;
    double logL;
    int imax, imin;
    
    imin = het->MN;
    imax = het->MM;
    
    AS = (double*)malloc(sizeof(double)* (dat->N));
    ES = (double*)malloc(sizeof(double)* (dat->N));
    
    ResponseFast(dat, ll, params, AS, ES);
    
    HH = (fourier_nwip2(AS, AS, dat->SN[0], imin, imax, dat->N)+fourier_nwip2(ES, ES, dat->SN[1], imin, imax, dat->N));
    HD = (fourier_nwip2(AS, dat->data[0], dat->SN[0], imin, imax, dat->N)+fourier_nwip2(ES, dat->data[1], dat->SN[1], imin, imax, dat->N));
    
    logL = HD-0.5*HH;
    
   //printf("%e %e %e\n", logL, HH, HD);
    
    free(AS);
    free(ES);
    
    return(logL);
    
}

double SNRstart(struct Data *dat, int ll, double *params)
{
    double *AS, *ES;
    double HH;
    double fstart;
    double SNRmin = 0.01;
    double HHmin;
    int i, j, imax;
    
    // This subroutine finds the frequency where we start to get some SNRb accumulation
    
    HHmin = SNRmin*SNRmin;
    
    
    imax = dat->N/2-1;
    
    AS = (double*)malloc(sizeof(double)* (dat->N));
    ES = (double*)malloc(sizeof(double)* (dat->N));
    
    ResponseFast(dat, ll, params, AS, ES);
    
    HH = 0.0;

    i = 0;
    do
    {
        i++;
        j = dat->N -i;
        HH += 4.0*( (AS[i]*AS[i]+AS[j]*AS[j])/dat->SN[0][i] + (ES[i]*ES[i]+ES[j]*ES[j])/dat->SN[1][i]);
        
    }while(i < imax && HH < HHmin);
    
    fstart = (double)(i)/dat->Tobs;
    
    free(AS);
    free(ES);
    
    return(fstart);
    
}

double Likelihood(struct Data *dat, int ll, double *params)
{
    double *AS, *ES;
    double HH, HD;
    double logL;
    
    AS = (double*)malloc(sizeof(double)* (dat->N));
    ES = (double*)malloc(sizeof(double)* (dat->N));
    
    ResponseFast(dat, ll, params, AS, ES);
    
    HH = (fourier_nwip(AS, AS, dat->SN[0], dat->N)+fourier_nwip(ES, ES, dat->SN[1], dat->N));
    HD = (fourier_nwip(AS, dat->data[0], dat->SN[0], dat->N)+fourier_nwip(ES, dat->data[1], dat->SN[1], dat->N));
    
    logL = HD-0.5*HH;
    
   //printf("%e %e %e\n", logL, HH, HD);
    
    free(AS);
    free(ES);
    
    return(logL);
    
}

double Likelihood_Slow(struct Data *dat, int ll, double *params)
{
    int i;
    double *AS, *ES;
    double HH, HD;
    double AA, EE;
    double logL;
    FILE *out;
    
    AS = (double*)malloc(sizeof(double)* (dat->N));
    ES = (double*)malloc(sizeof(double)* (dat->N));
    
    ResponseFreq(dat, ll, params, AS, ES);
    
    HH = (fourier_nwip(AS, AS, dat->SN[0], dat->N)+fourier_nwip(ES, ES, dat->SN[1], dat->N));
    HD = (fourier_nwip(AS, dat->data[0], dat->SN[0], dat->N)+fourier_nwip(ES, dat->data[1], dat->SN[1], dat->N));
    
    logL = HD-0.5*HH;
    
    free(AS);
    free(ES);
    
    return(logL);
    
}

void ResponseFreq(struct Data *dat, int ll, double *params, double *AS, double *ES)
{
    
    /*   Indicies   */
    int i,j, k, n, m, a, M, nn, nmin, nmax;
    
    /*   GW Source data   */
    double Mc, theta, phi, psi, D, iota, A, Aplus, Across, f0, fdot, phio;
    double costh, sinth, cosph, sinph, cosi, cosps, sinps;
    
    /*   Time and distance variables   */
    double xi, t;
    
    /*   Miscellaneous  */
    double xm, fstep, power, om, mx;
    
    double Amp, Phase, fonfs, f, x;
    
    double Aprime, Pprime, fi, fend;
    
    double HC, HS, hp, hc, Tobs;
    
    double m1, m2, chi1, chi2, phic, tc, distance, Mtot, eta, dm, fr, af;
    
    double *ta, *xia, *FF;
    
    double Fp, Fc, kdotx, delt, fmin, fmax;
    
    double XR, XI, YR, YI, ZR, ZI;
    
    double fstart, fstop, Tcut, fref;
    
    double sqrtTobs;
    
    int nfmin, nfmax, nf;
    
    int NA, N;
    
    clock_t start, end;
    double cpu_time_used;
    
    double m1_SI, m2_SI, deltaF;
    
    double *FpAR, *FpAI, *FcAR, *FcAI;
    double *FpER, *FpEI, *FcER, *FcEI;
    
    FILE *out;
    
    N = dat->N;
    Tobs = dat->Tobs;
    sqrtTobs = dat->sqrtTobs;

    
    // Have to generate full signal to get the merger phase correct
    // There are probably ways to work around this
    StartStop(ll, params, dat->Tstart, dat->Tend, dat->dt, &fstart, &fstop, &fr);
    
    // Because of the way the LDC phase is set at merger for massive BHs, we have to take the signal all the way out to
    // the merger frequency even if the obsevation time doesn't get us to merger. The signal is still
    // truncated by the window at Tend, so the SNRs etc will be correct
    if(sflag == 0) fstop = 2.0*fr;
    
   // printf("%e %e %e\n", fstart, fstop, fr);

    
    
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

    chi1 = params[2];  // Spin1
    chi2 = params[3];  // Spin2
    phic = params[4];  // merger phase
    tc = params[5];    // merger time
    cosi = params[10];  // inclination
    
    Aplus = 0.5*(1.+cosi*cosi);
    Across = -cosi;
    
    AmpPhaseFDWaveform *ap = NULL;
    double fRef_in;
    double *AF, *TF;
    int ret, flag1, flag2;
    
    AF = (double*)malloc(sizeof(double)* (N/2));
    TF = (double*)malloc(sizeof(double)* (N/2));

    
    nfmin = (int)(fstart*Tobs);
    if(nfmin < 1) nfmin = 1;
    nfmax = (int)(fstop*Tobs);
    if(nfmax <= nfmin) nfmax = nfmin+1;
    if(nfmax > N/2) nfmax = N/2;
    nf = nfmax-nfmin;
    
    fmin = (double)(nfmin)/Tobs;
    fmax = (double)(nfmax)/Tobs;
    
    deltaF = 1.0/Tobs;
    
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
    
    fRef_in = PDfref;
    
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
    
    
    // compute the Frequency time series
    timearray(params, freq, nf, TF, ap);
    
    RAantenna(params, nf, TF, FF, xia, FpAR, FpAI, FcAR, FcAI, FpER, FpEI, FcER, FcEI);
    
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
    
    // printf("%d %d\n", nfmin, nfmax);
    
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
            
            // Tukey filter to match what is done in time domain
            x = 1.0;
            if(t > dat->Tend || t < dat->Tstart)
            {
                x = 0.0;
            }
            else
            {
             if(t < dat->Tstart+t_tuke) x = 0.5*(1.0+cos(PI*((t-dat->Tstart)/t_tuke-1.0)));
             if(t > dat->Tend-t_tuke && t < dat->Tend) x = 0.5*(1.0-cos(PI*(t-dat->Tend)/t_tuke));
            }
             
            
            kdotx = t-xi;
            
            Amp = x*AF[m];
            Phase = ap->phase[m]+2.0*phic;
            
            HC = Amp*cos(2.0*PI*f*(Tobs-tc+dat->dt/2.0-kdotx)-Phase);
            HS = Amp*sin(2.0*PI*f*(Tobs-tc+dat->dt/2.0-kdotx)-Phase);
            

            AS[n] = FpAR[m]*Aplus*HC - FpAI[m]*Aplus*HS - FcAR[m]*Across*HS - FcAI[m]*Across*HC;
            AS[N-n] = FpAI[m]*Aplus*HC + FpAR[m]*Aplus*HS - FcAI[m]*Across*HS + FcAR[m]*Across*HC;
            
            ES[n] = FpER[m]*Aplus*HC - FpEI[m]*Aplus*HS - FcER[m]*Across*HS - FcEI[m]*Across*HC;
            ES[N-n] = FpEI[m]*Aplus*HC + FpER[m]*Aplus*HS - FcEI[m]*Across*HS + FcER[m]*Across*HC;
           
            
        }
        
        
    }
    
    
    /*   Deallocate Arrays   */
        
        DestroyAmpPhaseFDWaveform(ap);
        DestroyRealVector(freq);

    
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
    
    return;
}


void ResponseFast(struct Data *dat, int ll, double *params, double *AS, double *ES)
{
    double *AAmp, *EAmp, *APhase, *EPhase;
    double *APC, *APS, *EPC, *EPS;
    double af, fr, df, DT, fac, deltaF, f, fmax, fmin, x;
    double m1, m2, chi1, chi2, Mtot, Mc, eta;
    double Amp, Phase, tf, fguess, kxm;
    double m1_SI, m2_SI, distance, tc, phic;
    double cp, sp, cpx, spx;
    int i, NF;
    double A, P;
    double px, fnew;
    
    int NFmax = 100000;
    
    double *AF, *PF, *FF, *TF;
    
    FILE *out;
    
    FF = (double*)malloc(sizeof(double)* (NFmax));
    
    SetUp(dat, ll, params, NFmax, &NF, FF);
    
    AF = (double*)malloc(sizeof(double)* (NF));
    PF = (double*)malloc(sizeof(double)* (NF));
    TF = (double*)malloc(sizeof(double)* (NF));
    
    Intrinsic(ll, params, dat->Tobs, NF, FF, TF, PF, AF);
    
    AAmp = (double*)malloc(sizeof(double)* (NF));
    EAmp = (double*)malloc(sizeof(double)* (NF));
    APhase = (double*)malloc(sizeof(double)* (NF));
    EPhase = (double*)malloc(sizeof(double)* (NF));
    
    Extrinsic(params, dat->Tstart, dat->Tend, NF, FF, TF, PF, AF, AAmp, EAmp, APhase, EPhase, &kxm);
    
    APC = (double*)malloc(sizeof(double)* (NF));
    EPC = (double*)malloc(sizeof(double)* (NF));
    APS = (double*)malloc(sizeof(double)* (NF));
    EPS = (double*)malloc(sizeof(double)* (NF));
    
    
    /*
     out = fopen("inAP.dat","w");
     for (i=0; i< NF; i++)
     {
     fprintf(out,"%e %e %e %e\n", FF[i], TF[i], PF[i], AF[i]);
     }
     fclose(out);
     */
    
     /*
     out = fopen("fastAPr.dat","w");
     for (i=0; i< NF; i++)
     {
     fprintf(out,"%e %e %e %e %e\n", FF[i], AAmp[i], APhase[i], EAmp[i], EPhase[i]);
     }
     fclose(out);
      */
    
    // safer to interpolate the sine and cosine of the phases
    // because of difficulty catching all the phase wraps
    // Removing intrinsic phase evolution, to be restored later
        for (i=0; i< NF; i++)
        {
         APC[i] = cos(APhase[i]+PF[i]);
         APS[i] = sin(APhase[i]+PF[i]);
         EPC[i] = cos(EPhase[i]+PF[i]);
         EPS[i] = sin(EPhase[i]+PF[i]);
        }
    
    gsl_interp_accel *Pacc = gsl_interp_accel_alloc();
    gsl_spline *Pspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(Pspline, FF, PF, NF);
    
    gsl_interp_accel *PACacc = gsl_interp_accel_alloc();
    gsl_spline *PACspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(PACspline, FF, APC, NF);
    
    gsl_interp_accel *PECacc = gsl_interp_accel_alloc();
    gsl_spline *PECspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(PECspline, FF, EPC, NF);
    
    gsl_interp_accel *PASacc = gsl_interp_accel_alloc();
    gsl_spline *PASspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(PASspline, FF, APS, NF);
    
    gsl_interp_accel *PESacc = gsl_interp_accel_alloc();
    gsl_spline *PESspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(PESspline, FF, EPS, NF);
    
    gsl_interp_accel *AAacc = gsl_interp_accel_alloc();
    gsl_spline *AAspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(AAspline, FF, AAmp, NF);
    
    gsl_interp_accel *AEacc = gsl_interp_accel_alloc();
    gsl_spline *AEspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(AEspline, FF, EAmp, NF);
    
    
    
    phic = params[4];  // merger phase
    tc = params[5];    // merger time
    
    AS[0] = 0.0;
    ES[0] = 0.0;
    AS[dat->N/2] = 0.0;
    ES[dat->N/2] = 0.0;
    
    for (i=1; i< dat->N/2; i++)
    {
        f = (double)(i)/dat->Tobs;
        AS[i] = 0.0;
        ES[i] = 0.0;
        AS[dat->N-i] = 0.0;
        ES[dat->N-i] = 0.0;
        
        if(f > FF[0] && f < FF[NF-1])
        {
            // put in time/phase shift and restore intrinsic phase evolution
            px = 2.0*PI*f*(dat->Tobs-tc + dat->dt/2.0)-2.0*phic-gsl_spline_eval (Pspline, f, Pacc);
            cpx = cos(px);
            spx = sin(px);
            
            cp = gsl_spline_eval (PACspline, f, PACacc);
            sp = gsl_spline_eval (PASspline, f, PASacc);
            A = gsl_spline_eval (AAspline, f, AAacc);
            AS[i] = A*(cp*cpx-sp*spx);
            AS[dat->N-i] = A*(cp*spx+sp*cpx);
            
            cp = gsl_spline_eval (PECspline, f, PACacc);
            sp = gsl_spline_eval (PESspline, f, PASacc);
            A = gsl_spline_eval (AEspline, f, AAacc);
            ES[i] = A*(cp*cpx-sp*spx);
            ES[dat->N-i] = A*(cp*spx+sp*cpx);
        }
        
    }
    
    
    free(TF);
    free(FF);
    free(PF);
    free(AF);

    free(AAmp);
    free(EAmp);
    free(APhase);
    free(EPhase);
    free(APC);
    free(EPC);
    free(APS);
    free(EPS);
    
    gsl_spline_free(Pspline);
    gsl_spline_free(PACspline);
    gsl_spline_free(PECspline);
    gsl_spline_free(PASspline);
    gsl_spline_free(PESspline);
    gsl_interp_accel_free(Pacc);
    gsl_interp_accel_free(PACacc);
    gsl_interp_accel_free(PECacc);
    gsl_interp_accel_free(PASacc);
    gsl_interp_accel_free(PESacc);
    
    gsl_spline_free(AAspline);
    gsl_spline_free(AEspline);
    gsl_interp_accel_free(AAacc);
    gsl_interp_accel_free(AEacc);
    
    
}

void timearray(double *params, RealVector *freq, long N, double *TF, AmpPhaseFDWaveform *ap)
{
    
    int flag;
    int i, j;
    double v;
    double tc,  deltaF, fmax;
    
    tc = params[5];    // merger time

    
    for (i = 1; i < N-1; ++i) TF[i] = ( ap->phase[i+1]- ap->phase[i-1])/(2.0*PI*(freq->data[i+1]-freq->data[i-1])) + tc;
    
    TF[0] = TF[1];
    TF[N-1] = TF[N-2];
    
    j = N-1;
    flag = 0;
    for (i = 0; i < N-1; ++i)
    {
        // catch where time turns over
        if(TF[i+1] < TF[i] && flag == 0)
        {
            j = i;
            flag = 1;
        }
        // don't allow time to go too far into the past
        if(TF[i] < -Tpad) TF[i] = -Tpad;
    }
    
    // freeze time at turn over
    for (i = j; i < N; ++i)
    {
        TF[i] = TF[j];
    }
    
}



void StartStop(int ll, double *params, double Tstart, double Tend, double dt, double *fstart, double *fstop, double *frg)
{
    
    double m1, m2, m1_SI, m2_SI, chi1, chi2, tc, dm;
    double Mtot, eta, Mc, af, fr;
    double Amp, Phase, distance;
    double fmin, fmax;
    double fnew, tf, fny, Tseg;
    int i;
    
    Tseg = Tend-Tstart;
    
    // Here we space the frequency array to give approximately equal spacing in time
    // The dynamic frequency spacing is capped to be between 1e-6 and 1e-4 Hz

    fmin = FofT(ll, Tseg, params, &fr, dt, Tstart);

    // nan catcher
    if(fmin != fmin) fmin = 1.0/Tseg;
    if(fmin < 0.0) fmin = 1.0/Tseg;
    
    fmax = FofT(ll, Tseg, params, &fr, dt, Tend);
    
    // nan catcher
    if(fmax != fmax) fmax = 2.0*fr;
    if(fmax < fmin) fmax = 2.0*fmin;
    
    //printf("%e %e %e\n", fmin, fmax, fr);
    
    *frg = fr;
    *fstart = fmin;
    *fstop = fmax;
    
}

void Intrinsic(int ll, double *params, double Tobs, int NF, double *FF, double *TF, double *PF, double *AF)
{
    
    AmpPhaseFDWaveform *ap = NULL;
    RealVector *freq;
    
    double m1, m2, chi1, chi2, Mtot, Mc;
    double m1_SI, m2_SI, distance, tc, phic;
    double eta, dm, fx;
    
    double af, fr, dtdf, sqrtTobs;
    
    int i, ii, ret, flag;
    double fonfs, t, told, tx;
    
    sqrtTobs = sqrt(Tobs);

    double fRef_in;
    
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
    Mc = exp(params[0])*TSUN;
    Mtot = exp(params[1])*TSUN;
    distance = exp(params[6])*1.0e9*PC_SI; // distance
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
    
    fRef_in = PDfref;

    
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
    
    
    // need to capture turn-over in the SPA time (if it happens)
    flag = 0;
    told = ap->time[0]+tc;
    TF[0] = told;
    i = 1;
    do
      {
          t = ap->time[i]+tc;
          
          TF[i] = t;
          
          if(t < told && flag == 0)
          {
              flag = 1;
              ii = i-1;
              tx = told;
          }
          
          told = t;
          i++;
          
      } while( i< NF && flag == 0);
    
      // there was a turn-over
    if(flag == 1)
      {
          int NR=10;
          double df;
          RealVector *fref;
          AmpPhaseFDWaveform *ar = NULL;
          df = (FF[ii+1] - FF[ii])/(double)(NR-1);
          fref = CreateRealVector(NR);
          for (i=0; i< NR; i++) fref->data[i] = FF[ii]+(double)(i)*df;
          ret = IMRPhenomDGenerateh22FDAmpPhase(&ar,fref,0.0,fRef_in,m1_SI,m2_SI,chi1,chi2,distance);
          
          flag = 0;
          told = ar->time[0]+tc;
          i = 1;
          do
            {
                
                t = ar->time[i]+tc;
                
                if(t < told && flag == 0)
                {
                    flag = 1;
                    fx = fref->data[i-1];
                    tx = told;
                }
                
                told = t;
                i++;
                
            } while( i< NR && flag == 0);
          
          DestroyAmpPhaseFDWaveform(ar);
          DestroyRealVector(fref);
          
          for (i=ii; i< NF; i++)
            {
                TF[i] = tx + (FF[i]-fx)*dtdf;
            }
      }
   
    
    for (i=0; i< NF; i++)
    {
        PF[i] = ap->phase[i];
        AF[i] =  h22fac*ap->amp[i]/sqrtTobs;
        fonfs = freq->data[i]/fstar;
        AF[i] *= (4.0*fonfs*sin(fonfs));   // conversion to fractional frequency and leading order TDI transfer function
        
        //printf("%d %.12e %.12e\n", i, FF[i], TF[i]);
        
    }
    
    DestroyAmpPhaseFDWaveform(ap);
    DestroyRealVector(freq);

    
}



/*************************************************************************/
/*        Rigid approximation position of each LISA spacecraft           */
/*************************************************************************/
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

void lisaskyloc(double t, double *params, double *thetaL, double *phiL)
{
    double *x, *y, *z;
    double *r12, *r13, *yL;
    double *ov, *nv;
    double *sv;
    double nm;
    double costh, sinth, phi;
    double sinph, cosph, p;
    double ndo, xdo, ydo;
    int i;
    
    ov = double_vector(4); nv = double_vector(4); sv = double_vector(4);
    x = double_vector(4); y = double_vector(4); z = double_vector(4);
    r12 = double_vector(4); r13 = double_vector(4); yL = double_vector(4);
    
    costh = params[7];   sinth = sqrt(1.0-costh*costh);
    phi = params[8];   // EclipticLongitude
    sinph = sin(phi);
    cosph = cos(phi);
    
    ov[1] = sinth*cosph;
    ov[2] = sinth*sinph;
    ov[3] = costh;
    
    spacecraft(t, x, y, z);
    
    r12[1] = (x[2] - x[1])/Larm;   r13[1] = (x[3] - x[1])/Larm;
    r12[2] = (y[2] - y[1])/Larm;   r13[2] = (y[3] - y[1])/Larm;
    r12[3] = (z[2] - z[1])/Larm;   r13[3] = (z[3] - z[1])/Larm;
    
    // vector orthogonal to LISA plane
    nv[1] = r12[2]*r13[3]-r12[3]*r13[2];
    nv[2] = r12[3]*r13[1]-r12[1]*r13[3];
    nv[3] = r12[1]*r13[2]-r12[2]*r13[1];
    
    // magnitude
    nm = sqrt(nv[1]*nv[1]+nv[2]*nv[2]+nv[3]*nv[3]);
    
    // turn into unit vector
    for(i=1;i<=3;i++) nv[i] /= nm;
    
    // LISA frame is defined with nv as the z-axis and r12 as the x-axis.
    
    // y axis of LISA frame
    yL[1] = nv[2]*r12[3]-nv[3]*r12[2];
    yL[2] = nv[3]*r12[1]-nv[1]*r12[3];
    yL[3] = nv[1]*r12[2]-nv[2]*r12[1];
    
    ndo = 0.0;
    for(i=1;i<=3;i++) ndo += nv[i]*ov[i];
    
    *thetaL = asin(ndo);
    
    //project out component of ov in z direction
    for(i=1;i<=3;i++) ov[i] -= ndo*nv[i];
    
    nm = sqrt(ov[1]*ov[1]+ov[2]*ov[2]+ov[3]*ov[3]);
    // turn projected location vector into unit vector
    for(i=1;i<=3;i++) ov[i] /= nm;
    
    // angle of projected sky location wrt to x axis (chosen to be r12) in LISA frame
    xdo = 0.0;
    for(i=1;i<=3;i++) xdo += r12[i]*ov[i];
    
    // angle of projected sky location wrt to x axis (chosen to be r12) in LISA frame
    ydo = 0.0;
    for(i=1;i<=3;i++) ydo += yL[i]*ov[i];
    
    p = atan2(ydo,xdo);
    
    if(p < 0.0) p += TPI;
    
    *phiL = p;

    free_double_vector(ov); free_double_vector(nv); free_double_vector(sv);
    free_double_vector(x); free_double_vector(y); free_double_vector(z);
    free_double_vector(r12); free_double_vector(r13); free_double_vector(yL);
    
    
}


// merger time at guiding center of detector
double Tmerger(double *params, double t)
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
    
    spacecraft(t, x, y, z);
        
        // guiding center
        xa = (x[1]+x[2]+x[3])/3.0;
        ya = (y[1]+y[2]+y[3])/3.0;
        za = (z[1]+z[2]+z[3])/3.0;
        
    kdotx = (xa*kv[1]+ya*kv[2]+za*kv[3])/clight;
    
    td = t+kdotx;
    
    free_double_vector(kv);
    free_double_vector(x); free_double_vector(y); free_double_vector(z);
    
    return(td);
}

double f_at_t(double m1, double m2, double chi1, double chi2, double tc, double dt, double t)
{
    // 3PN f(t)
    int i;
    double f, fr, fny, af, M, eta, chi, theta;
    double gamma_E=0.5772156649; //Euler's Constant-- shows up in 3PN term
    double PN1, PN15, PN2, PN25, PN3, PN35;
    double theta2, theta3, theta4, theta5, theta6, theta7;
    
    fny = 1.0/(2.0*dt);
    
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
    
    if(f > fny) f = fny;
    
    return(f);
    
}



double FofT(int ll, double Tobs, double *params, double *frg, double dt, double tref)
{
    
    double m1, m2, m1_SI, m2_SI, chi1, chi2, tc, dm;
    double Mtot, eta, Mc, af, fr;
    double Amp, Phase, distance;
    double fref;
    double fnew, tf;
    int i;
    
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
    
    chi1 = params[2];
    chi2 = params[3];
    tc = params[5];

    af = FinalSpin0815(eta, chi1, chi2);
    fr = fring(eta, chi1, chi2, af)/Mtot;
    
    *frg = fr;
    
    //printf("%e %e %e\n", m1, m2, Mtot);
    
    
    // Here we space the frequency array to give approximately equal spacing in time
    // The dynamic frequency spacing is capped to be between 1e-6 and 1e-4 Hz
    
    if(tref > tc)
    {
        fref = 2.0*fr;
    }
    else
    {
    
    // guess at f of t
        
    fref =  f_at_t(m1, m2, chi1, chi2, tc, dt, tref);
    
    // find the frequency at t= tref.
    i = 0;
    do
    {
        getfreq(Tobs, &fnew, &tf, &Amp, &Phase, tref, fref, 0.0, PDfref, m1_SI, m2_SI, chi1, chi2, distance, tc);
        //printf("%e %e %e %e\n", fref, fnew, tf, tref);
        fref = fnew;
        i++;
    }while(i < 10 && fabs(tf-tref) > 1.0 && fref == fref);
        
    }
    
    return(fref);
    
}

void getfreq(double Tobs, double *fnew, double *tf, double *Amp, double *Phase, double t, double fguess, double phic, double fRef_in, double m1_SI, double m2_SI, double chi1, double chi2, double distance, double tc)
{
    AmpPhaseFDWaveform *ap = NULL;
    double ep, u, v, tnew, x;
    double delT, delF, dtdf, fonfs, sqrtTobs;
    int ret;
    double M_sec;
    RealVector *freq;
    
    sqrtTobs = sqrt(Tobs);
    
    M_sec = (m1_SI+m2_SI) * MTSUN_SI/MSUN_SI;
    
    if(fguess < 1.0/Tobs) fguess = 1.0/Tobs;
    
    ep = 1.0e-6/M_sec;
    v = (4.0*PI*ep);
    u = (2.0*PI*ep*ep);
    
    if(fguess-ep < 0.0) ep = 0.5*fguess;
    
    
    freq = CreateRealVector((3));
    
    freq->data[0] = fguess-ep;
    freq->data[1] = fguess;
    freq->data[2] = fguess+ep;
    
    ret = IMRPhenomDGenerateh22FDAmpPhase(&ap,freq,phic,fRef_in,m1_SI, m2_SI, chi1, chi2, distance);
    
    tnew = (ap->phase[2]-ap->phase[0])/v +tc;
    
    dtdf = (ap->phase[2]+ap->phase[0]-2.0*ap->phase[1])/u;
    
    delT = t-tnew;
    
    delF = delT/dtdf;
    
    *fnew = fguess + delF;
    
    *tf = tnew;
    
    x = h22fac*ap->amp[1]/sqrtTobs;
    fonfs = fguess/fstar;
    x *= 4.0*fonfs*sin(fonfs);   // conversion to fractional frequency and leading order TDI transfer function

    *Amp = x;
    
    *Phase = ap->phase[1];
    
    DestroyAmpPhaseFDWaveform(ap);
    DestroyRealVector(freq);
    
}




void map_params(int ll, double *params)
{
    
    double Mtot, Mc;
    double m1, m2, chi1, chi2;
 
    if(ll ==1) // switch to log for masses and distance
    {
    params[0] = log(params[0]);
    params[1] = log(params[1]);
    params[6] = log(params[6]);
    }
    
    if(ll ==2) // // switch to log for Mc, Mtot and distance
      {
       Mtot = params[0]+params[1];
       Mc = pow(params[0]*params[1], 3.0/5.0)/pow(Mtot,1.0/5.0);
       params[0] = log(Mc);
       params[1] = log(Mtot);
       params[6] = log(params[6]);
          
        //printf("%e %e\n", Mc, Mtot);
      }
    
    while(params[4] > PI)   params[4] -= PI;
    while(params[4] < 0.0)  params[4] += PI;
    while(params[8] > TPI)  params[8] -= TPI;
    while(params[8] < 0.0)  params[8] += TPI;
    while(params[9] > PI)   params[9] -= PI;
    while(params[9] < 0.0)  params[9] += PI;
    
    if(ll == 0 || ll == 1)
    {
    if(params[1] > params[0])  // catch if m2 > m1 and flip
    {
        m1 = params[1];
        chi1 = params[3];
        m2 = params[0];
        chi2 = params[2];
        params[0] = m1;
        params[1] = m2;
        params[2] = chi1;
        params[3] = chi2;
    }
    }
    
}

void RAantenna(double *params, int NF, double *TF, double *FF, double *xi, double *FpAR, double *FpAI, double *FcAR, double *FcAI,
               double *FpER, double *FpEI, double *FcER, double *FcEI)
{
    
    /*   Indicies    */
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

        FpAR[n] = (2.0*fpx-fpy-fpz)/3.0;
        FcAR[n] = (2.0*fcx-fcy-fcz)/3.0;
        
        FpER[n] = (fpz-fpy)/sq3;
        FcER[n] = (fcz-fcy)/sq3;
                   
        fpx = -0.5*( (dplus[1][2]*cosps+dcross[1][2]*sinps)*TI[1][2] - (dplus[1][3]*cosps+dcross[1][3]*sinps)*TI[1][3] );
        fcx = -0.5*( (-dplus[1][2]*sinps+dcross[1][2]*cosps)*TI[1][2] - (-dplus[1][3]*sinps + dcross[1][3]*cosps)*TI[1][3] );
                                         
        fpy = -0.5*( (dplus[2][3]*cosps+dcross[2][3]*sinps)*TI[2][3] - (dplus[2][1]*cosps+dcross[2][1]*sinps)*TI[2][1] );
        fcy = -0.5*( (-dplus[2][3]*sinps+dcross[2][3]*cosps)*TI[2][3] - (-dplus[2][1]*sinps + dcross[2][1]*cosps)*TI[2][1] );
                                                               
        fpz = -0.5*( (dplus[3][1]*cosps+dcross[3][1]*sinps)*TI[3][1] - (dplus[3][2]*cosps+dcross[3][2]*sinps)*TI[3][2] );
        fcz = -0.5*( (-dplus[3][1]*sinps+dcross[3][1]*cosps)*TI[3][1] - (-dplus[3][2]*sinps + dcross[3][2]*cosps)*TI[3][2] );
                                                                                     
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


