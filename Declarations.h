
struct MBH_Data
{
    int N;
    int Nch;
    double dt;
    double Tobs;
    double sqrtTobs;
    double Tstart;
    double Tend;
    double fmin;
    double fmax;
    double **SN;
    double **SM;
    double **data;
};

struct Het
{
    int M;
    int MM; //!<Maximum frequency bin of heterodyne likelihood
    int MN; //!<Minimum frequency bin of heterodyne likelihood
    int J;
    int NR;
    int Nch;
    int *fgrid;
    double SNR;
    double *freq;
    double *logLR;
    double *DD;
    double ***IP;
    double ***SL;
    double **amp;
    double **phase;
    double **rc;
    double **rs;
    double **dc;
    double **ds;
    double **aa;
    double ***lc;
    double ***ls;
    double ***ldc;
    double ***lds;
    double *pref;
};

void freehet(struct Het *het);
void lisaskyloc(double t, double *params, double *thetaL, double *phiL);
void SetUp(struct MBH_Data *dat, int ll, double *params, int *NFS, double *FF);
void StartStop(int ll, double *params, double Tstart, double Tend, double dt, double *fstart, double *fstop, double *frg);
void Intrinsic(int ll, double *params, double Tobs, int NF, double *FF, double *TF, double *PF, double *AF);
void ResponseFast(struct MBH_Data *dat, int ll, double *params, double *AS, double *ES);
void ResponseFreq(struct MBH_Data *dat, int ll, double *params, double *AS, double *ES);
double chisq(struct MBH_Data *dat, int ll, double *params, double *AR, double *ER);
double chisq_het(struct MBH_Data *dat, struct Het *het, int ll, double *params, double **ampR, double **phaseR);
void heterodyne(struct MBH_Data *dat, struct Het *het, int ll, double *params);
void legendre_maker(int J, int U, double **P);
void fullphaseamp(struct MBH_Data *dat, int ll, int K, double *params, double *freq, double *Aamp, double *Eamp, double *Aphase, double *Ephase);
double log_likelihood_het(struct MBH_Data *dat, struct Het *het, int ll, double *params, double *sx);
double Fstat_het(struct MBH_Data *dat, struct Het *het, int ll, double *params, double *sx, double tm);
double SNRstart(struct MBH_Data *dat, int ll, double *params);
void FisherHet(struct MBH_Data *dat, struct Het *het, int ll, double *params, double **Fisher);
void FisherSubHet(struct MBH_Data *dat, struct Het *het, int ll, int *pmap, double *params, double **Fisher);

void map_params(int ll, double *params);
int *int_vector(int N);
void free_int_vector(int *v);
double **double_matrix(int N, int M);
void free_double_matrix(double **m, int N);
int **int_matrix(int N, int M);
void free_int_matrix(int **m, int N);
double *double_vector(int N);
void free_double_vector(double *v);
void FisherFastPE(double *params);
double ***double_tensor(int N, int M, int L);
void free_double_tensor(double ***t, int N, int M);
void Inverse(double **M, double **IM, int d);
void spacecraft(double t,  double *x, double *y, double *z);
void RAantenna(double *params, int NF, double *TF, double *FF, double *xi, double *FpAR, double *FpAI, double *FcAR, double *FcAI,
               double *FpER, double *FpEI, double *FcER, double *FcEI);
void timearray(double *params, RealVector *freq, long N, double *TF, AmpPhaseFDWaveform *ap);
double FofT(int ll, double Tobs, double *params, double *frg, double dt, double tref);
void Extrinsic(double *params, double Tstart, double Tend, int NF, double *FF, double *TF, double *PF, double *AF, double *AAmp, double *EAmp, double *APhase, double *EPhase, double *kxm);
void efix(struct MBH_Data *dat, struct Het *het, int hr, int ll, double *params, double *min, double *max, double *eval, double **evec, double zs);
void het_space(struct MBH_Data *dat, struct Het *het, int ll, double *params, double *min, double *max);
void instrument_noise(double f, double *SAE);
void getfreq(double Tend, double *fnew, double *tf, double *Amp, double *Phase, double t, double fguess, double phic, double fRef_in, double m1_SI, double m2_SI, double chi1, double chi2, double distance, double tc);
void update(struct MBH_Data *dat, struct Het *het, int typ, int k, int ll, double *logLx, double **paramx, double **paramy, double **sx, double **sy, double *min, double *max, int *who, double *heat, double ***history, int NH, double **ejump, double ***evec, int **cv, int **av, gsl_rng *r);
void FisherEvec(double **fish, double *ej, double **ev, int d);
double fourier_nwip2(double *a, double *b, double *Sn, int imin, int imax, int N);
void FisherSub(struct MBH_Data *dat, int ll, int *pmap, double *params, double **Fisher);
void FisherFast(struct MBH_Data *dat, int ll, double *params, double **Fisher);
void de_jump(double *paramsx, double *paramsy, double **history, int m, int d, gsl_rng *r);
double det(double **A, int N);
double Tmerger(double *params, double t);
void get_component_masses(double *params, int flag, double *m1, double *m2);
void print_mbh_chain_file(struct MBH_Data *dat, struct Het *het, int *who, double **paramx, double *logLx, double **sx, int ll, int mc, FILE *chain);
void set_mbh_priors(struct MBH_Data *dat, int massFlag, double *min, double *max);

