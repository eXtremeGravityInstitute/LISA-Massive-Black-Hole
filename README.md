# Sangria

The first step is to read in the Sangria data and dump it to ASCII file. The code segmentSangria.c does this. It also produces a collection of shorter segments for the search code and outputs the data in both the time and frequency domain

gcc -o segmentSangria segmentSangria.c -lhdf5 -lgsl

./segmentSangria

The code also pads the data out to 16 times the length of one segment and FFTs the windowed and padded data. This can be used when running on the full data set. [TODO - change from radix 2 FFTs to general FFTs so less padding needed.] The segment code does not talk very long to run. Next up we need PSD estimates. This is done by invoking the script

source spec.sh 12

Here the "12" is the number of segments to work through. Each segment takes several minutes to run. Wavelet denoting is used in the initial PSD estimation. A short MCMC is used to improve the spline and line model. This code is taken from QuickCBC and adapted for LISA. The script also makes Qscans of the whitened A, E channel data for each segment.

The code SpecAverage.c is used to average and interpolate the PSD estimates from each segment to produce a PSD estimate for the full data set. Compile with

gcc -o SpecAverage SpecAverage.c -lgsl

And run with

./SpecAverage.

[TODO - understand where the errant factors of 2 are in the FFT and spectral estimate for the full, padded data set are coming from]

Next up is the search. This is invoked by calling the script

source search.sh 1 1

In the Sangria training data the first segment to include BH merger is number 1 (the numbering starts at 0). The arguments are the first and last segments searched. Here we are just searching segment 1 since the goal is to demonstrate how the code works, not repeat the entire analysis. To run the entire analysis type

source search.sh 0 11

The search code reads in the FFTed data for a segment and the PSD model derived by the SpecFit code. The code makes Scans along the way to demonstrate the signal removal. Note that the initial Qscan differs a little from the one produced by the SpecFit code since the former uses the median/line PSD with the latter uses the spline, Lorentizian fit.

The search code uses a PTMCMC and a maximized likelihood for each channel. It usually locks onto signals pretty quickly. For the first 1/8th of the intrinsic search, the code maximizes over time, phase and distance. For the rest of the intrinsic search it just maximizes over phase and distance, which is much faster (no FFTs). The next step is to find the sky location.

The search makes multiple passes over each segment (up to some pre-set maximum). It has a stopping criterion based on the maximum log likelihood reached. Once it stops finding signals or runs up against the maximum number of sweeps it stops and moves on the the next segment. The code returns the residual for the segment and the best fit values for the source parameters for each signal found. The sky search uses an F-statistic likelihood and the full response with no heterodyning, so it is relatively slow (even though it is only using about one month of data). On my laptop, each sweep of the intrinsic search takes about 15 minutes and each sky search takes about 15 minutes. The current settings are pretty lean - designed for speed rather than depth.

Once all the segments have been searched we extract the unique signals with SNR > 12 using the code unique.c

gcc -o unique unique.c -lgsl

./unique

This code produces the files search_sources.dat and source_info.dat. The first file contains the maximum likelihood solutions from the search. The second is a handy lookup table that tells you which source mergers in which segment.

Now we are ready to run a refined MCMC. This is compiled via

 clang -Xpreprocessor -fopenmp -lomp -w -o  PTMCMC PTMCMC.c Utils.c Response.c IMRPhenomD_internals.c IMRPhenomD.c -lgsl -lgslcblas  -lm

To run on segment 1 and source 0 you would type

./PTMCMC 1 0

To run on source 7 using the full data set you would type

./PTMCMC -1 7

You can either run on one of the 12 segments (each roughly 30 days) or on the full data set. The segments are label 0 through 11. Putting -1 for the segment number causes the full data set to be used. The search of the Sangria training data yields 15 unique sources numbered 0 through 14.

Using the full data set results in a longer set-up time (data read etc), but once the MCMC is off and running the evaluation time for the heterodyned likelihood is pretty much unchanged by the duration of the data set.

The PTMCMC outputs a single chain file, "chain.dat". [TODO - name the chain file by the source number and segment label]

