seed = -1 *seed for random numbers -1 for random

seqfile = ../../partitioned_alignment.phy  *alignment file
treefile = ../../calibrated_tree.nwk *topology tree file
mcmcfile = mcmc.mcmc *mcmc output
outfile = mcmc.out *stdout output

ndata = 4 *number of partitions in alignment file
seqtype = 2 *0-2 nucleotides, codons, amino acids
usedata = 2 in.BV *0-3 no data, seq, approximate in.bv
clock = 3 *1-3 strict clock, uncorrelated rates, autocorrelated rates

model = 0 *0-10 JC69, K80, F81, F84, HKY85, T92, TN93, REV, UNREST, REVu, UNRESTu
alpha = 0.5 *alpha value
ncatG = 5 *number discrete gamma categories
cleandata = 0 *0,1 leave ambiguous sites, remove them

BDparas = 1 1 0.1 * Birth, Death, Sampling
kappa_gamma = 1 1 * alpha, beta of kappa prior gamma - irrelevant here
alpha_gamma = 1 1 * alpha, beta of alpha prior gamma
rgene_gamma = 2 39.513197 1 * total tree length / root age (mean) scaled so that alpha = 2
sigma2_gamma = 1 1 * alpha, beta branch rate variation prior

finetune = 1: .1 .1 .1 .1 .1 .1 *magical trickery
print = 1 *0-2 no mcmc sample, everything except branch lengths, everything
burnin = 200000 *n samples discarded
sampfreq = 20000 *how often to sample chain
nsample = 20000 *number of samples before finish
