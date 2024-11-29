**Divergence time estimation using MCMCTree with multiple loci.**

1. Preprocess your gene families

You will need alignments of all the gene (protein) families you are using (trimmed if desired), a best fitting substitution model, and a gene tree with branch lengths for each.

I would use IQTree because of its fancy model selection and the subsequent bits have been set up for IQTree

1. Partition your genes into families based on rate.

First, run extract_tree_lengths.py – you will need python modules glob and ete3 installed. You will also need to modify this line – this isn’t ideal and I’ll probably update it soon but I’m a busy boy.

treefiles = glob.glob("Pinniped_phylogeny_files/gene_trees/trees/\*.treefile")

to point to where your treefiles are (as if listing them using bash.

Run using

python3 extract_tree_lengths.py

And it should make a file called tree_lengths.csv.

Included in this directory is an R script (cluster_fams.R) that I suggest you run line by line in Rstudio so that you can look modify it on the fly.

First you need to modify the working directory via the line

setwd("~/OneDrive - The University of Nottingham/seals/clock")

It then takes your gene trees and draws a histogram of their heights (the maximum branch length between two tips. Have a look and decide on the number of partitions you think fit the data. For example, in the professional figure below, you would decide 0, 1, 2 and 3 breaks from left to right respectively.

![](data:image/png;base64,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)![](data:image/png;base64,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)![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGEAAABgCAYAAAANWhwGAAAAAXNSR0IArs4c6QAAAIRlWElmTU0AKgAAAAgABQESAAMAAAABAAEAAAEaAAUAAAABAAAASgEbAAUAAAABAAAAUgEoAAMAAAABAAIAAIdpAAQAAAABAAAAWgAAAAAAAABIAAAAAQAAAEgAAAABAAOgAQADAAAAAQABAACgAgAEAAAAAQAAAGGgAwAEAAAAAQAAAGAAAAAA3a/5DgAAAAlwSFlzAAALEwAACxMBAJqcGAAAAqlJREFUeAHtndFO20AUBWnV/+ryZThfFvNldC/gaqWm8dmIMQ+eK62y4JOzMCMHpDzk6eljWn+49vU2rKXvnYMILP2cEf64LzHOAQRG6Lf27YCf4fRH3AI/fu96ekIwgJ9wv/UBASUEkOiIEmjCQb8SAkh0RAk04aBfCQEkOqIEmnDQr4QAEh1RAk046FdCAImOKIEmHPQrIYBER5RAEw76lRBAoiNKoAkH/UoIINERJdCEg34lBJDoiBJowkG/EgJIdEQJNOGgXwkBJDqiBJpw0K+EABIdUQJNOOhXQgCJjiiBJhz0KyGAREeUQBMO+pUQQKIjSqAJB/1KCCDRESXQhIN+JQSQ6IgSaMJBvxICSHRECTThoF8JASQ6ogSacNCvhAASHVECTTjoV0IAiY4ogSYc9CshgERHlEATDvqVEECiI0qgCQf9Sggg0REl0ISDfiUEkOiIEmjCQb8SAkh0RAk04aBfCQEkOqIEmnDQr4QAEh1RAk046FdCAImOKIEmHPQrIYBER5RAEw76lRBAoiNKoAkH/UoIINGRH/2A+sCKvamcAxHwToDAztQqYYYWlFUCBHamVgkztKCsEiCwM7VKmKEFZZUAgZ2pVcIMLSirBAjsTK0SZmhBWSVAYGdqlTBDC8oqAQI7U6uEGVpQtiSsQXcLMkYeJOCd8CC4r3zar1722lfbKT3jp5OvO0yK2zbLtnnksd4xa32dEfIjvPaecxkCa9/X2h0l7CL60sD62VZ30fK5//tQ7zO7jmfwLqLuhJp6OWq1cQ4n8Ox/R4cz/+fAl01CvUY530Rgezmq4+tvgnM8gct2J9TRz32ttXGOJTDeCdvJrW9qbbP2Ta0zTdv5ZcfrLzvZe5cv/eJyS8K9J3nt/wSW4dLvvm/D17e27wLqwh/20JI6/JGlCwAAAABJRU5ErkJggg==)![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGkAAABkCAYAAACFHB7kAAAAAXNSR0IArs4c6QAAAIRlWElmTU0AKgAAAAgABQESAAMAAAABAAEAAAEaAAUAAAABAAAASgEbAAUAAAABAAAAUgEoAAMAAAABAAIAAIdpAAQAAAABAAAAWgAAAAAAAABIAAAAAQAAAEgAAAABAAOgAQADAAAAAQABAACgAgAEAAAAAQAAAGmgAwAEAAAAAQAAAGQAAAAAzj8QUAAAAAlwSFlzAAALEwAACxMBAJqcGAAAA0JJREFUeAHtnAFu2zAQBJ2i/yr9srgvC/OylpeYAOGkMM/e87rAEFAkUac9eSYMHLTO4fA52ti9je3Psp3GMeNJCJzGc6xy1uMQxzAT+HGlfxvXY2OYCawr57tjVpNZ0LWVZH482geBl7HF6rk2oo5hIsBKMoHPtEVShpapFkkm8Jm2SMrQMtUiyQQ+0xZJGVqmWiSZwGfaIilDy1SLJBP4TFskZWiZapFkAp9pi6QMLVMtkkzgM22RlKFlqkWSCXymLZIytEy1SDKBz7RFUoaWqRZJJvCZtkjK0DLVIskEPtMWSRlaplokmcBn2iIpQ8tUiyQT+ExbJGVomWqRZAKfaYukDC1TLZJM4DNtkZShZapFkgl8pi2SMrRMtUgygc+0RVKGlqkWSSbwmbZIytAy1SLJBD7TFkkZWqbakNQ3ereNGkqKCLCSisAqY5GkpFmUhaQisMrYkPS+Edg2aigpIsBKKgKrjA1JfSPw10YNJUUEWElFYJWxSFLSLMpCUhFYZez8Y4P8YUIlVXEWK0kMtCIOSRVUxZlIEgOtiENSBVVxJpLEQCvikFRBVZw5JfWN3LZRQ0kBgSmpIJpIFQEkqUgW5iCpEK4qGkkqkoU5SCqEq4qekvgndBXRgpwpqSCaSBUBJKlIFuYgqRCuKhpJKpKFOUgqhKuKRpKKZGEOkgrhqqKRpCJZmIOkQriqaCSpSBbmIKkQrioaSSqShTlIKoSrikaSimRhDpIK4aqip6S+EcgHyTYgVZRMSRXZZIoIzI++RBwffxFBVcewktREC/KQVABVHflTHXhDXrvhnp1b+k7R/1DzCEltAfG6HK/zy7T8sI/E+b+hTvL0fwe2i0vra7+49OU0nvc0Z5VvHNo5NPbz7fqcO196il2/8ynanfdnbj+O4p6V1C86tItzTrUE+og7ZiVpH4G0HQJH3t3tYDLXIMksYKc9P+52KHlrXtaV1L3PQvdvCPyOufX3pPdx3mLywaMX9GsFmY+M7KNZCIr9Yf1xF+dvY2txIB79nBffCPN47s+XSnZtpMYWI/PL5Ocd93/tS8T62pfpL4frPR8XLyXF5Glst7ygPu6Lcfkw/WP2Ob400WN0Uc5WzF87ykAvHfK1mgAAAABJRU5ErkJggg==)

Modify this line with the number of clusters you have identified (it’s set at 4 now)

breaks = getJenksBreaks(dists$dist, 4, subset = NULL)

It will then draw boundaries between four populations of tree heights, which I use as partitions in mcmctree. (note – currently the lines are hardcoded with some old values and I might have forgotten to change this, so you might need to change these yourself otherwise the figure won’t make sense (the program will still work as intended – it’s just for your understanding))

It will then write “clusters.csv”, which has the gene family name and the height and the cluster.

Finally, write your partitioned alignment file in phyllip format. Run

python3 generate_alignment.py &lt;alignments_directory&gt;

, replacing the &lt;alignments directory&gt; with your own directory containing all the alignments with the extention “.faa”. This makes a file “partitioned_alignment.phy” which is what you need.

1. Prepare your calibrated species tree

You need to use as many calibrations as possible – which is outside the scope of this manual (thank God). It’s down to you what to use and how to implement them, but I usually impose a uniform distribution with a hard minimum bound youngest possible age of the fossil and the a soft boundary at the maximum bound, which is usually the age at which we would expect to find evidence for the clade if it existed, but that we don’t.

I would modify the newick string for the species tree in a text editor. First remove the branch lengths and any clade labels. Next, you need a line that has the number of species at the start (not sure why – ask Ziheng), then for each calibrated node you can add a string in the following format.

'B(&lt;X&gt;,&lt;Y&gt;,1e-300,0.025)'

Where &lt;X&gt; is the minimum age and &lt;Y&gt; the soft maximum (1e-300 is the probability of the divergence time being lower than the minimum, so basically 0, and 0.025 is the probability of it being older than the maximum). I find it easiest to keep copying the newick string to FigTree and showing ranch labels to ensure I’ve been calibrating the correct nodes and that I haven’t accidentally messed with the toplogy). Your calibrated tree will look something like this.

21

(((((((bearded-seal,(hooded-seal,((largha-seal,harbour-seal),(ringed-seal,(grey-seal,caspian-seal)))))'B(1.32,2.729,1e-300,0.025)',((weddell-seal,(southern-elephant-seal,northern-elephant-seal))'B(0.505,2.729,1e-300,0.025)',hawaiian-monk-seal)'B(0.71,2.729,1e-300,0.025)')'B(1.382,2.729,1e-300,0.025)',(walrus,(northern-fur-seal,(stellar-sea-lion,californian-seal-lion))'B(0.5333,2.729,1e-300,0.025)')'B(1.597,2.729,1e-300,0.025)')'B(2.045,2.729,1e-300,0.025)',(sea-otter,sable)'B(1.957,4.88,1e-300,0.025)')'B(2.48,4.88,1e-300,0.025)',(polar-bear,giant-panda)'B(1.25,4.88,1e-300,0.025)')'B(3.38,4.88,1e-300,0.025)',dog)'B(3.771,6.609,1e-300,0.025)',cat)'B(3.771,6.609,1e-300,0.025)';

1. Run mcmctree for proteins.

Running mcmctree for protein data is complicated because of the approximate likelihood estimation it employs. You basically have to run it with the simplest possible model for each gene family first, then generate a hessian matrix for each, using the best fitting model then run it properly. For more information and an alternative guide, see [this handy github vignette by Sandra Álvarez Carretero](https://github.com/sabifo4/Tutorial_MCMCtree). Thankfully, I’ve got a bash script that runs all the different steps sequentially.

In your current directory, add the models directory and a directory within which to run the analysis.

In this latter directory set up 6 directories called clock_1 through clock_6.

You can also modify mcmctree.ctl – the control file - and add it to each of these 6 directories. I’ll go though each line you need to think about in the control file, whose contents follows (everything after “\*” is a comment).

seqfile = ../../partitioned_alignment.phy \*alignment file

treefile = ../../calibrated_tree.nwk \*topology tree file

mcmcfile = mcmc.mcmc \*mcmc output

outfile = mcmc.out \*stdout output

ndata = 4 \*number of partitions in alignment file

seqtype = 2 \*0-2 nucleotides, codons, amino acids

usedata = 2 in.BV \*0-3 no data, seq, approximate in.bv

clock = 3 \*1-3 strict clock, uncorrelated rates, autocorrelated rates

model = 0 \*0-10 JC69, K80, F81, F84, HKY85, T92, TN93, REV, UNREST, REVu, UNRESTu

alpha = 0.5 \*alpha value

ncatG = 5 \*number discrete gamma categories

cleandata = 0 \*0,1 leave ambiguous sites, remove them

BDparas = 1 1 0.1 \* Birth, Death, Sampling

kappa_gamma = 1 1 \* alpha, beta of kappa prior gamma - irrelevant here

alpha_gamma = 1 1 \* alpha, beta of alpha prior gamma

rgene_gamma = 2 39.513197 1 \* total tree length / root age (mean) scaled so that alpha = 2

sigma2_gamma = 1 1 \* alpha, beta branch rate variation prior

finetune = 1: .1 .1 .1 .1 .1 .1 \*magical trickery

print = 1 \*0-2 no mcmc sample, everything except branch lengths, everything

burnin = 200000 \*n samples discarded

sampfreq = 20000 \*how often to sample chain

nsample = 20000 \*number of samples before finish

So the first four lines are your input files and output files – you need to make sure they are pointing to the right place, recalling that the mcmctree.ctl file will be in the clock_x directories).

ndata needs to equal the number of partitions you chose. seqtype is 2 for amino acids. usedata is how the program treats the alignment, which will change automatically where appropriate.

Clock = 2 for uncorrelated rates. Clock = 3 for autocorrelated – this is a decision you have to make - a lot of people do both and share the results for both, hoping there isn’t much difference. [There is a complex way of selecting a model too](https://dosreislab.github.io/2018/08/31/bppr-stepstones.html).

Model needs to start as 0. ncatG is the number of gamma cateogies, which governs how sites can vary in rate within a partition, cleandata is about how you treat ambiguous data. These can all stay the same unless you have a reason to change them

BDparas is the birth and death parameters and kappa_gamma is about the transition-tranversion ratio. alpha_gamma is how rate can vary between the clusters. Don’t change these.

rgene_gamma is the expectation of rate and is important. The gamma refers to a distribution with 2 paramentes, alpha and beta, which is what you specify. Keep an alpha as 2, so the distribution is fairly defuse. Beta needs to be specified so that the mean of the distribution is the expected substation rate. You have to do some maths!

Mean_rate = tree_height / root_age

Tree_height I specify as the average tree height of a gene tree in the analysis.

Root_age is the mean of the prior you put on the root (from the fossil calibrations).

Beta = 2/mean_rate = what you should put in the second field of this line.

Leave sigma2_gamma, finetune and print

Burnin, samplefreq, and n_sample refer to how long to run the Markov chain for. The burnin is the number of samples before saving the values. Samplefreq is how often to take results. N_sample is how many samples to take. Unfortunately, there is no way of knowing what these should be until the run is finished, but when the analyses have finished, you can load all the .mcmc files in tracer. If the effective sample sizes are too low and/or the dates don’t all match, you will need to run the analysis again with a larger samplefreq. N_sample should always be 10,000 or 20,000 – you shouldn’t neeed any more to draw representative distributions.

When you’ve changed the file, put a copy in each clock_x directory.

You also need to make a file with the best model for each partition in the order the appear in your partitioned alignment file, which you can get using IQTree with modelfinderonly, on the alignment for each partition. The only likely models in the directory supplied are JTT, WAG, LG, Dayhoff, and Blosum62, so you can just which of these does best (regardless of any +G +I etc.) in the file. Call it best_models.txt.

Also included is mcmctree_proteins.sh, which is a slurm script that does everything you need to do. You will need to modify it to work with your job scheduler if it is not slurm. And also change this line:

for i in {0001..XXXX}

Change the XXXX to the number of gene families with leading “0”s. and

cd /your/favourite/directory/clock_${SLURM_ARRAY_TASK_ID}

with your directory. Put this in the same directory that your clock_1 through clock_6 are contained in and submit the job.

If you have any issues, just make sure the scripts are pointing to files that exist – you might have to add or remove a ../ or something. You can also ask me, but have a go at it first.

Done – turns out that’s a lot more long-winded than expected.
