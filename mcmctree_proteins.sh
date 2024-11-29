#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32g
#SBATCH --time=7-00:00:00
#SBATCH --array 1-6


cd /your/favourite/directory/clock_${SLURM_ARRAY_TASK_ID}

rm in.BV
touch in.BV
sed -i "s/usedata = 2 in.BV/usedata = 3/g" mcmctree.ctl

j=1
/newhome/ab17362/paml4.9h/src/mcmctree
for i in {0001..XXXX}
do
    echo $i
    #Obtain the model for the partition
    family="$(head -n $j ../../order_of_genes_in_alignment.txt| tail -n 1)"
    model="$(grep -h "Best model according to BIC" ../../gene_alignments/${family}.fasta.phy.trimmed.prottest | cut -f 6 -d " " | cut -f 1 -d +)"
    j=$((j+1))
    #if statement to asign the model file
    if [ $model = "JTT" ]
    then
        modelfile="../../models/jones.dat"
    fi
    if [ $model = "WAG" ]
    then
        modelfile="../../models/wag.dat"
    fi
    if [ $model = "LG" ]
    then
        modelfile="../../models/lg.dat"
    fi
    if [ $model = "Blosum62" ]
    then
        modelfile="../../models/blosum62.dat"
    fi
    if [ $model = "Dayhoff" ]
    then
        modelfile="../../models/dayhoff.dat"
    fi



    #Modify tmp.ctl
    echo "tmp${i}.ctl before seds"
    cat tmp${i}.ctl
    sed -i "s/model = 0/model = 2/g" tmp${i}.ctl
    sed -i "s|aaRatefile =|aaRatefile = ${modelfile}|g" tmp${i}.ctl
    echo "fix_alpha = 0" >> tmp${i}.ctl
    echo "alpha = 0.5" >> tmp${i}.ctl
    echo "ncatG = 5" >> tmp${i}.ctl

    echo "tmp${i}.ctl after seds"
    cat tmp${i}.ctl

    #Run codeml
    /newhome/ab17362/paml4.9h/src/codeml tmp${i}.ctl
    echo "in.BV before concatination"
    cat in.BV

    cat rst2 >> in.BV
    mv rst2 rst2${i}

    echo "in.BV after concatination"
    cat in.BV
done


sed -i "s/usedata = 3/usedata = 2 in.BV/g" mcmctree.ctl
/newhome/ab17362/paml4.9h/src/mcmctree
