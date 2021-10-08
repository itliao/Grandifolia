#!/bin/bash
#
#SBATCH --mem=230000
#SBATCH --cpus-per-task=5
#SBATCH --job-name=fs_test
#SBATCH -p rausherlab,common,scavenger
#SBATCH --account=rausherlab
#SBATCH --mail-user=itl4
#SBATCH --mail-type=ALL

#running fastSTRUCTURE, 10 separate runs with different starting seeds
while read i j
do
        
echo "Processing $i $j"

        for K in {2..10}
        do
        python2 /hpc/group/biodept/itl4/Programs/fastStructure/structure.py \
                -K $K \
                --input=/work/itl4/Austinii/admixture/plink6_210110 \
                --output=/work/itl4/Austinii/fastStructure/test$i \
                --full --seed=$j
        done
done < /hpc/group/biodept/itl4/RI_Austinii/randomNumbers.txt

#for choosing complexity from each run
python /hpc/group/biodept/itl4/Programs/fastStructure/chooseK.py \
	--input=/work/itl4/Austinii/fastStructure/test 