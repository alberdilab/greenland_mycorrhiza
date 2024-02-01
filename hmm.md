# HMM-based analysis pipeline

## Download UNITE database

## Align ITS sequences

```sh
#Create the batch file
cat <<EOF > msa.sh
#!/bin/bash
#SBATCH --job-name=0_subset
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=128gb
#SBATCH --time=6:00:00

module load raxml-ng/1.2.0
raxml-ng --search1 --msa sh_general_release_dynamic_25.07.2023.fasta --model GTR+G --threads 24
EOF

#Launch the batch file
sbatch msa.sh
```

## Create HMMs
