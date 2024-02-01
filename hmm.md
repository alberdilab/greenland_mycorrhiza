# HMM-based analysis pipeline

## Download UNITE database

## Align ITS sequences

```sh
#Create the batch file
cat <<EOF > msa.sh
#!/bin/bash
#SBATCH --job-name=msa
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=128gb
#SBATCH --time=6:00:00

module load mafft/7.515
mafft --retree 1 --thread 24 sh_general_release_dynamic_25.07.2023.fasta > unite_align.fasta
EOF

#Launch the batch file
sbatch msa.sh
```

## Create HMMs

```sh
#Create the batch file
cat <<EOF > buildhmm.sh
#!/bin/bash
#SBATCH --job-name=0_subset
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=128gb
#SBATCH --time=6:00:00

module load hmmer/3.3.2
hmmbuild globins4.hmm globins4.sto
EOF

#Launch the batch file
sbatch msa.sh
```
