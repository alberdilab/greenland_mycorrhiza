# Greenland mycorrhiza
Bioinformatic respository of the shotgun metagenomic data analyses of the greenlandic mycorrhizal network project

This pipeline must be run in a high computation cluster. 

## Genome-resolved microbial metagenomics
The pipeline only relies on two software:

- Mamba
- Snakemake

The rest of the many required softwatre are downloaded and installed inside the conda environments.

### Prepare working environment

```sh
# Load dependencies mamba and snakemake (not needed if they are already installed in root)
module load mamba/1.5.6 snakemake/7.20.0

# Create working directory
mkdir greenland_shotgun

# Clone metagenomic assembly+binning pipeline repository
cd greenland_shotgun
git clone https://github.com/3d-omics/mg_assembly.git
cd mg_assembly

# Create screen session 
screen -S greenland_shotgun

# Create conda environments and run test data to validate them
# It might take 20-30 minutes to download and install all software
./run
```

### Prepare the input files

#### Data files
- Move the raw data files to the 'resources/reads/' directory or create soft links.

```sh
# In this example, download reads from ERDA
cd resources/reads/
wget https://sid.erda.dk/share_redirect/dcKtF82NjL/SH_IndPCR_1A_EKDL210009000-1a-AK4939-AK6653_HTF5CDSX2_L1_1.fq.gz
wget https://sid.erda.dk/share_redirect/dcKtF82NjL/SH_IndPCR_1A_EKDL210009000-1a-AK4939-AK6653_HTF5CDSX2_L1_2.fq.gz
cd ../../
```
