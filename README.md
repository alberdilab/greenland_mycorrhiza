# Greenland mycorrhiza
Bioinformatic respository of the shotgun metagenomic data analyses of the greenlandic mycorrhizal network project

This pipeline must be run in a high computation cluster.

## Root software dependencies
The pipeline only relies on two software:

- Mamba
- Snakemake

The rest of the many required softwatre are downloaded and installed inside the conda environments.

## Prepare working environment

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
screen -S

# Create conda environments and run test data to validate them
./run
```

## Prepare the input files

### Data files

### Sample input file

