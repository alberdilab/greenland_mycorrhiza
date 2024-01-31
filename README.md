# Greenland mycorrhiza
Bioinformatic respository of the shotgun metagenomic data analyses of the greenlandic mycorrhizal network project

The following genome-resolved metagenomic pipeline is run: https://github.com/3d-omics/mg_assembly

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
screen -S greenland_shotgun

# Create conda environments and run test data to validate them
# It might take 20-30 minutes to download and install all software
./run
```

## Prepare the input files

### Data files
- Move the raw data files to the 'resources/reads/' directory or create soft links.

```sh
# In this example, download reads from ERDA
cd resources/reads/
wget https://sid.erda.dk/share_redirect/dcKtF82NjL/SH_IndPCR_1A_EKDL210009000-1a-AK4939-AK6653_HTF5CDSX2_L1_1.fq.gz
wget https://sid.erda.dk/share_redirect/dcKtF82NjL/SH_IndPCR_1A_EKDL210009000-1a-AK4939-AK6653_HTF5CDSX2_L1_2.fq.gz
cd ../../
```
  
- Download the reference genome of the host to the 'resources/reference/' directory.

```sh
# In this example, download the genome of Carex myosuroides from Genbank
cd resources/reference/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/015/225/GCA_028015225.1_CAF_CaMyos_1.0/GCA_028015225.1_CAF_CaMyos_1.0_genomic.fna.gz
cd ../../
```
### Sample input file
Edit the 'config/samples.tsv' document and add samples names, library identifiers, file locations, adapter sequences and (co)assembly information.

- **sample_id:** unique name that will be used to identify files associated with each sample. Ideally, these should be short and contain only letters and numbers to avoid hiccups during the data processing.
- **library_id:** names of the different libraries that might be associated with each sample. If each sample has a single library, using 'lib1' in all cases will make the trick. This name will be appended to the sample_id.
- **forward_filename:** relative path to the forward file. Usually 'resources/reads/' and the forward file name.
- **reverse_filename:** relative path to the forward file. Usually 'resources/reads/' and the reverse file name.
- **forward_adapter:** adapter sequence to be used in the quality-filtering step. Usually the Illumina adapter 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'.
- **reverse_adapter:** adapter sequence to be used in the quality-filtering step. Usually the Illumina adapter 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'.
- **assembly_ids:** comma-separated list of (co)assembly groups.

```tsv
sample_id	library_id	forward_filename	reverse_filename	forward_adapter	reverse_adapter	assembly_ids
M10	lib1	resources/reads/sSH_IndPCR_1A_EKDL210009000-1a-AK4939-AK6653_HTF5CDSX2_L1_1.fq.gz	resources/reads/SH_IndPCR_1A_EKDL210009000-1a-AK4939-AK6653_HTF5CDSX2_L1_2.fq.gz	AGATCGGAAGAGCACACGTCTGAACTCCAGTCA	AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT	M10
```

### Add host reference genome to the config file
Edit 'config/features.yml' with reference databases:

```yaml
host:
  fasta: resources/reference/GCA_028015225.1_CAF_CaMyos_1.0_genomic.fna.gz #this is the line that has been changed

magscot:
  pfam_hmm: workflow/scripts/MAGScoT/hmm/gtdbtk_rel207_Pfam-A.hmm.gz
  tigr_hmm: workflow/scripts/MAGScoT/hmm/gtdbtk_rel207_tigrfam.hmm.gz

dram_database: "resources/mock_dram_db"
gtdbtk_database: "resources/mock_gtdbtk_db"
singlem_database: "resources/mock_singlem_db"
kraken2_database: "resources/kraken2_mock"
```

## Run pipeline
Once all the input files and data are properly set, it is time to launch the pipeline.
```sh
# Ensure the pipeline is launched from the working directory
cd greenland_shotgun/mg_assembly

# Resume the screen session if you closed it
screen -r greenland_shotgun

# Load mamba and snakemake dependencies if needed
module load mamba/1.5.6 snakemake/7.20.0

# Run the pipeline using slurm
./run_slurm
```
