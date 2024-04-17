# Computational meta'omics workflow examples

- Author: Sharon Grim
- Version: 0.1 
- Date: April 2024

---

*Do you have sequenced (meta')omic samples, and don't know where to start? Do you need inspiration for your workflow? Are you a HPC user, or have access to a computing cluster?*

*You've come to (one of) the right place(s)!*

This is an example workflow I used for processing and analyzing metagenomic samples in March and April 2024. Example scripts will be available as noted, and if you have access to the HPC, they will be in the ```/proj/omics/bioinfo``` space as well. Most of these tools are invoked via Singularity containers on the HPC, which are explicitly loaded in the wrapper scripts; if you choose to go off-cuff, please remember to wrap the job in either an ```srun``` or ```sbatch``` command.

---

1. Raw read QC through FastP
2. Metagenomic assembly through Megahit
3. Choose Your Own Adventure
    - 3A. Making MAG bins
      - i. index assembly and map samples
      - ii. use metabat
      - iii.  assess bin quality and completion
    - 3B. Gene calling
      - i. Prokka
      - ii. Prodigal
    
---

### 1. FastP to trim your reads.
   
After you get your raw sequences, the first step in an 'omics workflow is quality control of your sequenced reads. I implemented QC and trimming of raw fastq files through [FastP](https://github.com/OpenGene/fastp), which I invoke on the HPC through a Singularity container. I wrapped the trimming and QC in a shell script (mcgomics_fastp.sh), and submitted the job to the Slurm scheduler:

```sbatch mcgomics_fastp.sh (input directory path* (output directory path) (optional text file of sample names)```

For this example, I can submit the job from anywhere in my space on the server, but I can change directories to my working space and use relative paths thereafter by using the ```--chdir``` flag prior to the submission script. (You can also override any resource requests in the shell script, by including the appropriate flags in between sbatch and the name of the script. See my below examples where I gradually lengthen the command line submission with these flags...)

```
sbatch --chdir=$SCRATCH/projects/ /proj/omics/bioinfo/scripts/slurm/mcgomics_fastp.sh \
	$SCRATCH/project/fastqs \
	$SCRATCH/project/fastp \
	$SCRATCH/project/samplenames.txt
```

The samplenames.txt is a text file with one sample name per line, corresponding to the name of the FastQ file(s) (paired reads). 

Output from this wrapper was:
- trimmed reads (for forward and reverse sets, two fastq.gz files)
- a trim report in html format
- a trim report in JSON format

Planning resource allocation, I used 3GB of memory, 8 CPUs, and 45min of job time for a sample that was 20M reads of 151bp each.

---

### 2. MEGAHIT to assemble your metagenome

After QC and trimming, (meta)genomic reads get assembled. I like to use [MEGAHIT](https://github.com/voutcn/megahit/). Just as with FastP, I wrapped the asssembly step in a shell script to submit as a Slurm job on our HPC:

```sbatch mcgomics_megahit.sh (output directory) (directory with your trimmed reads) (optional: extra flags to pass to MEGAHIT such as "--continue")```

That output from FastP can go right into this assembly step. You can also use symbolic links (symlinks) in a new directory, if you want to separate out your trimmed reads for different assemblies. 

Optional flags include "--continue" if you want to pick up an assembly you already started. 

As before, I'm changing directory to my working directory before running the script, but in this example my arguments are provided *relative to that working directory*.

```
sbatch --chdir=$SCRATCH/projects /proj/omics/bioinfo/scripts/slurm/mcgomics_megahit.sh \
 ./coassembly \
 ./fastp \
 "--continue"
```

The MEGAHIT outputs will be a "final.contigs.fa" file (if assembly finished) and a checkpoints text file (for multiple kmer assemblies), some log information, etc.

---

### 3. Choose Your Own Adventure

After assembling, you have some choices of next steps. These are not exclusive, but you can make life easier downstream if you kind of know what you want to do now before you continue.

To call genes via Prokka or Prodigal, go to step 3B.

To make metagenome-assembled-genome bins from your assembly, go to step 3A.

To map metatranscriptome samples to your assembled metagnome, go to page XYZ. *not available yet, check later!*

---

### 3A. Make metagenome-assembled-genomic (MAG) bins

#### i. Index and map your samples to the assembly

In lieu of cultured representatives and their genomes, representative organisms from natural assemblages can be approximated through metagenome-assembled-genomic (MAG) bins from your metagenome assembly.  A standard way to create those is to use differential coverage and tetranucleotide frequencies in the metagenome. 

For this approach, first you need to map your sequences to the metagenomic assembly to get differential coverage (BAM) files. I use [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer) to index the assembly and map samples. *(You can use another index and mapping software like bwa-2, nothing against it, but I didn't write wrappers for those.)* Bowtie2 is installed on the Poseidon HPC via a module, which gets invoked in this shell script for indexing your assembly:

```sbatch mcgomics_bowtie2_indexing.sh (full path to assembled contigs file) (prefix for your indexed assembly) (directory for output of indexing) (optional: extra arugments such as "--large-index")```

and for mapping reads to that assembly:

```sbatch mcgomics_bowtie2_mapping.sh (full path to the directory with assembly index) (prefix for your indexed assembly) (output directory for mapping files) (full path to trimmed reads) (optional: extra arguments such as "-q")```

The (prefix) needs to be consistent. The directory of trimmed reads needs to have filenames (or symlinks) with "trim", and "rev" or "fwd", to be compliant with my argument parser. Not very glorious, but the BAM output from bowtie2 has to be sorted and indexed. I used samtools for that step, which is called via a module in this wrapper.

Examples:

```
sbatch --chdir=$SCRATCH/projects --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=12gb --partition=compute /proj/omics/bioinfo/scripts/slurm/mcgomics_bowtie2_indexing.sh \
 ./coassembly/final.contigs.fa \
 "sample" \
 ./coassembly/
```

Here I made a new directory, and within symlink'ed the trimmed reads I wanted to map to this assembly:

```
sbatch --chdir=$SCRATCH/projects --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=20gb --partition=compute /proj/omics/bioinfo/scripts/slurm/mcgomics_bowtie2_mapping.sh \
 ./coassembly/bowtie2-index \
 "sample" \
 ./coassembly/coassembly_mapping/${i} \
 ./mapping/${i}
```

Your output files are:
- build.log
- (prefix).0-9.bt2
- (prefix).rev.0-9.bt2
- (sample).sorted.bam
- (sample).sorted.bam.bai
- (sample).bowtie2.log
- (sample).name_sorted.bam

---

#### ii. Use metabat binning software

After getting those coverage files (sample).sorted.bam and (sample).sorted.bam.bai, you can use a binning software such as [metabat](https://bitbucket.org/berkeleylab/metabat/) to generate MAG bins. Here is the wrapper I wrote:

```sbatch mcgomics_metabat_extra.sh (path to assembly) (output directory for your bins) (path to your BAM files) (prefix for your bins) (optional: extra arguments like "--noBinOut")```

I made a clean directory with the mapping files symlinked within, in this example:

```
sbatch --chdir=$SCRATCH/projects/coassembly/coassembly_metabat --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=12Gb --time=8:00:00 --partition=compute /proj/omics/bioinfo/scripts/slurm/mcgomics_metabat_extra.sh \
 $SCRATCH/projects/coassembly/final.contigs.fa \
 $SCRATCH/projects/coassembly/coassembly_metabat/ \
 $SCRATCH/projects/coassembly/bam_files \
 "coassembly"
```

---

#### iii. Assess quality and completion of your MAG bins.

What to do with your bins (or your assembled genome, if it's not a metagenome)? Assess quality and phylogeny via tools like [checkm](https://github.com/Ecogenomics/CheckM/) (implemented in [dRep](https://github.com/MrOlm/drep)), [dastool](https://github.com/cmks/DAS_Tool), or [CheckM2](https://github.com/chklovski/CheckM2).

*CheckM2*

```sbatch mcgomics_checkm2.sh (output directory) (directory with your binned contigs) (optional but recommended: arguments such as "--extension=fa,--force")```

```
sbatch --chdir=$SCRATCH/projects --nodes=1 --ntasks=1 --cpus-per-task=12 --mem=32Gb --time=6:00:00 --partition=compute /proj/omics/bioinfo/scripts/slurm/mcgomics_checkm2.sh \
 ./coassembly/checkm2/ \
 ./coassembly/coassembly_metabat/ "--extension=fa,--force"
```

*dRep*

```sbatch mcgomics_drep.sh (output directory) (directory with your binned contigs) (optional: arguments such as "--skip_plots,--ignoreGenomeQuality" which skip checkm)```

```
sbatch --chdir=$SCRATCH/projects --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=24Gb --time=24:00:00 --partition=compute /proj/omics/bioinfo/scripts/slurm/mcgomics_drep.sh \
 $SCRATCH/projects/coassembly/drep \
 $SCRATCH/projects/coassembly/coassembly_metabat \
 "--skip_plots,--ignoreGenomeQuality"
```

*DAS Tool*

I didn't get around to wrapping arguments for dastool, nor did I try it out much, but you can invoke it via singularity all the same. Start a job via ```srun``` or ```sbatch```, then load singularity as a module.

```
module load singularity

singularity shell \
 --cleanenv \
 --bind="/proj/omics/,/scratch/${USER},/user/${USER}" \
 /proj/omics/bioinfo/databases/nfx_singularity_cache/dastool_1.1.7.sif
```

---

### 3B. Call genes

What if you want to look at genes in your assembly? We can use [Prokka](https://github.com/tseemann/prokka) like so:

```sbatch mcgomics_prokka.sh (output directory) (input contigs file) (kingdom for gene calling) (prefix for your results) (optional: provide the training files for gene calling)```


#### i. Prokka

```
sbatch --chdir=$SCRATCH/projects/coassembly/ --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=84Gb --time=128:00:00 --partition=compute --qos=unlim /proj/omics/bioinfo/scripts/slurm/mcgomics_prokka.sh \
 ./prokka \
 coassembly.contigs.fa \
 Bacteria \
 sample_coassembly
```

Outputs from Prokka are:
- (sample_coassembly).fna *is the contigs sequences*
- (sample_coassembly).fsa *is contigs sequences with some flavortext [gcode=11] [organism=Genus species] [strain=strain]*
- (sample_coassembly).gff *is list of contigs with the gene calls, and contig sequences*
- (sample_coassembly).tbl *is feature table with start, stop, CDS, and prediction for each gene*
- (sample_coassembly).tsv *is list of genes called with product information*
- (sample_coassembly).sqn *is a large json-like formatted file with genes, sequences, etc.*

#### ii. Prodigal

Or if you're like the majority of Prokka users and encounter unexplicable slowdowns at the tbl2asn step... call genes quickly with [prodigal](https://github.com/hyattpd/prodigal/wiki/cheat-sheet). No wrapper for this yet either but after you invoke a job ```srun``` or ```sbatch``` here's a start:

```
module load default-environment
module load bio
module load prodigal

prodigal -i (path to your input fasta file) -p meta -f gff -g 11 &> prodigal.log

```

---

Reach out if you have any questions or suggestions!

*[Sharon Grim](sharon.grim@whoi.edu), 2024 April*
