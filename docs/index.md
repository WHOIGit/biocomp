# BioComp: Bioinformatics protocols at WHOI

_Sharon Grim, December 2023_

This repository will hold documentation about using common 'omics analytical tools. Unless otherwise specified, these tools are installed on the **Poseidon HPC cluster**.

Access the Github Page here: https://whoigit.github.io/biocomp/

---

1. Available software on Poseidon
   * [Software via Singularity](https://github.com/WHOIGit/biocomp/wiki/Singularity-software-availability-20231027)
   * [RStudio on the HPC](https://github.com/WHOIGit/biocomp/wiki/Running-RStudio-on-the-HPC)

2. Common analytical workflows for MCG members
   * [DESeq for differential expression analysis](https://github.com/WHOIGit/biocomp/wiki/MCGomics-Bash-20230802)

3. Other bioinformatics protocols
  * [Submitting to NCBI's SRA](https://github.com/WHOIGit/biocomp/wiki/Uploading-to-NCBI-SRA-via-ASCP)
  * [Lab protocols](protocols.md)
    * [Apprill Lab MiSeq Library Preparation](protocols/ApprillLab_MiSeqLibraryPrep_SOP_Feb2024.pdf)
    * [Apprill Lab SOP for processing 16S rRNA gene amplicons through DADA2](protocols/DADA2_decontam_ApprillLabProtocol.html)

4. Code snippets
  * [Bandaging a metagenome-assembled-genomic bin](code-pages/snip-bandage.txt) (bash)
  * [Slurm job management, Bash profiles](https://github.com/WHOIGit/biocomp/wiki/MCGOmics-Bash-20231004)
  * [Finding duplicate files](https://github.com/WHOIGit/biocomp/wiki/Finding-your-duplicate-files-on-the-server-20231005)
  * [Parsing taxonomy of OTUs/ASVs](code-pages/taxonomy_parser.R) (R)
  * [Evaluating significant differential abundance of ASVs using Kruskal-Wallis testing](code-pages/kw_test.R) (R)
