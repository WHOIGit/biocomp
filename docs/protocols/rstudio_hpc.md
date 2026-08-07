# Running Interactive RStudio from the HPC

*20260807 Sharon Grim*

>Note: This protocol was written for using RStudio containers supported by Sharon at WHOI HPC located at `/proj/omics/bioinfo/`.
>
>This does not necessarily apply to the HPC installation of R (see `module avail`).
>
>Currently Apptainer and Singularity are interchangeable on the server, don't be alarmed if you type `singularity shell...` and get `Apptainer>` in response.
>
>If you go off script, please take time to research your path before reaching out to Sharon.* 

So you want to run (interactive RStudio in an) containerized R on the HPC to use the computational resources? Going forward, "container RStudio" refers to the Apptainer/Docker image of RStudio that comes from [Rocker Project](https://rocker-project.org/). This is different from your local installation on your computer, and different from the HPC installation of R. If you're still reading because this applies to you, please look below at these tips.

1. Using SLURM to run your job
2. Version controlling your packages
3. Command line installation of R packages
4. Other options

In this document, I have provided the command lines `within these highlights` unless it's a substantial chunk and not an elegant one-liner. In which case, they have blocks for code.

When I can, I provided the command feedback *as I went through the process* within blocks like below:
```
'getOption("repos")' replaces Bioconductor standard repositories, see
'help("repositories", package = "BiocManager")' for details.
Replacement repositories:
    CRAN: https://p3m.dev/cran/__linux__/noble/latest
Bioconductor version 3.23 (BiocManager 1.30.27), R 4.6.1 (2026-06-24)
Installing package(s) 'phyloseq'
trying URL 'https://bioconductor.org/packages/3.23/bioc/src/contrib/phyloseq_1.56.0.tar.gz'
Content type 'application/x-gzip' length 6387174 bytes (6.1 MB)

```

---

## Using SLURM to run your job

I've provided example slurm scripts to prepare and run container RStudio: `/proj/omics/bioinfo/scripts/slurm/`. Specifically, what you want to keep in mind is *which version of R you want* because these scripts are specific to versions:


`ls -1 /proj/omics/bioinfo/scripts/slurm/singularity_launch_rstudio_*.sh`
```
/proj/omics/bioinfo/scripts/slurm/singularity_launch_rstudio_20240624.sh
/proj/omics/bioinfo/scripts/slurm/singularity_launch_rstudio_20250728.sh
/proj/omics/bioinfo/scripts/slurm/singularity_launch_rstudio_20251104.sh
/proj/omics/bioinfo/scripts/slurm/singularity_launch_rstudio_20260131.sh
/proj/omics/bioinfo/scripts/slurm/singularity_launch_rstudio_20260310.sh
/proj/omics/bioinfo/scripts/slurm/singularity_launch_rstudio_20260807.sh
```

`diff --side-by-side /proj/omics/bioinfo/scripts/slurm/singularity_launch_rstudio_20260807.sh \
 /proj/omics/bioinfo/scripts/slurm/singularity_launch_rstudio_20260310.sh | head -n 20 | tail --lines=+12`

```
#20260807						                                          <
# Rstudio upgrade from R v. 4.6.1			                        <
# make sure we can install from source any packages you need  <
							                                                <
#20260310							                                        #20260310
#blizzard-caused server failure, now Hydra does not have pyth	#blizzard-caused server failure, now Hydra does not have pyth
#if 'python' not found, use python3				                    #if 'python' not found, use python3

```

Briefly, `singularity_launch_rstudio_20260807.sh` uses the container RStudio based on R v.4.6.1, whereas `singularity_launch_rstudio_20260310.sh` uses R v.4.5.2. 

>(`singularity_launch_rstudio_20260131.sh` also uses R v.4.5.2, but there was an error invoking python between January and March 2026 so that's why I patched it.)

---

## Version controlling your packages

If you use these slurm scripts, by default you will install R packages into /user/$USER/R/rocker-rstudio/ so as to not complicate the HPC's version of R. This process also makes use of /scratch/$USER/ so make sure you have access to that.

But sometimes even 'supported' packages installed via `install.packages()` or `BiocManager::install()` (*I'm looking at you, [Spiec-Easi](github.com/zdk123/SpiecEasi)*), can mess up your other packages. If you encounter issues when installing new packages, here is what I recommend:

1. **If you're upgrading from a base.major.minor to base.major.minor+1 such as from  4.5.1 to 4.5.2:** Make a backup of your currrent (even if it's broken) package directory. In this example, let's try to triage 4.5.2 when [phyloseq](https://joey711.github.io/phyloseq/) messed up [igraph](https://r.igraph.org/):

```
mkdir -p ~/R/rocker-rstudio/backup_4.5.2

sbatch --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=8gb --time=4:00:00 --partition=scavenger --qos=scavenger --wrap='rsync -zav ~/R/rocker-rstudio/4.5.2/* ~/R/rocker-rstudio/backup_4.5.2/'
```

2. **If you are using a different minor or major version of R than what you had originally,** make sure to edit your Renviron file accordingly. In this example, I'm going from R v.4.5.2 to v.4.6.1 so my Renviron file is edited as such:

 
`mv /scratch/$USER/rstudio-server/Renviron.site /scratch/$USER/rstudio-server/backup_Renviron.site`
`vi /scratch/$USER/rstudio-server/Renviron.site` (or `nano` if that's your flavor)
```
BIOCONDUCTOR_USE_CONTAINER_REPOSITORY=FALSE
OMP_NUM_THREADS=2
SINGULARITY_SLURM_JOB_CPUS_PER_NODE=2
SINGULARITY_SLURM_JOB_NUM_NODES=1
SINGULARITY_SLURM_TASKS_PER_NODE=1
SINGULARITY_RINGDIR="/scratch/${USER}/rstudio-server"
SINGULARITY_BIND="/proj/omics/,/scratch/${USER}:${SCRATCH},/user/${USER}:${HOME}"
SINGULARITYENV_R_PARALLELLY_AVAILABLECORES_FALLBACK=1
SINGULARITYENV_R_LIBS="${HOME}/R/rocker-rstudio/4.6.1"
SINGULARITYENV_R_LIBS_USER="${HOME}/R/rocker-rstudio/4.6.1"
SINGULARITYENV_R_ENVIRON_USER="${SCRATCH}/rstudio-server/Renviron.site"
```

3. Whether you had to backup your packages *(using same base.major.minor)* or just edit your Renviron file *(new version)*, here's where the workflow is the same. Start an interactive job and go into R (*not RStudio*) within the container. The below example uses R version 4.6.1, but edit it to reflect the version you are intending to use:

```
singularity exec \
--cleanenv \
--env-file "/scratch/${USER}/rstudio-server/Renviron.site" \
--env='R_LIBS=${HOME}/R/rocker-rstudio/4.6.1' \
--env='R_LIBS_USER=${HOME}/R/rocker-rstudio/4.6.1' \
--env='R_PARALLELLY_AVAILABLECORES_FALLBACK=1' \
--bind="${HOME}/R/rocker-rstudio/4.6.1:/usr/local/lib/R/host-site-library" \
--bind="/proj/omics/,/scratch/${USER}:${SCRATCH},/user/${USER}:${HOME}" \
/proj/omics/bioinfo/databases/nfx_singularity_cache/rstudio_rocker.4.6.1.sif R
```

You should be welcomed with something like this screen:
```
R version 4.6.1 (2026-06-24) -- "Happy Hop"
Copyright (C) 2026 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> 
```

4. Within R console, you can use typical R commands such as `install.packages()` or `BiocManager::install()`. But if those generate errors like below... keep reading.

`install.packages("igraph")`
```
Installing package into ‘/user/sharon.grim/R/rocker-rstudio/4.6.1’
(as ‘lib’ is unspecified)
trying URL 'https://p3m.dev/cran/__linux__/noble/latest/src/contrib/igraph_2.3.3.tar.gz'
Content type 'binary/octet-stream' length 6246009 bytes (6.0 MB)
==================================================
downloaded 6.0 MB

* installing *binary* package ‘igraph’ ...
* package ‘igraph’ successfully unpacked and SHA256 sums checked
* DONE (igraph)

The downloaded source packages are in
	‘/tmp/Rtmp1FlMfT/downloaded_packages’

```

`library(igraph)`
```
Error: package or namespace load failed for ‘igraph’ in dyn.load(file, DLLpath = DLLpath, ...):
 unable to load shared object '/user/sharon.grim/R/rocker-rstudio/4.6.1/igraph/libs/igraph.so':
  libglpk.so.40: cannot open shared object file: No such file or directory

```

---

## Command line installation of R packages

Great! Now R won't load a package that you saw was installed, within a container version of RStudio that you didn't even install, and you get this cryptic message.

What you need to do is **quit** ... this instance of R console (`q()`) because we're going to be invoking R commands *from the command line*.

```
singularity shell \
--cleanenv \
--env-file /scratch/sharon.grim/rstudio-server/Renviron.site \
--env='R_LIBS=${HOME}/R/rocker-rstudio/4.6.1' \
--env='R_LIBS_USER=${HOME}/R/rocker-rstudio/4.6.1' \
--env='R_PARALLELLY_AVAILABLECORES_FALLBACK=1' \
--bind="${HOME}/R/rocker-rstudio/4.6.1:/usr/local/lib/R/host-site-library" \
--bind="/proj/omics/,/scratch/${USER}:${SCRATCH},/user/${USER}:${HOME}" \
/proj/omics/bioinfo/databases/nfx_singularity_cache/rstudio_rocker.4.6.1.sif
```

You'll see this as a result:

`Apptainer> `

To install packages *from source* to your R_LIBS directory within this container, here are some examples:

`R -q -e 'install.packages("remotes")'`
```
> install.packages("remotes")
Installing package into ‘/user/sharon.grim/R/rocker-rstudio/4.6.1’
(as ‘lib’ is unspecified)
trying URL 'https://p3m.dev/cran/__linux__/noble/latest/src/contrib/remotes_2.5.0.tar.gz'
Content type 'binary/octet-stream' length 440847 bytes (430 KB)
==================================================
downloaded 430 KB

* installing *binary* package ‘remotes’ ...
* package ‘remotes’ successfully unpacked and SHA256 sums checked
* DONE (remotes)

The downloaded source packages are in
	‘/tmp/RtmpngvCOE/downloaded_packages’
```

Specifically, because **phyloseq** will try to update **igraph** from binary, and RStudio within the container will not permit that (see [this error](https://r.igraph.org/articles/installation-troubleshooting#libglpk-so-40-cannot-open-shared-object-file-no-such-file-or-directory)), we need to install from source before going into RStudio. The [archive of available igraph versions](https://cran.r-project.org/src/contrib/Archive/igraph/) has v.2.3.2 as the most recent, so we use that link as below:

`R -q -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.3.2.tar.gz", quiet = TRUE, lib = .libPaths()[1])'`
```
> install.packages("https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.3.2.tar.gz", quiet = TRUE, lib = .libPaths()[1])
inferring 'repos = NULL' from 'pkgs'
trying URL 'https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.3.2.tar.gz'
Content type 'application/x-gzip' length 5209216 bytes (5.0 MB)
==================================================
downloaded 5.0 MB

```

Now we have an installation of **igraph** that should work within R in container RStudio. You can check that by `singularity exec` as above, and within R console type `library(igraph)`.

---

When troubleshooting the installation of **phyloseq**, I saw that some packages were trying to be updated from binary (such as **igraph**) but others should not have had a problem. 

Let's install some of the [requirements](https://github.com/joey711/phyloseq/blob/master/DESCRIPTION) for phyloseq before doing the main package install. You can compare the requirements against your current available packages:

`ls -1 ~/R/rocker-rstudio/4.6.1/` *(or whichever version you're using)*

I saw that I didn't have **data.table** and **cluster** already, so in preparation for **phyloseq** I installed those packages from CRAN.

`R -q -e 'install.packages(c("cluster", "data.table"), lib = .libPaths()[1])'`

```
> install.packages(c("cluster", "data.table"), lib = .libPaths()[1])
trying URL 'https://p3m.dev/cran/__linux__/noble/latest/src/contrib/cluster_2.1.8.3.tar.gz'
trying URL 'https://p3m.dev/cran/__linux__/noble/latest/src/contrib/data.table_1.18.4.tar.gz'
* installing *binary* package ‘cluster’ ...
* package ‘cluster’ successfully unpacked and SHA256 sums checked
* DONE (cluster)
* installing *binary* package ‘data.table’ ...
* package ‘data.table’ successfully unpacked and SHA256 sums checked
* DONE (data.table)

The downloaded source packages are in
	‘/tmp/RtmpdLy9Zm/downloaded_packages’
> 
```

As before, you can check within the container R console that these can be attached. 

---

Now let's install phyloseq:

`R -q -e 'BiocManager::install(c("phyloseq"), lib = .libPaths()[1], ask = TRUE, update = FALSE)'`
```
>> BiocManager::install(c("phyloseq", "zdk123/SpiecEasi"), lib = .libPaths()[1], ask = TRUE)
'getOption("repos")' replaces Bioconductor standard repositories, see
'help("repositories", package = "BiocManager")' for details.
Replacement repositories:
    CRAN: https://p3m.dev/cran/__linux__/noble/latest
Bioconductor version 3.23 (BiocManager 1.30.27), R 4.6.1 (2026-06-24)
Installing package(s) 'phyloseq'
trying URL 'https://bioconductor.org/packages/3.23/bioc/src/contrib/phyloseq_1.56.0.tar.gz'
Content type 'application/x-gzip' length 6387174 bytes (6.1 MB)
==================================================
downloaded 6.1 MB

* installing *source* package ‘phyloseq’ ...
** this is package ‘phyloseq’ version ‘1.56.0’
** using staged installation
** R
** data
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (phyloseq)

The downloaded source packages are in
	‘/tmp/RtmpZE4rVi/downloaded_packages’
> 
```

Then install **SpiecEasi** from Bioconductor:

`R -q -e 'BiocManager::install(c("SpiecEasi"), lib = .libPaths()[1], ask = TRUE, update = FALSE)'`
```
> BiocManager::install(c("SpiecEasi"), lib = .libPaths()[1], ask = TRUE, update = FALSE)
'getOption("repos")' replaces Bioconductor standard repositories, see
'help("repositories", package = "BiocManager")' for details.
Replacement repositories:
    CRAN: https://p3m.dev/cran/__linux__/noble/latest
Bioconductor version 3.23 (BiocManager 1.30.27), R 4.6.1 (2026-06-24)
Installing package(s) 'SpiecEasi'
trying URL 'https://bioconductor.org/packages/3.23/bioc/src/contrib/SpiecEasi_2.0.0.tar.gz'
...
** R
** data
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** checking absolute paths in shared objects and dynamic libraries
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (SpiecEasi)

The downloaded source packages are in
	‘/tmp/RtmpJ3hvF8/downloaded_packages’
> 
```

All looks good so far. Doublecheck it by going into the R Console within the container (`singularity exec ... R`) and attaching the packages of interest:
`library(phyloseq)`
`library(SpiecEasi)`

---

## Other options

If you encounter issues with installing/ugprading packages within container RStudio, so long as the the above options for package installation are working (R console within the container, command line within the container), you should be able to navigate and troubleshoot those errors.

For example, by installing `remotes` in R we can use this to install any $TOOL at specific version ranges:

`R -q -e 'remotes::install_version($TOOL, version = ">= 1.0.1, < 2.2.1", lib = .libPaths()[1])'`

Or if you want to install from source or github:

`remotes::install_github("zdk123/SpiecEasi")`

---

>Do you have any follow-up or weird cases to discuss? Reach out via email, Slack, or even add an issue to Github.


