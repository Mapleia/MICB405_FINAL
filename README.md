# MICB405_FINAL
Measuring response to metal toxicity in Chlamydomonas.

## Setup
1. Install [R Language](https://www.r-project.org/) and [RStudio](https://rstudio.com/products/rstudio/download/).
2. Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
3. Install [Git](https://git-scm.com/downloads).
4. Maybe install [GitHub Desktop](https://desktop.github.com/).

1. `git clone https://github.com/Mapleia/MICB405_FINAL.git`


## To download SRAs

1. Install [SRA-Toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software).
2. Add the SRA-Toolkit to your environmental variable called PATH. See a tutorial [here for Windows](https://www.architectryan.com/2018/03/17/add-to-the-path-on-windows-10/).
3. Run [`./download_sra.sh`](./download_sra.sh) if you are in a bash shell, and [`download_sra.ps1`](./download_sra.ps1) if you are using powershell.

    > **Note:** you can change the list of SRA IDs by changing the [`.txt` file](./references/sra_id_list.txt) in the script.


## Copy the STAR parsed zip file
```bash
# scp [user@]SRC_HOST:]file1 [user@]DEST_HOST:]file2
user_name="your_user_name"
src_host="@orca1.bcgsc.ca"
archive="/projects/micb405/finalprojects2023W1/group11/alignments/PRJNA608616/Cd_exposure_STAR.zip"

dest_arch="Cd_exposure_STAR.zip"
scp "$user_name$src_host:$archive" $dest_arch
```