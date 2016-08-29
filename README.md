
# BARTI (Barcode Tracking In vivo)
## Pipeline for barcode sequence data analysis

<Hello - welcome to our Github home! Please feel free to peruse our project and share our interest in tuberculosis and infectious disease, as well as the joy of biology and CS. If you do use some of the code or other material from this project, please be sure to cite us (see the reference in the overview). Thanks! :)  Scripts written by Vivian W. Leung (1). 2016.>


Version 1.0 (stable, updated 26 August 2016) on <i>master</i> branch.
(<i>v1</i> and <i>v2</i> branches are in development.)



(Published (XX):
Digitally Barcoding Mycobacterium tuberculosis reveals in vivo infection dynamics in the macaque model of tuberculosis)

<extrasmall>Constance J. Martin (1, 2 5), Anthony M. Cadena (3, 5), Vivian W. Leung (1), Philana Ling Lin (4), Pauline Maiello (3), Nathan Hicks (1), Michael R. Chase (1), JoAnne L. Flynn (3) & Sarah M. Fortune (1, 2)<extrasmall>
<br/><br/>

<extrasmall> <small>(1) Department of Immunology and Infectious Diseases, Harvard T.H. Chan School of Public Health, Boston, Massachusetts, USA. (2) Ragon Institute of Massachusetts General Hospital, Massachusetts Institute of Technology and Harvard, Cambridge, Massachusetts, USA. (3) Department of Microbiology and Molecular Genetics, University of Pittsburgh School of Medicine, Pittsburgh, Pennsylvania, USA. (4) Department of Pediatrics, Children’s Hospital of Pittsburgh, University of Pittsburgh Medical Center, Pittsburgh, Pennsylvania, USA. (5) These authors contributed equally to this work. Correspondence should be addressed to J.L.F (joanne@pitt.edu) or S.M.F. (sfortune@hsph.harvard.edu).<small> </extrasmall>

-------------
## Overview
------------
### Dependencies


<tl;dr>For those more comfortable with programming and Python, skip to Dependencies/Packages.

<If you have Linux, you probably know what you are doing anyway>

##### I.  Basic Needs (Command-Line Interface (CLI) and Python)</b>
0. (For Windows users only). Download a command-line interface program (CLI). Recommended: [Cygwin (cygwin.com)](https://cygwin.com).

0. Install Python v2 (ver. 2.7.11 or later; not v3). (If you have OS X) Recommended: [Anaconda (MiniConda version; anaconda.org)](https://anaconda.org).

0. pip????

##### II.  Packages

The following packages are required (versions in parentheses).

* numpy (1.10)
* pandas (0.18)
* sqlalchemy
* regex  [see Note 1]
* jupyter (1.0.0)  [see Note 2]

<b>Note 1 (regex)</b>: this <code>regex</code> is NOT the native <code>re</code> package. Install <code>regex</code> via <code>binstar</code>. The package is available as conda for OS X and Linux; for Windows, only the pip type is available.

<b>Note 2 (jupyter)</b>: Currently, Jupyter is required to read and execute the scripts (.ipynb type) while we test out our command-line .py version. <Sorry!>

##### (Quick Notes)
<b>How to install packages</b>

(If you're new to using command-lines and want more detail, read the expanded instructions in "How to install packages".)

In your CLI (e.g. Terminal or Cygwin)

Option A (preferred): to install with donda, type:
> $ <code> conda install [your_package_name]</code> <br/>
Example: <code> conda install my_new_package</code> <br/>...<br/>
> <code>Proceed ([y]/n)?</code> (blinking cursor here) # type <code>y</code> to continue, <code>n</code> to cancel, and <code>Enter</code> to submit.

Option B: to install with pip (if conda doesn't work or package is unavailable on conda).

> <code> pip install [your_package_name]</code> <br/>
Example: <code> pip install conda_cant_find_me</code> <br/>

> If already installed, <br/>
> <code> pip install [your_package_name] --upgrade </code> <br/>
Example: <code> conda install im_already_here --upgrade</code> <br/>

-------
### Downloading the pipeline.

#### <CAN WE COMPILE A PIPELINE?????>

#### Option 1:  Download the scripts [here](https://github.com/sarahfortunelab/barcodetracking/archive/master.zip).

Or, on the [main repository page](https://github.com/sarahfortunelab/barcodetracking), select the green "Clone or download" button. A dropdown menu will appear, and select the "Download ZIP" option.

Open and extract the files in the downloaded zip, and ensure they are in your desired folder (so you may easily access your data from the scripts, i.e. file paths)

#### Option 2:  Clone repository to your Github.

If you would like to get the project in an existing folder, initialize the folder first with <code>$ git init</code>. (If not, skip this step.)

Then, clone and add remote for your repository.

><code>$ git clone https://github.com/sarahfortunelab/barcodetracking.git<br/>$ git remote add [your_remote_name]</code>, conventionally named  <code>origin</code> for the first remote.


#### Option 3:  Fork repository to also stay up-to-date with our developments.

On the Github repository page, click the button <b>Fork</b> in the top-right corner to create a fork of the original project.

Then, navigate to forked Git repository, click <b>Clone or download</b>, and copy the Git address.

Clone your repository to your local folder with the following, wherein "YOUR_USER_NAME" is where your own user name will be:

<code>$ git clone https://github.com/YOUR_USER_NAME/barcodetracking.git</code>

To sync your repository with our repository to get updates, add an upstream remote to our repository:

<code>$ git add remote upstream https://github.com/sarahfortunelab/barcodetracking.git</code>

To receive updates to your local branch, e.g. <code>master</code> through your remote <code>upstream</code>:

<code>$ git fetch upstream<br/></code>
<code>$ git checkout [local branch to sync]</code>, typically named the same as that of the main repository<br/>
<code>$ git merge upstream/master</code><br/>






--------
## Comments on Installation
----------
### Dependencies

##### I. Basic Needs
1. <b>Command-line interface.</b>
   * <b>OS X</b>. The application "Terminal" is factory-installed on your computer as "Terminal". (Skip this step).
   * <b>Windows</b>. A Linux-like terminal is strongly preferred and is unfortunately not included in the Windows system. A good application is  [Cygwin (https://cygwin.com)](https://cygwin.com).
   <br/>(The native "Command Prompt" interface could be used, but is not supported, i.e. use at your own risk)

2. Python.  
   * This pipeline is written in and supports <u>Python v2</u> versions 2.7.11 and later. (Python v3 can be used, though the pipeline has not been tested for changes between v2 and v3.)


  MiniConda (a lighter version offered by Anaconda) includes fewer packages but will still For a lighter version, MiniConda includes fewer packages but takes up less storage space.. (Anaconda also includes many standard packages used in pipeline.)

##### II.  Packages

  Note on <code>jupyter</code>: there are other <code>jupyter</code> components that are suffixed. These are installed and updated alongside <code>jupyter</code> as dependencies, so you don't have to worry about those.

  <b>How to install packages</b>

  <u><i>Using conda</i></u>

  Using your CLI (Terminal or Cygwin), type:

  > <code> conda install [your_package_name]</code> <br/>
  Example: <code> conda install my_new_package</code> <br/>

  The terminal will then print some information on what other dependencies are needed, and then ask:

  > <code>Proceed ([y]/n)?</code> (blinking cursor here)

  > (To make things explicit <in the zen of python>, type <code>n</code> (to cancel) or  <code>y</code> (to continue) and hit <code>Enter</code>.)

  <u><i>Using pip</i></u>

  If you don't have Anaconda, or in case conda can't find a package, it's possible to find and install it through <code>pip</code>. The syntax is similar to that for conda.

  > <code> pip install [your_package_name]</code> <br/>
  Example: <code> pip install conda_cant_find_me</code> <br/>

  If it's already installed, the terminal will print an error that says it's already installed and which version it is. In that case, you can choose to update it by typing:

  > <code> pip install [your_package_name] --upgrade </code> <br/>
  Example: <code> conda install im_already_here --upgrade</code> <br/>


  (To clarify, anaconda usually installs most of these packages while installing Python, but they are NOT a part of (i.e. native to) Python.) has many of these packages.) (versions in parentheses). (If already installed, check the version and update as necessary.)


### Starting up



_______________________
## Sample Data

Sample data have been provided in the [data/sample_data](https://github.com/sarahfortunelab/barcodetracking) folder. These data are generated from samples with known numbers (i.e. 1, 5, 25, and 120) of barcoded plasmids.

The data are presented are the raw compressed FASTQ (.fastq.gz) files generated in Illumina sequencing.

Results of an indexed sample in two files: one for forward read ("R1") and the second for reverse read ("R2") data. (Read more about Illumina's file nomenclature in [Illumina's CASAVA User Guide](http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm).

Below are samples included and the corresponding file pairs:

*  <b>1 barcode</b>:
   * NH001_S1_L001_R1_001.fastq.gz
   * NH001_S1_L001_R2_001.fastq.gz


*  <b>5 barcodes</b>:
  * NH005_S5_L001_R1_001.fastq.gz
  * NH005_S5_L001_R2_001.fastq.gz


*  <b>25 barcodes</b>:
   * NH025_S11_L001_R1_001.fastq.gz
   * NH025_S11_L001_R2_001.fastq.gz


*  <b>120 barcodes</b>:
   * NH120_S14_L001_R1_001.fastq
   * NH120_S14_L001_R2_001.fastq

<The following sample outputs are provided:

*  Filtered XXX
*  stats XXX
*  thresholded XXX

Note that the sample .db output file containing all parsed data is not included due to file size restrictions.>
