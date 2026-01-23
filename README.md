# Lock_2026_SI
Supporting information for Lock 2026: Tectonics on early Earth driven by Earth’s changing shape during tidal recession of the Moon

Hello, and thank you for accessing this repository. With it, you can reproduce the calculations, figures, and movies presented in the Lock 2026. The  Below are details of the contents of the repository. With any luck, you should be able to clone this repository and run it on your own computer. If you have any issues, please contact the author at s.lock@bristol.ac.uk

Requirements:
You will need an installation of python (including the packages given in the environment.yml file) and jupyter notebooks. For some figures you will also need to download an additional data repository from XXXX. By default this assumed to be one file level up from the scripts in this repo, but this can be changed in the parameters cell of the relevant scripts. 

Contents:

calculate_change_in_lenth_data_*.ipynb
Notebooks used to calculate data presented in Lock (2016) that is not authorwise contained in the othere plotting scipts. By default these scripts will make duplicates of some of the data in the Data directory or additional data repo (see above) in the relevant sub-directories of the Data directory with files indicated by *recalc*. 

Data
Directory containing python data binary files, hd5 files, and orbital evolution calculation output from previous works used in some plotting scripts. These files can be reproduced from the plotting scripts, or from the balculate* scripys, but are provided here for ease of reproducing the figures.

Earth_correct_params_S3.20c
This directory contains the output from the HERCULES planetary structure code used in many of the calculations of Lock 2026, which describe the structure of Earth at different angular momenta. Each sub-directory includes the input file for HERCULES and the associated output file which is in a custom binary format. A python structure for reading in such output is provided in the Supporting_scripts directory, but for more information the user is referred to the documentation of the HERCULES code at: https://github.com/sjl499/HERCULESv1_user

enivonment.yml
Conda environment file necessary for running the scripts in this repository. The scripts may run with different versions of the packages than those listed here, but there is no garauntee. You can load this environment file into conda by running XXXXXXX

Figures*.ipynb
Jupyter notebooks that reproduce each of the plots in the paper. In most cases, the script also calculates the data necessary to reproduce the plot. For more basic plots, this is done on the fly automatically in the code. In others, the code blocks that perform the calculations are not run by default, there being an optional flag to activate these blocks, and the default is instead to load data from the relevant file in the "Data" directory. For Figures 7-9 and S1-S4, the data is too large to be included in this git repository. The user must therefore first download the accompanying data repository and link to it as either 'Data_repo' in the base directory of this repo, or give the correct path in the plotting notebooks.

Figures*.pdf
pdf's of each of the figures in Lock 2026 as produced by the above notebooks. 

Figure_10
Directory containing the adobe illustrator and pdf for Figure 10.

Helvetica.ttc
Font file for Helvetica that is used in all plots.

LICENSE
The license under which you agree to abide if you use, adapt, etc. any of the information in this repository. Please read if you intend to use any of the resources provided here. 

Movies*.ipynb
Jupyter notebooks that produce the slides used to construct the supplemental movies. To produce the movie slides, the user must therefore first download the accompanying data repository and link to it as either 'Data_repo' in the base directory of this repo, or give the correct path in the plotting notebooks. Movie slides are produced and saved in the Movie_slides directory when the script is run. The slides are not provided by default due to space restrictions.

.Movie_slides
Note, this directory is not included automatically in this repository. This directory will be produced if you run any of the Movies*.ipynb notebooks.

README.md
Your reading it... spooky right?

Support_scripts
This directory contains a few python scripts (functions, colormaps etc.) that are needed to perform the calculations and make the plots. Most of these functions are detailed in the manual for the HERCULES code at: https://github.com/sjl499/HERCULESv1_user