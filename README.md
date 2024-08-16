# AB-guilds_model

The AB method is a statistical dimension reduction approach for identifying microbial metabolic guilds from binary data sets of functions (e.g., gene annotations). This method can be applied to any type of binary data and determines the underlying groups of features that are driving the distribution of feature presence/absence in the overall data set.

# Current Release

We currently provide an early release that enables users to access the core AB functionality. However, this is a nascent project and in development, so if you are interested in applying this method to your data and have interest in specific functionalities and analyses that are not currently present, please feel free to reach out to github user ryanreyn directly to assist you. Release v0.1.0 contains an interactive R script with defined arguments that a user can manipulate to run the core AB functionality on their own specific datasets. This script, along with helper functions are provided in the ***Codes*** folder. This release does not at the moment include full end-to-end support for annotation and plot analyses. This release has several R package dependencies that can be viewed within the code itself but are listed here as well: argparser, gtools, dplyr, magrittr, gplots, reshape2, and tidyverse.

# Running the AB script

The AB script can be run by using the ***Rscript*** terminal function to call the main script, ***Run_AspectBernoulli.R***. To see all of the arguments that can be passed to the main script, navigate to the ***Codes*** folder and run the following line of code:

```
Rscript Run_AspectBernoulli.R --help
```

Running the main script will generate 3 output files in the defined output folder representing the full set of functions ordered by the specificity score of each function in each guild, the top scoring functions in each guild (number of functions determined by the user), and the number of mapback genomes per guild (based on the top X functions).

# Publication Reproduction

We provide a self-contained R markdown file and all of the requisite data objects necessary to reproduce the major figures and tables from the associated publication. These files can be found in the ***reproduce-publication_release*** folder.

# Citation

The publication associated with this work is Reynolds, R., Hyun, S., Tully, B., Bien, J. and Levine, N.M., Identification of Microbial Metabolic Functional Guilds from Large Genomic Datasets. Frontiers in Microbiology, 14, p.1197329.
