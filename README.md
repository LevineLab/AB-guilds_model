# AB-guilds_model

The AB method is a statistical dimension reduction approach for identifying microbial metabolic guilds from binary data sets of functions (e.g., gene annotations). This method can be applied to any type of binary data and determines the underlying groups of features that are driving the distribution of feature presence/absence in the overall data set.

# Publication Reproduction

We provide a self-contained R markdown file and all of the requisite data objects necessary to reproduce the major figures and tables from the associated publication. We currently provide an early release that enables users to access the core AB functionality. However, this is a nascent project and in development, so if you are interested in applying this method to your data and have interest in specific functionalities and analyses that are not currently present, please feel free to reach out to github user ryanreyn directly to assist you.

# Current Release

Currently, we have a pre-release of the AB guilds method described by an interactive R script with defined arguments that a user can manipulate to run the core AB functionality on their own specific datasets. This release does not at the moment include full end-to-end support for annotation and plot analyses. This release has several R package dependencies that can be viewed within the code itself but are listed here as well: argparser, gtools, dplyr, magrittr, gplots, reshape2, and tidyverse.

# Citation

The publication associated with this work is Reynolds, R., Hyun, S., Tully, B., Bien, J. and Levine, N.M., Identification of Microbial Metabolic Functional Guilds from Large Genomic Datasets. Frontiers in Microbiology, 14, p.1197329.
