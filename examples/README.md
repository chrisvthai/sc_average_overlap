# sc_average_overlap - Example code
We have published a manuscript of this work, which can be found [here](https://www.biorxiv.org/content/10.1101/2025.05.06.652497v1).

In this folder, we provide 3 main pieces of example code using average overlap.

The first two pertain to the main figures presented in the manuscript. First, an example script containing code to produce the data for average overlap benchmarking against other correlation-based methods in hierarchical clustering of groups in scRNA-seq analysis is provided, as well as the job script for running on a compute cluster. These are shown for demonstration purposes, and pertains to the first half of the paper, or Figure 1.

Second, a Python notebook that reproduces all figures for the thymocyte analysis is given. Code and instructions are provided for running Piccolo normalization and feature selection within an R environment using `rpy2`, however, the bulk of analysis implemented in this example notebook uses Piccolo counts saved to an external CSV file. This pertains to Figure 2 in our paper.

Lastly, an example tutorial that goes through a simple workflow for scRNA-seq analysis of PBMCS incorporating average overlap can be found in `3kpbmcs_tutorial.ipynb`. 