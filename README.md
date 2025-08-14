Code and data for the submitted manuscript "Reconstructing enterovirus D68 clade-specific transmission dynamics using phylodynamics and association with acute flaccid myelitis."

## Downloading the required sequence data
The included data and code allow a user to recreate the analyses performed in the manuscript. The sequence data are not included here but are publicly available on genbank with the accession numbers provided in the [metadata](./data/metafile.csv). Code to perform a genbank download in R and prepare a metadata file is provided in the [code](./code/1_genbank_pull.R) folder.

## Subsampling the data
In this analysis we perform multiple samples of the data under different schemes to test the sensitivity of the model to the included sequences and the sampling scheme used. The code used to generate these samples is provided in [2_subsampling.R](./code/2_subsampling.R) and code to plot histograms of the samples is provided in [3_plot_samples.R](./code/3_plot_samples.R).

## Analyses with BEAST2
Samples generated were used to run a phylodynamic model with [BDMM-Prime](https://github.com/tgvaughan/BDMM-Prime) in BEAST2. The BEAST code required for these analyses is provided in [BEAST_scripts](./code/BEAST_scripts). These files use the [feast](https://github.com/tgvaughan/feast) package to load sequence data and parse dates, allowing them to be run on any of the pairs of alignments generated. Feast needs to be installed via the beauti package manager or command line to run the scripts. The xml files can be added to any of the sample folders generated in the previous section, a change will need to be made to the file path for the sequence data to match the path on your system. The included scripts log the population trajectory every 250000 iterations. This can be commented out to generate only the Rt trajectories which is faster and less memory intensive to run. 

To test sensitivity to infectious period we ran the model three times on the same sample, changing only the time to recovery. The scripts for this sensitivity analysis can be found in [sensitivity_infectious_period](./code/BEAST_scripts/sensitivity_infectious_period)

## Generating plots
The Rt trajectories presented in Figure 3 can be generated with the code in [4_rt_trajectory.R](./code/4_rt_trajectory.R). This plot also requires the metadata file for the full dataset which can be found in [data](./data/metafile.csv). Plotting the population trajectories and testing the association with AFM cases is shown in [5_population_trajectory.R](./code/5_population_trajectory.R).
