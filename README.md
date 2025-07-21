# Simulating spatially explicit community data and analyising species co-occurrince

**A Snakemake workflow to simulate community data with a spatially and environmentally explicit generalized Lotka–Volterra model, to analyse the inference potential of species co‑occurrences, to reconstruct simulated interactions and environmental preference similarity**

---

## Table of Contents

1. [Background](#background)
2. [Code featured](#features)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Configuration](#configuration)
6. [A brief introduction: What is Snakemake?](#snakemake)
7. [Workflow overview](#workflow-overview)
8. [Plotting](#plotting)
9. [config.yaml settings](#config-settings)
10. [Repository structure](#repo-structure)
11. [Citation](#citation)
12. [License & Contact](#license--contact)

---

## Background

This workflow allows:

- The simulation of community data:
  - With any pairwise interaction coefficients
  - With a novel fundamental niche dynamic in which species respond to multiple environmental factors and can have varying degrees of specialisation
  - In a spatially explicit 3D network of habitats:
    - With dispersal dynamics constrained by habitat connectivity
    - With variable spatial distributions of environmental factors determining the extent of spatial autocorrelation of environmental drivers
    - With the ability to sample either individual habitats or aggregate community data in composite measurements at different spatial scales

- The analysis of simulated community data:
  - To determine richness and communities Bray-Curtis dissimilarities
  - To observe positive and negative co-occurrence via significant correlations and coefficient thresholds
  - To match co-occurrences with initial interaction and environmental preference inputs and determine the inference potential of co-occurrence data under given settings
  - To trace how the inference of either driver is influenced by the spatial scale of a sampling unit

## Code featured

- **Python** scripts for simulating community data
- **R** scripts for data processing, analysis, compiling data and plotting
- **Snakemake** workflow for reproducible, dependency‑aware execution with dynamic computational resource utilisation

## Requirements

- **Conda** (run with conda 25.5.1)
  - For python-, R-, module- and package-versions see conda.yaml 
- **Snakemake** (run with 7.32.4)

## Installation

1. **Clone .git repository**
(this will produce a folder in your terminals current directory)
   ```bash
   git clone https://github.com/rmSeanDarcy/SAMSARA.git
   cd SAMSARA
   ```
2. **Create and activate Conda environment**
   ```bash
   conda env create -f SAMSARA_env.yaml
   conda activate SAMSARA_env
   ```

## Configuration

We provide three config files:
- `config_test.yaml` -> Demo simulates scenarios covering all snakemake rules (few replicates; short execution <5min) 
- `config_fig2.yaml` -> Simulates all data required for producing Fig.2 (intermediate run time)
- `config_full.yaml` -> Simulates all data utilised in publication (very long runtime, execution on computing cluster recommended)

**Instructions:**

A. Copy, rename or symlink your chosen config as `config.yaml`:
   ```bash
   cp config-full.yaml config.yaml  # or config-test.yaml
   ```
B. Create your own `config.yaml` and set parameters:
   - `sp_num`: Number of species in simulation
   - `ratio_cf`: Ratio of negative to positive interactions
   - `radius`: Physical distance between nodes within which edges are established
   - **For list and all explanations of all parameters open simulation_scripts/simulation_functions/set_standard_initial_conditions.py. This also allows you to create a fixed set of parameters yourself.**

(You can find specific meaning of inputs for config files below)


## A brief introduction: What is Snakemake?
Snakemake is a workflow management system that helps you define and execute a series of data processing steps in a reproducible and efficient way. You write a Snakefile describing rules, each of which specifies input files, the code or shell commands to run, and output files. Snakemake then automatically builds a directed acyclic graph (DAG) of all steps, determines which ones need to run based on file timestamps, and schedules paralell execution leveraging the number of cores as specified. It can run either locally or on a computing cluster. It ensures that:
- Dependencies are handled automatically: Steps are only rerun when inputs change
- Parallelization is automatically applied: use the --cores flag to run independent rules simultaneously
- Reproducibility is enforced: the entire workflow is version-controlled and self-documenting (see hidden .snakemake folder for automatic documentation)
This makes it ideal for complex pipelines like our simulation and analysis where many intermediate files and language environments are orchestrated.


## Workflow overview

Execute the snakemake workflow (within your active conda environment):
  ```bash
  snakemake --snakefile Snakefile --configfile config.yaml --cores 6
  ```

Snakemake sequentially executes the following steps:
1. **Simulation**
Generation of community data in our spatially explicit generalised Lotka-Volterra system 
- `simulation_scripts/main_simulation.py` Main simulation script
- `simulation_scripts/multi_simulation.py` Simulates community data for Fig.5 in which interactions and species resource preferences are identical in the base (s1) and two control scenarios (c1, c2)
2. **Processing** 
Subsamples habitats/samples (N = 25) and calculates significant species co-occurrences and species-resource associations
- `analysis_scripts/process_habitat.R` N habitat samples ->  Single output
- `analysis_scripts/process_composite_fixed.R` Samples N composite samples containing fixed number of habitats -> Multiple outputs (for every composite sample side lengths l) 
- `analysis_scripts/process_composite_unequal.R` Samples N composite samples allowing variation in number of habitats -> Multiple outputs (for every composite sample side lengths l)
3. **Analysis**
Matches co-occurrence matrices with interaction matrices and species environmental preference similarity. Analyses species richness and Bray-Curtis dissimilarity
-`analysis_scripts/analysis.R`
4. **Collection**
Compiles data for every Figure into two dataframes (one for matching data, one for community metrics)
-`analysis_scripts/compile_data.R` Output data is now found in analysis_data/"figure ID"

Documentation on what individual steps do can be found in the scripts listed above. Each script referenced above calls functions from corresponding folders within either `simulation_scripts/' or 'analysis_scripts/'. Open the snakefile (in a text editor) to find the documented workflow code. For every simulation defined in our config data, snakemake runs the simluation (python), processing and analysis (R) of community data sequentially. Each simulation is fully independent, meaning all simulations can be run in parallel. Only in the final step, where data is compiled, snakemake has to wait until all data within a figure ID is finished before it can create the final dataset. 


## Plotting

All compiled result data for all simulations within a figure ID is stored under 'analysis_data/'. Plotting result data is not integrated into our workflow. If the full dataset (config_full.yaml) has been generated users can execute the following bash command to generate all .svg's used in the result figures. Depending on your setup, the working directory might need to be set to the working directory in the script.
```bash
Rscript analysis_scripts/plot_figures.R
```
Otherwise, or if only a part of the data was generated (f.ex. when running config_fig2.yaml), we recommend a line by line execution in an IDE like Rstudio. All functions called in `analysis_scripts/plot_figures.R` can be found in the corresponding scripts in 'analysis_scripts/plot_functions'.


## config.yaml settings

Every simulation is located in a hierarchy of folders with a unique path within 'simulation_data/'. Each is characterised by an 'experiment', a 'treatement' and a 'simulation' ID, which is also give folder names. 'simulation_data/fig2/s1/s1' f.ex. contains all simulated and result data from out base scenario. Each simulation listed under experiments receives four inputs:
- simulation_parameter_string: must take as a first element a set of standard initial conditions (see `simulation_scripts/simulation_functions/set_standard_initial_conditions.py` for more information). Additional parameters can be passed after this that override the standard initial conditions (separate by spaces)
- simulation_script: Mostly `simulation_scripts/main_simulation.py`. `simulation_scripts/multi_simulation.py` only called once for fig5
- analysis_parameter_string: Currently only argument for setting extent of noise (xi) implemented
- analysis_script: Specifies script for sampling either habitat level or composites (and here fixed or unequal number of habitats) 

Besides the first section under the header 'experiments:', there is a second section containing 'clones:'. As multiple downstream analyses can use the same simulated data (think analysis of sampling noise), we add a step in which this data is copied (or cloned) into new 'simulation' folders. For these simulations, the three level hierarchy is identical and all cloned and analysed data will appear in folders within the 'experiments/' folders also defined in the first config section. To specify which data is copied into these new paths 'simulations' under clones receive two additional inputs:
- source_treatment: Specifies the ID of the treatment of the data to be copied
- source_simulation: Specifies the ID of the simulation of the data to be copied
As with the simulations in the first config section, the path of a cloned file will follow the IDs given in the config file


## Repository structure

# Folder structure in .git repository
```text
├── README.md                         # You are here
├── analysis_scripts/                 # R processing, analysis, compiling and plotting code
│   ├── analysis_functions/         
│   │   ── ...
│   ├── analysis.R
│   ├── plot_figures.R                # Run manually for plotting
│   └── ...
├── simulation_scripts/               # Python simulation code
│   ├── simulation_functions/         
│   │   ── ...
│   ├── main_simulation.py
│   └── main_simulation.py
├── config_test.yaml                  # Demo settings
├── config_fig2.yaml                  # Allows recreation of Fig. 2
├── config_full.yaml                  # Simulates all data
├── SAMSARA_env.yaml                  # File from which conda environment is installed 
├── Snakefile                         # Workflow code
└── .gitignore                        # excludes .Rproj, .snakemake/, analysis_data/, simulation_data/ 
```

# Additional folders/files added in workflow (examples)
```text
├── simulation_data                   # Location where all simulated data is stored 
│   ├── fig2/                         # 'experiment ID'
│       └── s1/                       # 'treatment ID'
│           └── s1/                   # 'simulation ID'
│   └── ...
├── analysis_data                     # Contains data sets compiled on the 'experiment ID' level
│   ├── fig2/
│       ├── infm_res_compiled.csv     # Called in `plot_figures.R`
│   └── ...
└── figures                           # Output location for plots
```


## Outputs & Interpretation

All simulations are compiled for the various sections in the manuscript. Two main result files which are used for plotting:
- `analysis_data/fig2/infm_res_compiled.csv`: Data on the matching of co-occurrences with initial drivers 
- `analysis_data/fig2/full_res_compiled.csv`: Typical community data metrics (diversity indices)
A list of what column names represent can be found in `analysis_data/readme.csv`


## Citation

If you use this workflow, please cite:
> **[MISSING INFO: Your Name et al. (2025). Title. Journal. DOI]**


## License & Contact

This project is licensed under the MIT License.\
Questions or issues? Please open an issue or contact **[MISSING INFO: **[**your.email@institution.edu**](mailto\:your.email@institution.edu)**]**.

