configfile: "config.yaml"

# Extract experiments, treatments, and parameters from config
experiments = config["experiments"]
treatments = {exp: list(experiments[exp].keys()) for exp in experiments}

# Define final outputs for all combinations of experiments and treatments
#rule all:
#    input:
#        # Manually create a list of paths for each combination of exp and treatment
#        *["Result_master_dir/{exp}/done_collecting.txt".format(exp=exp, treatment=treatment)
#          for exp in experiments.keys()
#          for treatment in treatments[exp]]
rule all:
    input:
        ["Result_master_dir/{exp}/done_plotting.txt".format(exp=exp) for exp in experiments.keys()]

### I simulate Lotka-Volterra dynamics in a spatially explicit system.
# 'Outputs are a bunch of .csv's which are saved into my parent/subdirectory folder 'in-script'. 
# Currently I simulate the entire dataset before saving anything. In future I will append every run to result files.
# So I resort to generating a pseudo-output 'done.txt' which I update with every rule that is completed.
rule simulate_community_data:
    output:
        "Result_master_dir/{exp}/{treatment}/done_simulation.txt" # A dummy output file
    params:
        workdir_param=config["workdir"],
        param_string=lambda wildcards: experiments[wildcards.exp][wildcards.treatment]["param_string"],
        extra_params=lambda wildcards: " ".join(f"{k}={v}" for k, v in experiments[wildcards.exp][wildcards.treatment].get("extra_params", {}).items()) if experiments[wildcards.exp][wildcards.treatment].get("extra_params") else ""
    shell:
        """
        python3 simulation_code/master_slurm_simulation.py {params.workdir_param} {wildcards.exp} {wildcards.treatment} {params.param_string} {params.extra_params}
        touch {output}
        """

### Here I process the simlated community data in the follwing ways: 
# 1. The most important part of this step is determining species co-occurrences through null-model testing (1000 permutations of random communities).
# 2. I compare my inputs (interaction matrix and environmental preference similarity) with the calculated co-occurrences.
# 3. I calculate some other metrics (like diversity indices). A lot has been deleted from here
# 4. I create two files that contian all means of my n replicates.
# -> The first 'full_res.csv' contains basic metrics (co-occurrence network properties, diversity indices, etc.)
# -> The second 'infm_res.csv' contains the results from the input x co-occurrence comparisons
### This is currently a massive function and in future I would like to take it apart and run most components individually (especially the null-model part which could easilty be paralellised)
rule processing_community_data:
    input:
        "Result_master_dir/{exp}/{treatment}/done_simulation.txt"  # Ensure simulation is complete
    output:
        "Result_master_dir/{exp}/{treatment}/done_processing.txt" # A dummy output file
    params:
        workdir_param=config["workdir"],
        param_string=lambda wildcards: experiments[wildcards.exp][wildcards.treatment]["param_string"],
        extra_params=lambda wildcards: " ".join(f"{k}={v}" for k, v in experiments[wildcards.exp][wildcards.treatment].get("extra_params", {}).items()) if experiments[wildcards.exp][wildcards.treatment].get("extra_params") else ""
    shell:
        """
        Rscript simulation_code/master_slurm_analysis.R {params.workdir_param} {wildcards.exp} {wildcards.treatment} {params.extra_params}
        touch {output}
        """

### Here I take all results from my different treatments and plot the results of an entire 'experiment'
# Quite simple: Just takes all the results files of treatments in an experiment and creates one final .csv
# in the experiment folder that contains all results in one.
# It also adds some minor things like labels etc. for plotting.
rule collecting_treatment_data:
    input:
        lambda wildcards: [
            f"Result_master_dir/{wildcards.exp}/{treatment}/done_processing.txt"  # Ensure processing is complete
            for treatment in treatments[wildcards.exp]
        ]
    output:
        "Result_master_dir/{exp}/done_collecting.txt" # A dummy output file
    params:
        workdir_param=config["workdir"],
    shell:
        """
        Rscript simulation_code/collect_treatments.R {params.workdir_param} {wildcards.exp} 
        touch {output}
        """

### Here I take all results from my different treatments and plot the results of an entire 'experiment'
# Quite simple: Just takes all the results files of treatments in an experiment and creates one final .csv
# in the experiment folder that contains all results in one.
# It also adds some minor things like labels etc. for plotting.
rule plot_fig:
    input:
        lambda wildcards: [
            f"Result_master_dir/{wildcards.exp}/done_collecting.txt"  # Ensure processing is complete
        ]
    output:
        "Result_master_dir/{exp}/done_plotting.txt" # A dummy output file
    params:
        workdir_param=config["workdir"],
    shell:
        """
        Rscript simulation_code/plotFig.R {params.workdir_param} {wildcards.exp} 
        touch {output}
        """
