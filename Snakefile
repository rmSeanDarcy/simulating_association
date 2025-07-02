configfile: "config.yaml"

# Extract experiments, treatments, and parameters from config
experiments = config["experiments"]
treatments = {exp: list(experiments[exp].keys()) for exp in experiments}

# Define final outputs for the plotting rule
rule all:
    input:
        expand("Result_master_dir/{experiment}/report_{experiment}.pdf", experiment=config["experiments"].keys())  # Add LaTeX output here

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

### Here I collect all results from my different treatments of an entire 'experiment'
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
        workdir_param=config["workdir"]
    shell:
        """
        Rscript simulation_code/collect_treatments.R {params.workdir_param} {wildcards.exp} 
        touch {output}
        """

### Here I take all the collected results and plot the results of an entire 'experiment'
# This is currently just a barplot with each treatment as a bar. A total bar represents all co-occurrences,
# Different parts of the bar are coloured by how many co-occurrences corresponded to interactions and/or
# environmental preference similarity
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

### Now I will load my generated figures into a simple Rmarkdown document, which in a more extensive
# version could contain a whole range of analyses compiled into a single .pdf, which the user could
# inspect after running the code
rule rmarkdown:
    input:
        "Result_master_dir/{exp}/done_plotting.txt",  # Ensure the plotting is done
        #"Result_master_dir/{exp}/Fig2_B.png",  # The generated figure
        #"Result_master_dir/{exp}/Fig2_C.png",  # The generated figure
    output:
        "Result_master_dir/{exp}/report_{exp}.pdf",  # The generated PDF report
    params:
        rmd_file="report_template.Rmd",  # Path to the R Markdown template
    shell:
        """
        Rscript -e "rmarkdown::render('report_template.Rmd', output_file='Result_master_dir/{wildcards.exp}/report_{wildcards.exp}.pdf', params = list(experiment = '{wildcards.exp}'))"
        """
