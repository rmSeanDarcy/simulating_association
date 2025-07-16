# Make sure community data is simulated first before it is cloned
ruleorder: clone_simulation_data > generate_simulation_data

import yaml

### Load config file
with open("config_test.yaml") as f:
    config = yaml.safe_load(f)
experiments = config["experiments"]
clones = config.get("clones", {})

### List of all runs = (experiment, treatment, simulation)
all_runs = []
for fig, treatments in experiments.items():
    for treatment, treatment_info in treatments.items():
        for sim in treatment_info["simulations"]:
            all_runs.append((fig, treatment, sim))
for fig, treatments in clones.items():
    for treatment, treatment_info in treatments.items():
        for sim in treatment_info["simulations"]:
            all_runs.append((fig, treatment, sim))


rule all:
    input:
        expand(
            "analysis_data/{experiment}/infm_res_compiled.csv",
            experiment=config["experiments"].keys()
        ),
        expand(
            "analysis_data/{experiment}/full_res_compiled.csv",
            experiment=config["experiments"].keys()
        ),
        # Also add clones outputs to analysis if you want them compiled
        # Or at least ensure simulation data exists for clones so analysis can run


### Missing so far ###
#rule all:
#  input:
#    expand("simulation_data/{experiment}/{treatment}/{simulation}/infm_res.csv", zip, 
#           experiment=[r[0] for r in all_runs],
#           treatment=[r[1] for r in all_runs],
#           simulation=[r[2] for r in all_runs])


### Simulate community dynamics and output final timepoint community table 'n_pop_log_log.csv'
rule generate_simulation_data:
    output:
        "simulation_data/{experiment}/{treatment}/{simulation}/n_pop_log_log.csv"
    run:
        is_clone = (
            wildcards.experiment in clones
            and wildcards.treatment in clones[wildcards.experiment]
            and wildcards.simulation in clones[wildcards.experiment][wildcards.treatment]["simulations"]
        )
        if is_clone:
            print(f"[generate_simulation_data] Skipping simulation for clone: {wildcards.experiment}/{wildcards.treatment}/{wildcards.simulation}")
            # do nothing here, cloning handled by separate rule
            return

        # Non-clone simulation
        script = experiments.get(wildcards.experiment, {}) \
                            .get(wildcards.treatment, {}) \
                            .get("simulations", {}) \
                            .get(wildcards.simulation, {}) \
                            .get("simulation_script", None)
        sim_args = experiments.get(wildcards.experiment, {}) \
                              .get(wildcards.treatment, {}) \
                              .get("simulations", {}) \
                              .get(wildcards.simulation, {}) \
                              .get("simulation_parameter_string", "")

        if script == "multi_simulation.py":
            print("[generate_simulation_data] Skipping multi_simulation.py — handled separately.")
            return

        shell(f"""
            echo "[generate_simulation_data] Running simulation for: {wildcards.experiment}/{wildcards.treatment}/{wildcards.simulation}"
            python3 simulation_scripts/{script} {wildcards.experiment} {wildcards.treatment} {wildcards.simulation} {sim_args}
        """)

### Clone simulated community data where necessary
rule clone_simulation_data:
    input:
        src_npop = lambda wc: f"simulation_data/{wc.experiment}/" +
                             f"{clones[wc.experiment][wc.treatment]['simulations'][wc.simulation]['source_treatment']}/" +
                             f"{clones[wc.experiment][wc.treatment]['simulations'][wc.simulation]['source_simulation']}/n_pop_log_log.csv"
    output:
        "simulation_data/{experiment}/{treatment}/{simulation}/n_pop_log_log.csv"
    run:
        import os
        import shutil

        # source and destination directories
        src_treatment = clones[wildcards.experiment][wildcards.treatment]["simulations"][wildcards.simulation]["source_treatment"]
        src_sim = clones[wildcards.experiment][wildcards.treatment]["simulations"][wildcards.simulation]["source_simulation"]
        src_dir = f"simulation_data/{wildcards.experiment}/{src_treatment}/{src_sim}"
        dst_dir = f"simulation_data/{wildcards.experiment}/{wildcards.treatment}/{wildcards.simulation}"

        if not os.path.exists(src_dir):
            raise ValueError(f"[clone_simulation_data] Source directory does not exist: {src_dir}")

        os.makedirs(dst_dir, exist_ok=True)

        skip_files = {"sign_ajd_res_log.csv"}  # processed files to skip

        copied_files = 0
        for filename in os.listdir(src_dir):
            if filename in skip_files:
                print(f"[clone_simulation_data] Skipping processed file: {filename}")
                continue
            src_path = os.path.join(src_dir, filename)
            dst_path = os.path.join(dst_dir, filename)
            if os.path.isfile(src_path):
                shutil.copy2(src_path, dst_path)
                copied_files += 1
                print(f"[clone_simulation_data] Copied file: {filename}")

        if not os.path.exists(os.path.join(dst_dir, "n_pop_log_log.csv")):
            raise ValueError(f"[clone_simulation_data] Clone succeeded but output missing: {os.path.join(dst_dir, 'n_pop_log_log.csv')}")

        print(f"[clone_simulation_data] Total files copied: {copied_files}")


### Rule runs the script 'multi_simulation.py', which takes no inputs and is 'hard coded' to produce the data for fig5 
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!change the xfig5 to fig5!
rule simulate_fig5_multi:
  output:
    "simulation_data/xfig5/nw_s1/nw_s1/n_pop_log_log.csv",
    "simulation_data/xfig5/nw_c1/nw_c1/n_pop_log_log.csv",
    "simulation_data/xfig5/nw_c2/nw_c2/n_pop_log_log.csv"
  shell:
    """
    echo "Running multi_simulation.py for fig5"
    python simulation_scripts/multi_simulation.py
    """


### Here I process the simlated community data: 
# Here we sample the habitats either individually or sampling at different volumes (depends on the analysis script that is called)
# We attain the community data for subsequent analysis and determining species co-occurrences through null-model testing (1000 permutations of random communities)
rule process_community_data:
    input:
        "simulation_data/{experiment}/{treatment}/{simulation}/n_pop_log_log.csv"
    output:
        "simulation_data/{experiment}/{treatment}/{simulation}/sign_ajd_res_log.csv"
    params:
        script=lambda wc: (
            clones.get(wc.experiment, {}).get(wc.treatment, {})
                  .get("simulations", {})
                  .get(wc.simulation, {})
                  .get("analysis_script")
            or
            experiments.get(wc.experiment, {}).get(wc.treatment, {})
                       .get("simulations", {})
                       .get(wc.simulation, {})
                       .get("analysis_script", "MISSING_SCRIPT.R")
        ),
        analysis_args=lambda wc: (
            clones.get(wc.experiment, {}).get(wc.treatment, {})
                  .get("simulations", {})
                  .get(wc.simulation, {})
                  .get("analysis_parameter_string")
            or
            experiments.get(wc.experiment, {}).get(wc.treatment, {})
                       .get("simulations", {})
                       .get(wc.simulation, {})
                       .get("analysis_parameter_string", "")
        )
    shell:
        """
        Rscript analysis_scripts/{params.script} {wildcards.experiment} {wildcards.treatment} {wildcards.simulation} "{params.analysis_args}"
        """


### Next I analyse the processed community data:
# Significant co-occurrences are tested for their concordance with initial settings 
# (i.e. interaction matrix and environmental preference similarity)
# We also extract some simple community metrics such as richness and Bray-Curtis dissimilarity
rule analyse_simulations:
    input:
        "simulation_data/{experiment}/{treatment}/{simulation}/sign_ajd_res_log.csv"
    output:
        "simulation_data/{experiment}/{treatment}/{simulation}/infm_res.csv",
        "simulation_data/{experiment}/{treatment}/{simulation}/full_res.csv"
    params:
        script="analysis_scripts/analysis.R"
    shell:
        """
        Rscript {params.script} {wildcards.experiment} {wildcards.treatment} {wildcards.simulation}
        """


### Here we compile our simulated and analysed data into two dataframes per experiment
# The first contains all of the matching results (columns) for each simulation (row)
# The second gives us the community data metrics in the same format
# We use these dataframes for plotting
# Helper function to collect all sim output files for one experiment
import os

# Extract all experiments from the config
EXPERIMENTS = list(config["experiments"].keys())

# Collect all infm_res.csv and full_res.csv paths from simulation_data
def get_simulation_results_per_experiment(experiment):
    treatments = config["experiments"][experiment]
    files = []
    for treatment in treatments:
        simulations = treatments[treatment]["simulations"]
        for simulation in simulations:
            base = f"simulation_data/{experiment}/{treatment}/{simulation}"
            files.append(f"{base}/infm_res.csv")
            files.append(f"{base}/full_res.csv")
    return files

rule compile_data:
    input:
        lambda wildcards: get_simulation_results_per_experiment(wildcards.experiment)
    output:
        "analysis_data/{experiment}/infm_res_compiled.csv",
        "analysis_data/{experiment}/full_res_compiled.csv"
    shell:
        """
        Rscript analysis_scripts/compile_data.R {wildcards.experiment}
        """
