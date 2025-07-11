# Make sure community data is simulated first before it is cloned
ruleorder: simulate_community_data > clone_simulation_data

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


### Missing so far ###
rule all:
  input:
    expand("analysis_data/{experiment}/{treatment}/{simulation}/sign_ajd_res_log.csv", zip, 
           experiment=[r[0] for r in all_runs],
           treatment=[r[1] for r in all_runs],
           simulation=[r[2] for r in all_runs])


### Simulate community dynamics and output final timepoint community table 'n_pop_log_log.csv'
rule simulate_community_data:
  output:
    "simulation_data/{experiment}/{treatment}/{simulation}/n_pop_log_log.csv"
  params:
    is_clone=lambda wc: wc.experiment in clones and wc.treatment in clones[wc.experiment] and wc.simulation in clones[wc.experiment][wc.treatment]["simulations"],
    script=lambda wc: (
      experiments.get(wc.experiment, {}).get(wc.treatment, {})
        .get("simulations", {}).get(wc.simulation, {})
        .get("simulation_script", "None")
    ),
    sim_args=lambda wc: (
      experiments.get(wc.experiment, {}).get(wc.treatment, {})
        .get("simulations", {}).get(wc.simulation, {})
        .get("simulation_parameter_string", "")
    )
  run:
    if params.is_clone:
        print(f"Skipping simulation for clone: {wildcards.experiment}/{wildcards.treatment}/{wildcards.simulation}")
        return

    if params.script == "multi_simulation.py":
        print("Skipping multi_simulation.py — handled separately.")
        return

    shell(f"""
      python3 simulation_scripts/{params.script} \
        {wildcards.experiment} \
        {wildcards.treatment} \
        {wildcards.simulation} \
        {params.sim_args}
    """)


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


### Rule clones runs already conducted for alterations in down stream analysis. This is mainly the case when we add sampling noise
rule clone_simulation_data:
  output:
    "simulation_data/{experiment}/{treatment}/{simulation}/n_pop_log_log.csv"
  run:
    # Check if it's a clone
    try:
      source = clones[wildcards.experiment][wildcards.treatment]["simulations"][wildcards.simulation]
      source_treatment = source["source_treatment"]
      source_simulation = source["source_simulation"]
    except KeyError:
      raise ValueError(f"No source defined for cloned simulation {wildcards.experiment}/{wildcards.treatment}/{wildcards.simulation}")
    
    import shutil
    import os

    src_dir = f"simulation_data/{wildcards.experiment}/{source_treatment}/{source_simulation}"
    dst_dir = f"simulation_data/{wildcards.experiment}/{wildcards.treatment}/{wildcards.simulation}"
    # Ensure destination exists
    os.makedirs(os.path.dirname(dst_dir), exist_ok=True)
    # Copy tree (cloning)
    shutil.copytree(src_dir, dst_dir, dirs_exist_ok=True)


### Here I process the simlated community data: 
# Here we sample the habitats either individually or sampling at different volumes (depends on the analysis script that is called)
# We attain the community data for subsequent analysis and determining species co-occurrences through null-model testing (1000 permutations of random communities)
rule process_community_data:
  input:
    "simulation_data/{experiment}/{treatment}/{simulation}/n_pop_log_log.csv"
  output:
    "analysis_data/{experiment}/{treatment}/{simulation}/sign_ajd_res_log.csv"
  params:
    script=lambda wc: (
      experiments.get(wc.experiment, clones.get(wc.experiment, {}))
                .get(wc.treatment, {})
                .get("simulations", {})
                .get(wc.simulation, {})
                .get("analysis_script", "MISSING_SCRIPT.R")
    ),
    analysis_args=lambda wc: (
      experiments.get(wc.experiment, clones.get(wc.experiment, {}))
                .get(wc.treatment, {})
                .get("simulations", {})
                .get(wc.simulation, {})
                .get("analysis_parameter_string", "")
    )
  shell:
    """
    mkdir -p analysis_data/{wildcards.experiment}/{wildcards.treatment}/{wildcards.simulation}
    Rscript analysis_scripts/{params.script} \
      {wildcards.experiment} \
      {wildcards.treatment} \
      {wildcards.simulation} \
      "{params.analysis_args}"
    """





