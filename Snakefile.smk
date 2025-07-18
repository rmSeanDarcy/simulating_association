############################################################################
##### Snakefile for simulating and analysing community data            #####
############################################################################

############################################################################
### Load config files and modules 
configfile: "config_test.yaml"
import os, shutil

### Define labels for data that needs to be:
# 1. To be simulated (experiments in config.yaml)
# 2. To be cloned from simulated data (clones in config.yaml). This reduces
#    the need to simulate data with the same settings multiple times
experiments = config["experiments"]
clones      = config["clones"]

############################################################################
##### Enforce rule order (important)
ruleorder: simulate > clone > process_sim > process_clone > analyse_simulations > compile_data


############################################################################
##### 1. Define outputs and 'rule all'                                 ##### 
############################################################################

############################################################################
##### Define outputs for rule all
# Due to some conflicts with downstream analysis we resort to defining all 
# intermediate outputs as outputs for the total workflow.
### Output: rule simulate
SIM_MARKERS = [
    f"simulation_data/{exp}/{tr}/{sim}/.simulated"
    for exp, tr_d in experiments.items()
    for tr, sim_d in tr_d.items()
    for sim in sim_d["simulations"]
]
### Output: rule clone
CLONE_MARKERS = [
    f"simulation_data/{exp}/{tr}/{sim}/.cloned"
    for exp, cl_d in clones.items()
    for tr, sim_d in cl_d.items()
    for sim in sim_d["simulations"]
]
### Output: rule process
PROC_TARGETS = []
# from each simulate‐marker
for m in SIM_MARKERS:
    PROC_TARGETS.append(m.replace("/.simulated", "/sign_ajd_res_log.csv"))
# from each clone‐marker
for m in CLONE_MARKERS:
    PROC_TARGETS.append(m.replace("/.cloned", "/sign_ajd_res_log.csv"))
### Output: rule analyse
ANALYSIS_TARGETS = []
for m in SIM_MARKERS + CLONE_MARKERS:
    prefix = m.rsplit("/",1)[0]
    ANALYSIS_TARGETS += [
        f"{prefix}/infm_res.csv",
        f"{prefix}/full_res.csv",
    ]
### Output: rule compile
COMPILE_TARGETS = []
for exp in experiments:
    d = f"analysis_data/{exp}"
    COMPILE_TARGETS += [f"{d}/infm_res_compiled.csv", f"{d}/full_res_compiled.csv"]

############################################################################
##### Define rule all
rule all:
    input:
        SIM_MARKERS,
        CLONE_MARKERS,
        PROC_TARGETS,
        ANALYSIS_TARGETS,
        COMPILE_TARGETS


############################################################################
##### 2. Simulate community data 
############################################################################

############################################################################
##### Main rule that simulates community data
# In this rule we simulate assembly and produce the community data
# We strictly separate this from the next cloning rule. We produce a dummy output
# '.simulated' in stead of actual data to differentiate these data from cloned data
rule simulate:
    input:
        script=lambda wc: (
            "simulation_scripts/"
            + experiments[wc.exp][wc.tr]["simulations"][wc.sim]["simulation_script"]
        )
    output:
        marker="simulation_data/{exp}/{tr}/{sim}/.simulated"
    params:
        simargs=lambda wc: experiments[wc.exp][wc.tr]["simulations"][wc.sim]["simulation_parameter_string"]
    run:
        out = f"simulation_data/{wildcards.exp}/{wildcards.tr}/{wildcards.sim}"
        os.makedirs(out, exist_ok=True)
        
        ### This includes an exception for when 'multi_simulation.py' is called (fig5)  
        script_file = os.path.basename(input.script)
        if script_file == "multi_simulation.py":
            print(f"[simulate] SKIP {script_file} — handled by simulate_fig5_multi.")
            return

        if not os.path.exists(output.marker):
            shell(f"python3 {input.script} {wildcards.exp} {wildcards.tr} {wildcards.sim} {params.simargs}")
            open(output.marker, "w").close()
        else:
            print(f"[simulate] SKIP {wildcards.exp}/{wildcards.tr}/{wildcards.sim}")

############################################################################
##### Rule runs the script 'multi_simulation.py' (no inputs, 'hard coded' for fig5) 
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!change the xfig5 to fig5!
rule simulate_fig5_multi:
    output:
        "simulation_data/xfig5/nw_s1/nw_s1/.simulated",
        "simulation_data/xfig5/nw_c1/nw_c1/.simulated",
        "simulation_data/xfig5/nw_c2/nw_c2/.simulated"
    run:
        # run your multi‐simulation
        shell('''
            echo "Running multi_simulation.py for fig5"
            python simulation_scripts/multi_simulation.py
        ''')
        # now explicitly create each marker
        for marker in output:
            dirpath = os.path.dirname(marker)
            os.makedirs(dirpath, exist_ok=True)
            open(marker, "w").close()


############################################################################
##### 3. Clone community data for runs with repeat analyses
############################################################################

############################################################################
##### Rule clones community data
# config.yaml contains a section of under 'clones'. Here we define analyses
# that can build on already generated data. F.ex. the analysis of noise can be
# conducted on data simulated in basic scenarios. Again we give a dummy output
# '.cloned' due to many downstream snakemake dependency issues 
rule clone:
    input:
        src_marker=lambda wc: (
            f"simulation_data/{wc.exp}/"
            + clones[wc.exp][wc.tr]["simulations"][wc.sim]["source_treatment"]
            + f"/{clones[wc.exp][wc.tr]['simulations'][wc.sim]['source_simulation']}/.simulated"
        )
    output:
        marker="simulation_data/{exp}/{tr}/{sim}/.cloned"
    run:
        exp, tr, sim = wildcards.exp, wildcards.tr, wildcards.sim
        src = os.path.dirname(input.src_marker)
        dst = f"simulation_data/{exp}/{tr}/{sim}"
        os.makedirs(dst, exist_ok=True)
        if not os.path.exists(output.marker):
            print(f"[clone] RUN  {exp}/{tr}/{sim} ← copying from {src}")
            for fn in os.listdir(src):
                if fn in (".simulated", ".cloned", "sign_ajd_res_log.csv"):
                    continue
                shutil.copy2(os.path.join(src, fn), os.path.join(dst, fn))
            open(output.marker, "w").close()
        else:
            print(f"[clone] SKIP {exp}/{tr}/{sim}")


############################################################################
##### 4. Process community data
############################################################################
# 'Processing' data means that here we:
# 1. Sample either a fixed number of habitats or a certain number of cubes in
# different levels of volume
# 2. Calculate significant species co-occurrences / species x resource associations
#
# This implementation is not pretty but it works. Due to dependency issues we 
# split the processing of simulated and cloned community data into two rules

############################################################################
##### Rule processes simulated community data
rule process_sim:
    input:
        ### Only run for runs in 'experiments'
        _check=lambda wc: (
            []
            if wc.exp in experiments
               and wc.tr in experiments[wc.exp]
               and wc.sim in experiments[wc.exp][wc.tr]["simulations"]
            else
               (_ for _ in ()).throw(KeyError(f"Not an experiment: {wc}"))
        ),
        marker="simulation_data/{exp}/{tr}/{sim}/.simulated"
    output:
        "simulation_data/{exp}/{tr}/{sim}/sign_ajd_res_log.csv"
    params:
        script=lambda wc: experiments[wc.exp][wc.tr]["simulations"][wc.sim]["analysis_script"],
        args=lambda wc: experiments[wc.exp][wc.tr]["simulations"][wc.sim]["analysis_parameter_string"]
    shell:
        "Rscript analysis_scripts/{params.script} "
        "{wildcards.exp} {wildcards.tr} {wildcards.sim} '{params.args}'"

############################################################################
##### Rule processes cloned community data
rule process_clone:
    input:
        ### Only run for runs in 'clones'
        _check=lambda wc: (
            []
            if wc.exp in clones
               and wc.tr in clones[wc.exp]
               and wc.sim in clones[wc.exp][wc.tr]["simulations"]
            else
               (_ for _ in ()).throw(KeyError(f"Not a clone: {wc}"))
        ),
        marker="simulation_data/{exp}/{tr}/{sim}/.cloned"
    output:
        "simulation_data/{exp}/{tr}/{sim}/sign_ajd_res_log.csv"
    params:
        script=lambda wc: clones.get(wc.exp,{}).get(wc.tr,{}).get("simulations",{}).get(wc.sim,{}).get("analysis_script")
                             or experiments[wc.exp][wc.tr]["simulations"][wc.sim]["analysis_script"],
        args=lambda wc:   clones.get(wc.exp,{}).get(wc.tr,{}).get("simulations",{}).get(wc.sim,{}).get("analysis_parameter_string","")
                             or experiments[wc.exp][wc.tr]["simulations"][wc.sim]["analysis_parameter_string"]
    shell:
        "Rscript analysis_scripts/{params.script} "
        "{wildcards.exp} {wildcards.tr} {wildcards.sim} '{params.args}'"


############################################################################
##### 5. Analyse co-occurrence and community data
############################################################################

############################################################################
##### Main analysis function
# This funciton is the same for every run and gives us two outputs
# 1. 'infm_res.csv' gives us the data matching co-occurrences with drivers
# (interactions and enviornmental preference similarity) 
# 2. 'full_res.csv' gives us basic data on community diversity etc.
rule analyse_simulations:
    """
    From sign_ajd_res_log.csv produce infm_res.csv and full_res.csv
    """
    input:
        "simulation_data/{exp}/{tr}/{sim}/sign_ajd_res_log.csv"
    output:
        infm="simulation_data/{exp}/{tr}/{sim}/infm_res.csv",
        full="simulation_data/{exp}/{tr}/{sim}/full_res.csv"
    params:
        script="analysis_scripts/analysis.R"
    shell:
        """
        Rscript {params.script} {wildcards.exp} {wildcards.tr} {wildcards.sim}
        """


############################################################################
##### 6. Compile all data per 'experiment'
############################################################################
# Here we compile simulated, cloned, processed and analysed data into two 
# dataframes per experiment. The first contains all of the matching results
# (columns) for each simulation (row). The second gives us the community data
# metrics in the same format. Resulting dataframes are ready to be plotted.
# Use analysis_scripts/plot_figures.R to reproduce the plots used in our
# manuscript
rule compile_data:
    # only match real experiment names
    wildcard_constraints:
        exp="|".join(experiments.keys())
    input:
        infm=lambda wc: [
            f"simulation_data/{wc.exp}/{tr}/{sim}/infm_res.csv"
            for tr in experiments[wc.exp]
            for sim in experiments[wc.exp][tr]["simulations"]
        ],
        full=lambda wc: [
            f"simulation_data/{wc.exp}/{tr}/{sim}/full_res.csv"
            for tr in experiments[wc.exp]
            for sim in experiments[wc.exp][tr]["simulations"]
        ]
    output:
        infm="analysis_data/{exp}/infm_res_compiled.csv",
        full="analysis_data/{exp}/full_res_compiled.csv"
    run:
        os.makedirs(f"analysis_data/{wildcards.exp}", exist_ok=True)
        shell(f"Rscript analysis_scripts/compile_data.R {wildcards.exp}")