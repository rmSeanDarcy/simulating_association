Shell commands to make my homework 4 work:

# Move to working directory of interest
$cd ???				# Replace with your working dir

# 1. Retrieve folder from github 
git clone https://github.com/rmSeanDarcy/SAMSARA.git
cd SAMSARA

# 2. Set up conda environment
conda env create -f homework4.yaml
conda activate homework4

# 3. Run the snakemake file
snakemake - 6 		# Or whatever number of cores you want ('--cores all' also works)

# 4. If you want to make the DAG (already included is)
snakemake --forceall --dag | dot -Tpdf > dag.pdf
