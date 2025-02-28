Shell commands to make my homework 4 work:

# Move to working directory of interest
$cd ???				# Replace with your working dir

# 1. Retrieve folder from github 
git clone https://github.com/rmSeanDarcy/SAMSARA.git
cd SAMSARA

# 2. Set up conda environment
conda env create -f homework4.yaml
conda activate homework4

### 3. !VERY IMPORTANT! Change the workdir in 'config.yaml' to your working directory!!!
workdir: "/your/path/workdir/here"
# This is an oversight on my part and I should make all of this independent of setting working dirs
# I would have had to change a lot of my code, but relative paths is the way to go!!!

# 3. Run the snakemake file
snakemake -c 6 		# Or whatever number of cores you want ('--cores all' also works)

# 4. If you want to make the DAG (already included is)
snakemake --forceall --dag | dot -Tpdf > dag.pdf
