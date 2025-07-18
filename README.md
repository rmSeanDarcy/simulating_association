############################################################################
##### Readme: SAMSARA workflow                                         #####
############################################################################

############################################################################
##### Simulating spatially explicit community data and analyising with multiple drivers of community assembly



Motivation:
This project is a project I have been working on for the last two years. We are currently finishing up the manuscript and I will not be revisiting the analysis. I think this is a great project to integrate into a snakemake workflow, because of the many interdependent steps in data simulation and analysis. The main issue is that if I were to implement my entire (so far manual) workflow, this would take way too much time. Choosing another project might be easier. But I would actually like to, one day make this simulation and analysis pipeline available and easy to use. To make your and my life easier I will therefore only integrate the data simulation and a single example for the analysis. I will write the code so the simulation can be run on the cluster. I am not a bioinformatician so I do not work with massive datasets where many files require the same downstream manipulation.


Project description:
The goal of this project was to evaluate the inference potential species 'co-occurrences' (essentially strong correlations) hold. We simulate population dynamics using generalised Lotka-Volterra differential equations (gLVs) in spatially explicit habitats. Species have pairwise interaction coefficients (suppressing or enhancing each others growth) and specific resource preferences. Depending on how well their preferences match the provided resources in a habitat, which are informed by there position in 'resource landscapes' they can reach higher or lower abundances. So both the resulting carrying capacitiy from matching preferences for and abundance of resources as well as their interactions dictate a species final abundance. The analysis simply foucsses on reconstructing our known inputs/drivers. So do Sp.1 and Sp.2 co-occurr because they have a mutualistic interaction coefficient, or is it because they share a similar resource preference? We also have additional confounding dynamics such as dispersal, implemented as simple diffusion dynamic between connected nodes in our habitat network. In this workflow I will only be simulating our baseline scenario recreating Fig.2b, c).


My previous workflow:
I have been working without git. I have simulated data and done preliminary analysis on the cluster. For this I generated slurm scripts in a python script and had to execute them by hand. One issue here is that due to stringing a bunch of different scripts behind one another some required heavy resources and some could be parallelised, while others required few resources but needed to be in series. In retrospect, these should be separate scripts, which could be perfectly run with individually specified resources in a workflow manager. For everything post simulation and preliminary analysis (null model testing etc) I did a lot of explorative data analysis for which the workflow manager is less optimal. This means that I had to adapt my final plotting function to work with any run experiments. This part I would not normally add to my snakemake workflow as users of this model would have their own intentions for analysing the data. At least for this automated snakemake workflow makes less sense.
One thing, which I messed up here also is that I often refer to working directories by paths. So for simulating and processing the data my first input is always the working directory, which seems quite unnecessary as I could always refer to the relative paths from my main working directory. This is why it is important to change the working directory in the 'config.yaml' file for the workflow to work...


My workflow in snakemake:
It will be a little rough for you to understand my code and the huge amount of parameters that are implemented in this analysis. 
Just a short run through of how it works:

The config file:
The basic logic of my simulations and analysis follows a separation into 'experiments' and 'treatments'. An experiment will often have a similar set of parameters, while a treatment alters a single one or few. In the config.yaml both are specified and treatement folders are generated within an experiment folder. Post simulation all treatment results are then combined in a single output file in the results folder. See the snakemake file for more info but in short: Every runs gets many inputs, which are informed by the config file: The working directory, the name of the experiment folder, the name of the treatment folder, the name of a predefined set of parameters (see ./simulation_code/Load_simulation_functions/set_standard_initial_conditions.py) and additional parameters that can be overwritten (often specified within treatments).

Rule 1.  simulates the community dynamics producing a large amount of .csv's contianing original settings, randomly generated parameters (such as interaction matrices) and the final simulated population data.

Rule 2. takes this data and calculates significant co-occurrences through null model testing, including a correlation coefficient threshold. This function also analyses the how well out inputs predicted these co-occurrences and some diversity metrics.

Rule 3. collects all of these data generated in 'treatments' into a single file that can be found in the 'experiment' folder.

Rule 4. then plots all treatments in an experiment in a barplot. The bars show the mean number of co-occurrences with additional colouring for the number of co-occurrences that could be either matched with interactions or/and similarity in environmental preference.

Rule 5. I added this last minute to fill up my 5 mandatory rules... Generates a very basic markdown file which can be seen as a report for the analysed experiments. The output is stored in the experiment folder.

Rule all for the final output.

(I only have five rules in this project right now because this is the way I set up my project. I know 5 was the required minimum... It would make most sense to split Rule 2. into two steps, one where I calculate the significant co-occurrences and another one in which I can do the analysis of matching inputs and co-occurrences. Sadly I lacked the time to do this because the code is very large and the way it is now it would have taken too long to separate it. But I will definately do this once I have more time!)


My opinion on this homework:
I really think my work has greatly benefitted from this exercise. Especially the fact that I can create a config file with all of my experiments and treatments and simply run a single script to generate all of the data in parallel is a game changer. I aim to publish my model at some time down the road and will most likely do so using snakemake. If the only thing a user has to do is create a config file, this really benefits user friendliness. Such a config file could also easily be generated by having users specify their parameters on an interactive website or so. There are many exciting possibilities I hadn't even thought of. Sadly due to time constraints I have also yet to run log/benchmark tests and optimise resource use. I have also not yet run my scripts on the cluster, which is of course the long term goal. For now, when running locally it makes sense to only simulate few replicates for simulation times not to be too long. I spread my documentation to here, which contains long winded explanations and a concise INSTALL.md file to recreate the workflow.
