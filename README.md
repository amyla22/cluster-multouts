# cluster-multouts

## Simulations
### To run a single simulated data set with estimated parameters:
MO_simulated_data.R = main file 
 - This file sources the remaining three .R files to run the method on a single simulated data set and output results
 
MO_simulated_data.R 
 - This file contains the values and parameters defining the simulation setting (a) on 20 outcomes with 7 domains from the paper
 
SimulateData.R 
 - This file contains a function to simulate data similar the SCDS data
 
clustMultOutwImp_forsims.R
 - This file contains a function to implement the clustering procedure outlined in the paper

### To run a full simulation similar to paper:
Submitting sim20.500.7.a.sbatch to SLURM will source the following files:
 - SetParam.R = file containing true parameter values and domain arrangements
 - SimRun.R = file reading in arguments from .sbatch to run 100 simulated datasets
 - SimulateData.R 
 - clustMultOutwImp_forsims.R
This will produce 100 separate files summarizing results of each of the 100 simulated datasets

Running CompileFiles.R will compile the 100 tables of results into a single table. The tables printed at the end of the script match the summary of the simulations in Section 4. Specifically, the output in the first table will match the summary of setting (a) in Table 4.

