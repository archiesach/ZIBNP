#!/bin/sh
#SBATCH --job-name=step1_sim2   # Job nam
#SBATCH --mail-type=ALL             # Mail events (NONE, BEGIN, END, FAIL, ALL
#SBATCH --mail-user=archiesachdeva@ufl.edu   # Where to send mail
#SBATCH --qos=s.guha #biostat-dept  #s.guha
#SBATCH --account=s.guha #biostat-dept #s.guha
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --mem-per-cpu=1gb           # Memory per processor
#SBATCH --time=01:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=array_%A-%a.out    # Standard output and error log
#SBATCH --array=1-30%15              # Array range

pwd; hostname; date

module load R

echo "Start ---- Data simulation"

R CMD BATCH dataGen_sim2.R

echo "End ---- Data simulation"

date
