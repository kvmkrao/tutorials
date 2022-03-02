# Tutorials

On HPC cluster, load singularity module 
\
`module load singularity `

# pull a docker image 

command to pull Anaconda3-2021.11-Linux-x86_64 version (The World's Most Popular Data Science Platform)

`docker pull docker.io/kvmkrao/anaconda:latest

`singularity pull anaconda.sif docker://docker.io/kvmkrao/anaconda:latest`

command to pull MPICH3.2 (High-Performance Portable MPI)

`docker pull docker.io/kvmkrao/mpich:latest`

`singularity pull mpich.sif docker://docker.io/kvmkrao/mpich:latest
 
# run a docker container
`singularity run anaconda.sif `

`docker run -it  --rm -v $PWD:/data -w /data kvmkrao:anaconda `

Slurm script for a serial application:
```
#!/bin/bash
#SBATCH --job-name=serial        # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=4G                 # total memory per node (4 GB per cpu-core is default)
#SBATCH --time=00:05:00          # total run time limit (HH:MM:SS)

module load singularity
singularity pull anaconda.sif docker://docker.io/kvmkrao/anaconda:latest
singularity run anaconda.sif
```

Slurm script for parallel MPI applications 
```
#!/bin/bash
#SBATCH --job-name=parallel      # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=4               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G per cpu-core is default)
#SBATCH --time=00:05:00          # total run time limit (HH:MM:SS)

module load singularity
module load openmpi/gcc/3.1.5/64
srun singularity exec path/to/anaconda.sif
```
