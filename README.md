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

`singularity pull mpich.sif docker://docker.io/kvmkrao/mpich:latest`
 
 command to pull a opensource CFD software - OPenFOAM 
 \ 
(OpenFOAM is a C++ toolbox for the development of customized numerical solvers, and pre-/post-processing utilities for the solution of continuum mechanics problems, most prominently including computational fluid dynamics.) 
 
`docker pull docker.io/kvmkrao/kvmkrao/openfoam9:latest`

`singularity pull openfoam9.sif docker.io/kvmkrao/kvmkrao/openfoam9:latest`
 
# run a docker container
`singularity run anaconda.sif `

`docker run -it  --rm -v $PWD:/data -w /data kvmkrao:anaconda `

`docker run -it  openfoam9  /bin/bash`

`docker run -it  --rm -v $PWD:/data -w /data openfoam9  /bin/bash `

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

You can specify additional directories to bind mount into your container with the --bind option. 
For example, the data directory on the host system is bind mounted to the /mnt directory inside the container.

`$ singularity exec --bind /data:/mnt anaconda.sif `
