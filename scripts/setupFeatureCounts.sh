# setupFeatureCounts.sh
# Conda setup!

srun --partition=express --nodes=1 --cpus-per-task=2 --pty --time=00:60:00 /bin/bash
module load anaconda3/2021.05
conda create --prefix /work/geisingerlab/conda_env/subread python=3.9
# Deal with CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'.
conda init bash
source ~/.bashrc

conda activate /work/geisingerlab/conda_env/subread

# Biopython tweaks
# Make .condarc file, then tell conda to check bioconda.
conda config
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda install bioconda::subread