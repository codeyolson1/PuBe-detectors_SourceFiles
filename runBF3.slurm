#!/bin/tcsh
#SBATCH --time=12:00:00
#SBATCH -J runG4bay-BF3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH -o geant4BF3bay_%A.out
#SBATCH --account=coe-np
#SBATCH --partition=coestudent-np
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=codey.olson@utah.edu

# Set up my directories
setenv SCRATCHDIR /scratch/general/lustre/$USER/$SLURM_JOB_ID
setenv CURRENTDIR $HOME/Geant4/gitRepos/PuBe-detectors_SourceFiles
mkdir -p $SCRATCHDIR
mkdir -p $CURRENTDIR/outputsBF3/outputBF3bay-$SLURM_JOB_ID

# Setup the environment
cd
source custom1.csh
cd $CURRENTDIR

# Copy my data into the scratch directory:
cp -r $CURRENTDIR/build $SCRATCHDIR
cp -r $CURRENTDIR/macros $SCRATCHDIR

# Move into my scratch directory and run the simulation
cd $SCRATCHDIR/build
./pu-be runBF3.mac bf3

# Copy over the output of my file to my home directory
cp BF3Response.root $CURRENTDIR/outputsBF3/outputBF3bay-$SLURM_JOB_ID
cp BF3Response.root-conv.txt $CURRENTDIR/outputsBF3/outputBF3bay-$SLURM_JOB_ID

# Copy the stdouput to the output folder
cd $CURRENTDIR
cp geant4BF3bay_$SLURM_JOB_ID.out $CURRENTDIR/outputsBF3/outputBF3bay-$SLURM_JOB_ID
rm geant4BF3bay_$SLURM_JOB_ID.out
