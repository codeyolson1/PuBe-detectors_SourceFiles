#!/bin/tcsh
#SBATCH --time=0:30:00
#SBATCH -J runG4bay-He3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH -o geant4He3bay_%A.out
#SBATCH --account=coe-np
#SBATCH --partition=coestudent-np
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=codey.olson@utah.edu

# Setup the environment
source custom1.csh

# Set up my directories
setenv SCRATCHDIR /scratch/general/lustre/$USER/$SLURM_JOB_ID
setenv CURRENTDIR $HOME/Geant4/gitRepos/PuBe-detectors_SourceFiles
mkdir -p $SCRATCHDIR
mkdir -p $CURRENTDIR/outputHe3bay-$SLURM_JOB_ID

# Copy my data into the scratch directory:
cp -r $CURRENTDIR/build $SCRATCHDIR
cp -r $CURRENTDIR/macros $SCRATCHDIR

# Move into my scratch directory and run the simulation
cd $SCRATCHDIR/build
./pu-be runHe3.mac he3

# Copy over the output of my file to my home directory
cp He3Response.root $CURRENTDIR/outputHe3bay-$SLURM_JOB_ID
cp He3Response.root-conv.txt $CURRENTDIR/outputHe3bay-$SLURM_JOB_ID

# Copy the stdouput to the output folder
cd
cp geant4He3bay_$SLURM_JOB_ID.out $CURRENTDIR/outputHe3bay-$SLURM_JOB_ID
rm geant4He3bay_$SLURM_JOB_ID.out

# Delete my data from SCRATCHDIR
rm -rf $SCRATCHDIR
