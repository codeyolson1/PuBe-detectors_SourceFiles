#!/bin/tcsh
#SBATCH --time=18:00:00
#SBATCH -J runG4bay-He3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH -o geant4He3bay_%A.out
#SBATCH --account=coe-np
#SBATCH --partition=coestudent-np
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=codey.olson@utah.edu

# Set up my directories
setenv SCRATCHDIR /scratch/general/vast/$USER/$SLURM_JOB_ID
setenv CURRENTDIR $HOME/Geant4/gitRepos/PuBe-detectors_SourceFiles
mkdir -p $SCRATCHDIR
mkdir -p $CURRENTDIR/outputsHe3/outputHe3bay-$SLURM_JOB_ID

# Copy data to the scratchdirectory:
cd $HOME/Geant4/
cp -r data $SCRATCHDIR
cd $CURRENTDIR

# Setup the environment
cd
source g4ComputeEnv.csh
cd $CURRENTDIR

# Copy my data into the scratch directory:
cp -r $CURRENTDIR/build $SCRATCHDIR
cp -r $CURRENTDIR/macros $SCRATCHDIR

# Move into my scratch directory and run the simulation
cd $SCRATCHDIR/build
./pu-be runHe3-4ft.mac he3

# Copy over the output of my file to my home directory
cp He3Response.root $CURRENTDIR/outputsHe3/outputHe3bay-$SLURM_JOB_ID
cp He3Response.root-conv.txt $CURRENTDIR/outputsHe3/outputHe3bay-$SLURM_JOB_ID

# Copy the stdouput to the output folder
cd $CURRENTDIR
cp geant4He3bay_$SLURM_JOB_ID.out $CURRENTDIR/outputsHe3/outputHe3bay-$SLURM_JOB_ID
rm geant4He3bay_$SLURM_JOB_ID.out

# Delete my data from SCRATCHDIR
rm -rf $SCRATCHDIR
