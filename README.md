# toehold-exchange-walker-simulation
Kinetic Monte Carlo simulation of smFRET measurements of a DNA toehold-exchange walker as described in Li, Johnson-Buck, et al. Nature Nanotechnology 13, p. 723–729 (2018)

# Extended Description
Performs kinetic Monte Carlo simulation of single-molecule Forster resonance energy transfer (smFRET) measurements of a DNA toehold-exchange walker as described in Li, 
Johnson-Buck, et al. Nature Nanotechnology 13, p. 723–729 (2018)

The model models toehold binding and dissociation as well as branch migration processes using a kinetic Monte Carlo approach.

The output is a FRET-versus-time trajectory for a single simulated walker over a user-specified measurement interval.

The user can specify the following parameters:
nruns = the number of independent simulations to run
nfootholds = the number of footholds in the system
a = size of toehold domain of walker, in nucleotides
b = size of branch migration domain of walker, in nucleotides
tfinal = time of simulation, in seconds
tsegmentlength = segment size to break simulation into
itime = integration time of camera/measurement, in seconds

# License
This project is licensed under the BSD 3-Clause License.

# Author

Alex Johnson-Buck, 2018
