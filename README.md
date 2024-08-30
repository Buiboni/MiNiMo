# MiNiMo
Compute growth rates of an organism through its metabolic niche.
## Environment
To be able to run the script your environment need to have numpy and netCDF4 python packages.
Numpy and netCDF4 are available through pip, or conda.
## Usage
This script is aimed to be used with environmental conditions and a metabolic niche description
```
python minimo.py --help
usage: minimo.py [-h] [-d DIRECTORY] [-s STRESS] [-a AUX] environmental_fluxes niche_file flux_name carbon

positional arguments:
  environmental_fluxes  File with environmental fluxes in HDF4 format (.nc), temporal unit is second
  niche_file            File with the niche description in .ine (h-format), see lrslib for detail on the format. Unit of flux should be in molXX.gDW-1.h-1
  flux_name             File where flux names are stored, one per line, corresponding to the niche fluxes.
  carbon                Carbon composition of the considered system in mole. In a GSM, this is the stoichiometry of C in the biomass reaction.

options:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        Directory where result will be written, if not provided nothing will be written.
  -s STRESS, --stress STRESS
                        Nutrient on which to compute the stress
  -a AUX, --aux AUX     Auxiliary metabolite niche file computation


```
The outputs of the script will be in DIRECTORY: growth.npy, stress.npy and auxiliary.npy.