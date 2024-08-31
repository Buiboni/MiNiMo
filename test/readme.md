# Data
Download data at https://zenodo.org/records/13621818 (DOI: 10.5281/zenodo.13621818) and unzip the the nc files.
The script cannot read .ext file, one need to convert it through the command line:
```bash
lrs file.ext > file.ine
```
lrs can be install through ```apt-get```, or visit lrslib home page (http://cgm.cs.mcgill.ca/~avis/C/lrs.html).
# Prochlorococcus MED4
Carbon composition of the metabolic model was 37 mmolC
## Growth rate
To compute growth rate of prochlorococcus MED4 use this command line:
``` bash
python minimo.py -d {result_dir} {data_paper_dir}/diad_T.nc {data_paper_dir}/proc/bio_wo_iron/proc_orig_cake_lrs.ine {MiniMo_git_dir}/test/flux_proc 37
```

This will output growth.npy file at {result_dir} containing the growth of prochlorococcus MED4.

## Stress
To compute stress on PO4 fluxes:
``` bash
python minimo.py -s PPO4P -d {result_dir} {data_paper_dir}/diad_T.nc {data_paper_dir}/proc/bio_wo_iron/proc_orig_cake_lrs.ine {MiniMo_git_dir}/test/flux_proc 37
```

This will output growth.npy, and PPO4P_stress.npy file at {result_dir}.

## Auxiliary
To compute overproduction of amino acid by prochlorococcus MED4:
``` bash
python minimo.py -s PPO4P -a {data_paper_dir}/proc/niche_biomass/amino/proc_orig_cake_amino_lrs.ine -d {result_dir} {data_paper_dir}/diad_T.nc {data_paper_dir}/proc/bio_wo_iron/proc_orig_cake_lrs.ine {MiniMo_git_dir}/test/flux_proc 37
```

This will output growth.npy, PPO4P_stress.npy, and auxiliary.npy file at {result_dir}.

# Phaeodactylum Tricornutum
Carbon composition of the metabolic model was 45 mmolC:
``` bash
python minimo.py -s PPO4D -a {data_paper_dir}/phaeo/dmsp/Model_iLB1034_cnp_pisces_cake_gam_lrs.ine -d {result_dir} {data_paper_dir}/diad_T.nc {data_paper_dir}/phaeo/dmsp/Model_iLB1034_cnp_pisces_cake_gam_lrs_bio.ine {MiniMo_git_dir}/flux_diat 45

```

This will output growth.npy, PPO4D_stress.npy, and auxiliary.npy file at {result_dir}. Here the auxialiary is the DMSP production, and the growth needs to be corrected by a coefficient of 0.03333, as the niche is computed on the GAM rather than the biomass reaction.

# Thalassiosira Pseudonana
Carbon composition of the metabolic model was 40 mmolC:
``` bash
python minimo.py -s PPO4D -a {data_paper_dir}/thal/niche_biomass/faa/Thaps_HL_clean_cake_si_bio_comp_lrs.ine -d {result_dir} {data_paper_dir}/diad_T.nc {MiniMo_git_dir}/flux_diat_si 40
```

This will output growth.npy, PPO4D_stress.npy, and auxiliary.npy file at {result_dir}. Here the auxialiary is the free amino acid production.