# Darks elicit gamma-band activity in the human retina and visual cortex with narrower bandwidth and faster propagation time than lights

Britta U. Westner & Sarang S. Dalal


## Repository
This repository contains the data analysis scripts for Darks elicit gamma-band activity in the human retina and visual cortex with narrower bandwidth and faster propagation time than lights.
_This paper is submitted, a pre-print can be found_ [_here_](https://www.biorxiv.org/content/10.1101/153551v1).


The data analysis scripts are organized as follows:

### Configuration file
* `project_settings.m`

### Processing files
Preprocess and clean the data:
* `preprocess_data.m`
* `clean_erg.m`
* `clean_meg.m`

Prepare source localization:
* `compute_headmodels.m`
* `compute_leadfields.m`

### Generate figures
Figure 1
* `plot_evoked.m`

Figure 2
* `plot_evoked_hfo.m`
* `plot_erg_itc.m` together with `hilbert_erg.m`

Figure 3
* `plot_evoked_hfo.m`

Figure 4
* `hilbert_meg.m` together with `prep_hilbert_meg.m`

Figure 5
* `plot_itc_traces.m`

Figure 6
* `plot_itc_times.m`

### Processing functions
The folder `functions` contains several functions and scripts that are called by the scripts listed above.


## Dependencies
A full list of dependencies can be found at the end of `project_settings.m`.
This work uses [FieldTrip](https://github.com/fieldtrip/) for MEG and ERG data analysis.
* FieldTrip, version 2016-05-18, for data analysis
* FieldTrip, version 2016-07-18 with updated module `nutmegtrip` for plotting
