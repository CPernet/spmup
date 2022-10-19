# system tests

Those can also be used as demos.

## SPM face repetition dataset

Requires bids-matlab to run.

The main script will install it using git if bids-matlab is not in the matlab
path.

## ds000117

URL: https://openneuro.org/datasets/ds000117/versions/1.0.5/

Install with [datalad](http://handbook.datalad.org/en/latest/).

```bash
datalad install https://github.com/OpenNeuroDatasets/ds000117.git
```

Download the MRI data (T1w anat, all fieldmaps and functional data) of the 3
first subjects.

```bash
cd ds000117
datalad get sub-0[1-3]/ses-mri/anat/*T1w*
datalad get sub-0[1-3]/ses-mri/f*
```
