# uFWD
A collections of scripts to run unsupported forward models in SPM.

### Motivaiton
SPM forward modelling relies on what is available in FieldTrip, this framework allows to test/design forward models not included in Fieldtrip whilst not interfereing with the functionality of SPM or FT.

### Available models

+ uBEM: an M/EEG 3-shell Boundary Element Model. Based on the model used in MNE-Python, it gives numerically identical reults. (EEG support in uFWD coming immenently).

## Usage

To use after initialising SPM with `spm eeg` add the repository to your MATLAB path with `addpath /path/to/uFEM`
