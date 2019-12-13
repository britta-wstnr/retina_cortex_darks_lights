% Project settings for the retina_cortex_onoff code repository

%% Paths and computers

[~, comp] = system('hostname');
comp = strtrim(lower(comp));
if strcmp(comp, 'mylocalname')
    host = '/path/on/local'
else
    host = '/path/on/cluster'
end

base_dir = fullfile(host, '/path/to/project/dir');
fig_dir = fullfile(base_dir, 'figures_out');

% add the subdirectory for all the functions
addpath ./functions

%% Toolbox paths

toolbox_path = '/path/to/toolboxes';
fieldtrip_path =  fullfile(toolbox_path, 'fieldtrip-2016-05-18');
fieldtrip_path_nmt = fullfile(toolbox_path, '/fieldtrip-20160718');

%% Filenames

data_on_erg_prep = 'dataon_erg_1500.mat';
data_off_erg_prep = 'dataoff_erg_1500.mat';

data_on_meg_fname = 'dataon_meg_1500_cleaned.mat';
data_off_meg_fname = 'dataoff_meg_1500_cleaned.mat';

data_on_erg_fname = 'dataon_erg_1500_cleaned2.mat';
data_off_erg_fname = 'dataoff_erg_1500_cleaned2.mat';

source_on_meg_fname = 'sourcetf_medcov_ON.mat';
source_off_meg_fname = 'sourcetf_medcov_OFF.mat';

vol_fname = 'vol_ft.mat';
lf_fname = 'leadfield_mags_openmeeg.mat';
sourcem_fname = 'sourcemodel_10mm.mat';
coreg_fname =  'nuts_nii.mat';

%% Templates

source_template_fname = fullfile(fieldtrip_path, ...
                                 'template/sourcemodel/', ...
                                 'standard_sourcemodel3d10mm.mat');
t1_template_fname = fullfile(fieldtrip_path_nmt, ...
                               'template/anatomy/single_subj_T1.nii');
atlas_fname = (fullfile(fieldtrip_path_nmt, ...
                                        '/template/atlas/aal/ROI_MNI_V4.nii'));

%% Set variables etc
subjs = {'0001','0002', '0003', '0004', '0005',  '0006', '0007', '0008', ...
         '0009','0010'};

erg_channel_sel = {'EOG001', 'EOG002'};
meg_channel_sel = {'MEGMAG' '-MEG2511'};

stim_defs = {'on', 'off'};

freq_bands = [
     55,  75
     75,  95
    105, 125
    125, 145
    155, 195];

% filter specifications for oscillatory potential, ON
filter_specs_on.dtf_freqs = [50 100 150];
filter_specs_on.hp_freq_on = 55;
filter_specs_on.transit_width = 6;

% filter specifications for oscillatory potential, OFF
filter_specs_off.dtf_freqs = [50 100 150];
filter_specs_off.bp_freq_off = [75 95];
filter_specs_off.transit_width = 6;

%% Colors for plotting

colors_rb = [   % colorbrewer2.org -- colorblind and photocopy friendly
    239, 138,  98
    178,  24,  43
    103, 169, 207
     33, 102, 172] / 255;

colors_granger = [
    118,  42, 131  % for the Granger spectra
     90, 174,  97
    ] / 255;

light_grey = [.9, .9, .9];
grey = [.6, .6, .6];
%% Load Fieldtrip if not yet loaded

try
    ft_defaults
catch
    warning('Fieldtrip is not on your path yet, adding it.');
    addpath(fieldtrip_path)
    ft_defaults
end

[ft_ver, ft_path] = ft_version;
display(sprintf('You are using Fieldtrip on path %s', ft_path));

%% Add other toolboxes needed

% add path for NutMEGtrip
addpath(fullfile(fieldtrip_path_nmt, 'contrib/nutmegtrip'));

% add path to boundedline for IQR shadings
addpath(fullfile(toolbox_path, 'boundedline'));

% add path to circstat toolbox
addpath(fullfile(toolbox_path, 'circstat'));

% toolbox for reslicing nifti files
% and creating the bem segmentation  - based on SPM8
addpath(fullfile(toolbox_path, 'smp8newseg/NIfTI_20140122'));

% add path to nutmeg and spm for coregistration
% ONLY ENABLE WHEN NEEDED! CAN CLASH WITH FIELDTRIP
% addpath(fullfile(toolbox_path, 'nutmeg'))
% addpath(fullfile(toolbox_path, 'spm8'))

% add path to iso2mesh
addpath(fullfile(toolbox_path, 'iso2mesh179'))