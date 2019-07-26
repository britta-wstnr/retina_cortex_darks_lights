% Run a beamformer on preloaded data

% this script should be called from within another script

%% set subject specific paths

in_dir      = fullfile(base_dir, subjs{nn});
headmodel_dir = fullfile(base_dir, subjs{nn}, 'bem');

mri_dir = fullfile(base_dir, 'mri', ['s', subjs{nn}]);
rmri_fname = ['s', subjs{nn}, '_resl.nii'];

%% load precomputed steps

% load the headmodel
load(fullfile(headmodel_dir, 'vol_ft.mat'));

% load template source-models (MNI and individual) and convert to mm
load(source_template_fname);
template_grid = ft_convert_units(sourcemodel, 'mm');

load(fullfile(in_dir, sourcem_fname));
sourcemodel = ft_convert_units(sourcemodel, 'mm');

% load coregistration
load(fullfile(in_dir, 'nuts_nii.mat'));
coreg = nuts.coreg;

% load leadfield
load(fullfile(headmodel_dir, lf_fname));

%% coregister the data

% we are dealing with neuromag data for this project:
orig_coordsys = 'neuromag';

[mri2meg_tfm, coordsys] = ft_headcoordinates(coreg.fiducials_mri_mm(3, :), ...
                                             coreg.fiducials_mri_mm(1, :), ...
                                             coreg.fiducials_mri_mm(2, :), ...
                                             orig_coordsys);

meg2mri_tfm = inv(mri2meg_tfm);

% convert grad structure to MRI coordinate system, make a copy to not overwrite
% and coregister same data twice
mri_grad = data_stim.grad;
mri_grad = ft_convert_units(mri_grad, 'mm');

mri_grad = ft_transform_sens(meg2mri_tfm, mri_grad);
mri_grad.coordsys = 'spm';

% plug into the datasets:
data_coreg = data_stim; data_coreg.grad = mri_grad;

% also select channels:
cfg = [];
cfg.channel = meg_channel_sel;
data_coreg = ft_selectdata(cfg, data_coreg);

%% data covariance matrix and spatial filter

% this is a dummy:
cfg = [];
cfg.covariance = 'yes';
cov_mat = ft_timelockanalysis(cfg, data_coreg);

% compute robust covariance matrix:
[med_cov, ~] = compute_median_cov(data_coreg, leadfield, meg_channel_sel);
cov_mat.cov = med_cov;

% spatial filter for virtual channels
cfg = [];
cfg.method = 'lcmv';
cfg.grid = leadfield;
cfg.vol = vol;
cfg.lcmv.reducerank = 'no';
cfg.lcmv.fixedori  = 'yes';  % orientation selection (max. power)
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.weightnorm  = 'nai';  % unit-noise-gain BF with noise normalization
cfg.lcmv.keepfilter = 'yes';
source_out = ft_sourceanalysis(cfg, cov_mat);
