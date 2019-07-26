% Clean ERG data manually (data browser to look at trials)

project_settings;

% manual cleaning, so no loop:
nn = 10;

proc_dir = fullfile(base_dir, subjs{nn});

%% Load the preprocessed data

load(fullfile(proc_dir, data_on_erg_prep));
load(fullfile(proc_dir, data_off_erg_prep));

%% Filter for timelock

cfg = []
cfg.channel = {'EOG001', 'EOG002'};
cfg.hpfilter = 'yes';
cfg.hpfreq =5;
cfg.hpfilttype = 'firws';
cfg.hpfiltwintype = 'kaiser';
cfg.hpfiltdir = 'onepass-zerophase';
cfg.hpfiltdf = 4;
cfg.plotfiltresp = 'no';
cfg.demean = 'yes';
data_clean_filt = ft_preprocessing(cfg, dataon);
data_clean_filt_off = ft_preprocessing(cfg, dataoff);

%% Clean data again

dummy = ft_databrowser([], data_clean_filt);
datatrl = ft_rejectartifact(dummy, data_clean_filt);

%%  Data browser

cfg = [];
cfg.method = 'butterfly';
if(strcmp(subjs{nn}, '0008'))
    cfg.channel  = {'EOG002'};
else
    cfg.channel  = {'EOG001', 'EOG002'};
end
cfgartf = ft_databrowser(cfg, data_clean_filt);
dataon_clean2 = ft_rejectartifact(cfgartf, data_clean_filt);

%% Now find which trials to cut in the off data

trl_good = dataon_clean2.cfg.trl;
trl_old = datatrl.cfg.trl;

for ii = 1:length(trl_good)
    tmp =    find(trl_good(ii,1)==trl_old(:,1));
    if(size(tmp,1)>0)
        goodtrials(ii) = tmp;
    end
end

%%  "Clean" off data

cfg = [];
cfg.trials = goodtrials;
dataoff_clean2 = ft_selectdata(cfg, data_clean_filt_off);

%% save cleaned data

% "undo" the filtering:
cfg = [];
cfg.trials = goodtrials;
dataon_clean2 = ft_selectdata(cfg, dataon);
dataoff_clean2 = ft_selectdata(cfg, dataoff);

% save:
savepath = fullfile(proc_dir, data_on_erg_fname);
save(savepath, 'dataon_clean2', '-v7.3');

savepath = fullfile(proc_dir, data_off_erg_fname);
save(savepath, 'dataoff_clean2', '-v7.3');
