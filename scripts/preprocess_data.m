% Loading and preprocessing of all data

project_settings;

for nn = 1:length(subjs)

    proc_dir = 'long_subj_dir_identifier'

    meg_file_base = ['proj_identifier', subjs{nn}];
    meg_file_ending = '.fif';
    meg_file_num = {'001', '002', '003'};

    time = [];
    trial = [];

    % concatenate MEG files
    for jj = 1:3

        meg_file = fullfile(proc_dir, ...
                            [meg_file_base, meg_file_num{jj}, ...
                            meg_file_ending]);

        cfg = [];
        cfg.dataset = meg_file;
        cfg.demean = 'yes';
        datafull = ft_preprocessing(cfg);
        time = [time, datafull.time{1}];
        trial = [trial, datafull.trial{1}];
    end

    datafull.time{1}  = time;
    datafull.trial{1} = trial;

    % trigger business
    switch subjs{nn}
        case {'0001' '0002' '0003' '0004' '0005' '0006' '0010'}

            %% Photodiode stuff
            clear cfg
            cfg.channel = {'MISC001'};
            photodiode = ft_selectdata(cfg, datafull);

            photodiode.trial{1} = smooth(photodiode.trial{1},20)';

            indices = photodiode.trial{1} >= 0.02;
            strindices = num2str(indices,'%-d');
            onind = strfind(strindices, '01');
            onind = onind + 1;
            offind = strfind(strindices, '10');

            % correct onind for presentation screen at the end
            onind = onind(1:250);

            %% Cut data in trials - homemade

            clear cfgtrig
            cfgtrig.trl = [onind' - 5000 ...
                           onind' + 5000 ...
                           repmat(-5000, ...
                                  length(onind),1) repmat(1, length(onind),1)];
            dataon = ft_redefinetrial(cfgtrig, datafull);

        case {'0007' '0008' '0009'}

            %% Trigger stuff
            clear cfg
            cfg.channel = {'STI101'};
            trigger = ft_selectdata(cfg, datafull);

            indices = trigger.trial{1} >= 150;
            strindices = num2str(indices,'%-d');
            onind = strfind(strindices, '01');
            onind = onind + 1;
            offind = strfind(strindices, '10');

            % correct onind for presentation screen at the end
            onind = onind(1:250);
            % correct onind for the trigger delay present in the data:
            onind = onind + round(0.01787*5000);

            %% Cut data in trials - homemade
            clear cfgtrig
            cfgtrig.trl = [onind' - 5000 ...
                           onind' + 5000 ...
                           repmat(-5000, ...
                                  length(onind),1) repmat(1, length(onind),1)];
            dataon = ft_redefinetrial(cfgtrig, datafull);
    end

    %% Clean this

    % some quick filter just for the databrowser display
    cfg = []
    if(strcmp(subjs{nn}, '0001'));
        cfg.channel = {'EOG001', 'EOG002', 'EOG003'};
    else;
        cfg.channel = {'EOG001', 'EOG002', 'EOG003', 'EOG004'};
    end
    cfg.resamplefs = 1000;
    datares = ft_resampledata(cfg, dataon);

    cfg = removefields(cfg, 'resamplefs');

    cfg.hpfilter = 'yes';
    cfg.hpfreq =1;
    cfg.hpfilttype = 'but';
    cfg.demean = 'yes';
    datadisplay = ft_preprocessing(cfg, datares);
    clear datares

    %% Databrowser
    cfg = [];
    cfg.method = 'butterfly';
    if(strcmp(subjs{nn}, '0001'))
        cfg.channel  = {'EOG001', 'EOG002', 'EOG003'};
        cfg.mychan   = {'EOG003'};
    elseif(strcmp(subjs{nn}, '0008'))
        % for subject 0007: EOG001 came loose; discarded from further analysis
        cfg.channel  = {'EOG002', 'EOG003', 'EOG004'};
        cfg.mychan   = {'EOG003', 'EOG004'};
    else
        cfg.channel  = {'EOG001', 'EOG002', 'EOG003', 'EOG004'};
        cfg.mychan   = {'EOG003', 'EOG004'};
    end

    cfg.mychanscale = 5;
    cfgartf = ft_databrowser(cfg, datadisplay);
    dataon_cl_eog = ft_rejectartifact(cfgartf, dataon);

    %% additionally check the MEG
    % downsample just for display:

    cfg = []
    cfg.resamplefs = 1000;
    datares = ft_resampledata(cfg, dataon_cl_eog);

    cfg  = [];
    cfg.method = 'butterfly';
    cfg.channel = {'MEG'};

    cfgartf = ft_databrowser(cfg, datares);
    data_clean = ft_rejectartifact(cfgartf, dataon_cl_eog);

    savepath = fullfile(base_dir, subjs{nn}, 'data_clean');
    save(savepath, 'data_clean', '-v7.3');

    %% get the OFF data

    switch subjs{nn}
      case {'0001' '0002' '0003' '0004' '0005' '0006' '0010'}

            clear cfg
            cfg.channel = {'MISC001'};
            photodiode = ft_selectdata(cfg, data_clean);

            for jj = 1:numel(data_clean.trialinfo)

                photodiode.trial{jj} = smooth(photodiode.trial{jj}, 20)';

                indices = photodiode.trial{jj} >= 0.025;
                strindices = num2str(indices,'%-d');
                tmp = strfind(strindices, '10');
                offind(jj) = tmp(1);
            end

            %% Cut data in trials - homemade

            offindreal = data_clean.cfg.trl(:,1)+offind'-1;

            clear cfgtrig
            % cfgtrig.demean= 'yes';
            % cfgtrig.baselinewindow = [-0.4 -0.05];
            cfgtrig.trl = [offindreal - 5000 ...
                           offindreal + 2000 ...
                           repmat(-5000, length(offindreal), 1) ...
                           repmat(1, length(offindreal), 1)];
            dataoff = ft_redefinetrial(cfgtrig, data_clean);

        case {'0007' '0008' '0009'}

            clear cfg
            cfg.offset = -.48*5000;
            dataoff = ft_redefinetrial(cfg, data_clean);
    end

    cfg = [];
    cfg.channel = {'EOG001', 'EOG002'};
    dataon = ft_selectdata(cfg, data_clean);
    dataoff = ft_selectdata(cfg, dataoff);

    cfg = []
    cfg.channel = {'EOG001', 'EOG002'};
    cfg.resamplefs = 1500;
    datares = ft_resampledata(cfg, dataon);
    dataoffres = ft_resampledata(cfg, dataoff);

    dataon = datares;
    dataoff = dataoffres;

    save(fullfile(procdir, data_on_erg_prep), 'dataon', '-v7.3');
    save(fullfile(procdir, data_off_erg_prep), 'dataoff', '-v7.3');

    clear datares dataoffres

end
