% Clean MEG data manually (data browser to look at trials)

project_settings;
for nn = 1:length(subjs)

    project_settings;
    proc_dir = fullfile(base_dir, subjs{nn});

    % Load ERG data
    load(fullfile(proc_dir, data_on_erg_fname));

    % Load data containing MEG data
    load(fullfile(proc_dir, 'data_clean.mat'));

    %% get off data
    switch subjs{nn}
        case {'0001' '0002' '0003' '0004' '0005' '0006' '0010'}
            clear cfg
            cfg.channel = {'MISC001'};
            photodiode = ft_selectdata(cfg, data_clean);

            for jj=1:numel(data_clean.trialinfo)

                photodiode.trial{jj} = smooth(photodiode.trial{jj},20)';

                indices = photodiode.trial{jj} >=0.025;
                strindices = num2str(indices,'%-d');
                tmp = strfind(strindices, '10');
                offind(jj) = tmp(1);
            end

            %% Cut data in trials - homemade

            offindreal = data_clean.cfg.trl(:, 1) + offind' - 1;

            clear cfgtrig
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

    %% Downsample the MEG data

    cfg = []
    cfg.resamplefs = 1500;
    dataonres = ft_resampledata(cfg, data_clean);
    dataoffres = ft_resampledata(cfg, dataoff);

    %% Now find which trials to cut in the off data
    if(length(data_clean.trial) == ...
       size(dataon_clean2.cfg.previous.previous.previous.trl, 1))

        if(0)
            %% Get the trl from the MEG data

            % get the old trl (downsampled data)
            dummy = ft_databrowser([], dataonres);
            datatrl = ft_rejectartifact(dummy, dataonres);


            trl_good = dataon_clean2.cfg.trl;
            trl_old = datatrl.cfg.trl;

            for ii=1:length(trl_good)
                tmp =    find(trl_good(ii, 1) == trl_old(:, 1));
                if(size(tmp, 1) > 0)
                    goodtrials(ii) = tmp;
                end
            end
        else
            goodtrials = dataon_clean2.cfg.trials;
        end
    else
        error('Something with the data is wrong, check your trials')
    end

    %%  Clean MEG data

    cfg = [];
    cfg.trials = goodtrials;
    dataon_clean = ft_selectdata(cfg, dataonres);
    dataoff_clean = ft_selectdata(cfg, dataoffres);

    % save:
    savepath = fullfile(base_dir, subjs{nn}, data_on_meg_fname);
    save(savepath, 'dataon_clean', '-v7.3');

    savepath = fullfile(base_dir, subjs{nn}, data_off_meg_fname);
    save(savepath, 'dataoff_clean', '-v7.3');

    keep subjs nn
end
