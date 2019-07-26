function filtered_data = filter_osc_pot(data, condition, timelock, filter_specs)
    % Filtering for oscillatory potentials

    cfg = [];
    cfg.dftfilter = 'yes';
    cfg.dftfreq = filter_specs.dtf_freqs;
    cfg.demean = ' yes';

    if strcmp(condition , 'on')
        if(0)
            cfg.hpfilter = 'yes';
            cfg.hpfreq = hp_freq_on;
        else
            cfg.hpfilter = 'yes';
            cfg.hpfreq = filter_specs.hp_freq_on;
            cfg.hpfilttype = 'firws';
            cfg.hpfiltdf = filter_specs.transit_width;
            cfg.hpfiltdir = 'onepass-zerophase';
            cfg.hpfiltwintype = 'hann';
        end
    elseif strcmp(condition, 'off')
            cfg.hpfilter = 'yes';
            cfg.hpfilttype = 'firws';
            cfg.hpfreq = filter_specs.bp_freq_off(1);
            cfg.lpfiltdf = filter_specs.transit_width;
            cfg.lpfilter = 'yes';
            cfg.lpfilttype = 'firws'
            cfg.lpfreq = filter_specs.bp_freq_off(2);
            cfg.hpfiltdir = 'onepass-zerophase';
            cfg.hpfiltdf = filter_specs.transit_width;
            cfg.hpfiltwintype = 'hann';
            cfg.lpfiltdir = 'onepass-zerophase';
            cfg.lpfiltwintype = 'hann';
    else
        error(sprintf('Supported conditions are on and off, got %s', ...
                      condition));
    end

    tmp_data = ft_preprocessing(cfg, data);

    % avaerage over trials if needed:
    if strcmp(timelock, 'yes')
        filtered_data = ft_timelockanalysis([], tmp_data);
    else
        filtered_data = tmp_data;
    end
