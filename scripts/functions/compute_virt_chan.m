function virt_chan = compute_virt_chan(source, data, voxel)
    % Compute a virtual channel (apply spatial filter)

    if isempty(voxel)
        % get the spatial filter
        spatial_filter = cat(1, source.avg.filter{:});
        virt_chan  =[];

        for ii = 1:length(data.trial)
            virt_chan.trial{ii} = spatial_filter * data.trial{ii};
        end


    else
         % get spatial filter for this voxel
        spatial_filter = source.avg.filter{voxel};

        for ii = 1:length(data.trial)
            % get the virtual sensor
            virt_chan.trial{ii} = spatial_filter * data.trial{ii};
        end
    end

    % create dummy labels
    for ii = 1:size(virt_chan.trial{1}, 1)  % first dim is number of channels
        label{ii} = num2str(ii);
    end

    % make it a FieldTrip-like structure
    virt_chan.time       = data.time;
    virt_chan.fsample    = data.fsample;
    virt_chan.trialinfo  = data.trialinfo;
    virt_chan.label      = label;
