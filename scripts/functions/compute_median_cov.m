function [med_cov, chan_id] = compute_median_cov(data, leadfield, channel_sel)
    % Computes a robust covariance matrix based on the Minimum Covariance
    % Determinant estimator

    % take care of channel selections
    meg_labels = ft_channelselection(channel_sel, leadfield.label);
    [~, chan_id] = intersect(data.label, meg_labels);

    % prepare computation
    data_out = [];
    for jj = 1:length(data.trial)
        data_out.b(:, :, jj) = real(data.trial{jj}(chan_id, :));
    end

    n_samples = size(data_out.b,2);
    data_out.b = reshape(data_out.b, size(data_out.b, 1), ...
                         size(data_out.b, 2) * size(data_out.b, 3));

    % compute covariance
    med_cov(chan_id, chan_id) = robustcov(data_out.b');
