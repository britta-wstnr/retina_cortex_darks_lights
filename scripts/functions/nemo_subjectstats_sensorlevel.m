% Extends the Hilbert stats to group level

% this is a variant of
% https://github.com/meeg-cfin/nemolab/blob/master/hilbertsourceloc/nemo_hilbertstats_sensorlevel.m

% should best be used with nemo_hilbert_sensorlevel


% this is a hack that only works if all trials have same time axis! FIX THIS.
baselinewindow_idx = dsearchn(data{1}.time{1}',baselinewindow');

testidx = 1 : length(data{1}.time{1});

%% run subject level stats

disp(freq_bands(ii ,:));
tvals = ones(length(data{ii}.trial), numel(data{ii}.label), ...
             length(data{ii}.time{1}));

for jj=1:length(data{ii}.trial)
    tvals(jj, :, :) = data{ii}.trial{jj};
end

% create a dummy tlk structure to use for output storing
data_tested{ii} = ft_timelockanalysis([], data{ii}); % initialize source_hilb

p = ones(size(tvals, 2), size(tvals, 3));
p(:) = NaN;
zval = zeros(size(tvals, 2), size(tvals, 3));

ft_progress('init','etf');
% loop over channels (jj) and timepoints (kk)
for jj=1:size(tvals,2)
    tval_bl = squeeze(tvals(:, jj, ...
                            baselinewindow_idx(1):baselinewindow_idx(2)));
    tval_bl = tval_bl(:);
    for kk=1:length(testidx)
        tt = testidx(kk);
        % all subjects, jj=channel, tt=timepoints that really should be tested
        [p(jj, kk), ~, stats] = ranksum(tvals(:, jj, kk), tval_bl);
        zval(jj, kk) = stats.zval;
    end
    ft_progress(jj/size(tvals, 2),'%d of %d', jj, size(tvals, 2));
end
ft_progress('close');

data_tested{ii}.stat = zval;
data_tested{ii}.pval = log10(p);
