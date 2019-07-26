% File containing peak times and propagations times

% For peak identification, check plot_itc_traces.m

%% peaks plotting

on_retina = [
    NaN, 0.022, 0.024, NaN, 0.029, 0.034, 0.020, 0.039, 0.050, 0.066;
    0.048, 0.028, 0.034, 0.026, NaN, 0.032, 0.027, 0.037, 0.035, 0.048;
    0.029 0.027 NaN 0.027 0.025 0.026 0.037 0.026 0.024 0.031;
    0.013 0.029 0.029 0.022 0.030 0.027 0.023 0.027 0.027 0.025
    ] * 1000;

on_brain = [
    0.081, 0.069, 0.053, 0.058, 0.071, 0.074, 0.056, 0.047, NaN, 0.114;
    0.070, 0.066, 0.059, 0.096, 0.065, 0.084, 0.096, 0.046, 0.073, 0.087;
    0.055, 0.078, 0.065, 0.044, 0.084, 0.079, 0.082, 0.053, 0.076, 0.089;
    0.045, 0.071, 0.058, 0.057, 0.059, 0.054, 0.068, 0.051, NaN, 0.071;
    ] * 1000;

off_retina = [
    0.034, NaN, NaN, 0.040, 0.029, 0.033, 0.045, 0.026, NaN, 0.039] * 1000;

off_brain = [
    0.055, 0.059, 0.052, 0.047, 0.045, 0.070, 0.063, 0.067, 0.043, 0.065
            ] * 1000;

%% propagation times ON

for ii = 1:4
    prop_times_on(ii, :) = on_brain(ii, :) - on_retina(ii, :);
end

nanmedian(prop_times_on, 2);
iqr(prop_times_on, 2);

%% propagation times OFF

prop_times_off = off_brain - off_retina;

nanmedian(prop_times_off, 2);
iqr(prop_times_off, 2);

%% IQR of other data
display('ON retina')
nanmedian(on_retina, 2)
iqr(on_retina, 2)

display('OFF retina')
nanmedian(off_retina)
iqr(off_retina)

display('ON brain')
nanmedian(on_brain, 2)
iqr(on_brain, 2)

display('OFF brain')
nanmedian(off_brain)
iqr(off_brain)
