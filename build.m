clear; clc;

%% Při jakych úrovních bylo změřeno
levelsDb = -19:-1:-40;   % -15, -16, ..., -40 dBFS

kernelFiles = cell(1, numel(levelsDb));
for i = 1:numel(levelsDb)
    kernelFiles{i} = sprintf('H3_%ddBFS.mat', levelsDb(i));
end

fRefHz  = 5000;    % referencni kmitocet
fs_meas = 44100;   % fs 
OSF     = 4;       % oversampling 

%% 1) Vytvoření modelu(bez kalibrace vystupu) na zakladu měření 
model = build_ph_level_model(levelsDb, kernelFiles, fRefHz, fs_meas);

%% 2) Kalibrace
model = calibrate_ph_gain_vs_level_multi(model, levelsDb, OSF);

%% 3) Uloženi 
save('FenderJunior_PH_interp_model_M3_10kernels more.mat', 'model');



