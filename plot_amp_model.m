%% Porovnání výstupu reálného kytarového zesilovače a výstupu modelu
clear; clc; close all;

%% ===== SOUBORY =====
file_real  = 'Guitar real OUT M3 copy.wav';      % reálný zesilovač
file_model = 'M3 pro přílohu time var.wav';      % výstup modelu

% Soubor, do kterého bude uložen zarovnaný a gain-matched výstup modelu
outFile_model_gain_matched = 'model_gain_matched_aligned.wav';

%% ===== NAČTENÍ AUDIA =====
[x_real, fs_real]   = audioread(file_real);
[x_model, fs_model] = audioread(file_model);

% stereo -> mono
if size(x_real,2) > 1
    x_real = mean(x_real, 2);
end

if size(x_model,2) > 1
    x_model = mean(x_model, 2);
end

%% ===== PŘEVZORKOVÁNÍ V PŘÍPADĚ POTŘEBY =====
if fs_real ~= fs_model
    x_model = resample(x_model, fs_real, fs_model);
end
fs = fs_real;

%% ===== ODSTRANĚNÍ DC SLOŽKY =====
x_real  = x_real  - mean(x_real);
x_model = x_model - mean(x_model);

%% ===== VOLITELNÁ ŠPIČKOVÁ NORMALIZACE PŘED OŘEZEM =====
if max(abs(x_real)) > 0
    x_real = x_real / max(abs(x_real));
end

if max(abs(x_model)) > 0
    x_model = x_model / max(abs(x_model));
end

%% ===== OŘEZÁNÍ TICHA NA ZAČÁTKU A NA KONCI =====
% Práh ticha
trimThreshold = 0.01;              % lze nastavit např. 0.005 ... 0.02
minKeepSamples = round(0.02 * fs); % minimálně 20 ms aktivního signálu

x_real  = trim_silence_edges(x_real,  trimThreshold, minKeepSamples, fs);
x_model = trim_silence_edges(x_model, trimThreshold, minKeepSamples, fs);

%% ===== ZAROVNÁNÍ POMOCÍ KŘÍŽOVÉ KORELACE =====
% Pro signály s různou délkou se používá 'none', nikoliv 'coeff'
maxLagSec = 1.0;
maxLagSmp = round(maxLagSec * fs);

[c, lags] = xcorr(x_real, x_model, maxLagSmp, 'none');

% Ruční normalizace pro získání hodnoty podobné korelačnímu koeficientu
c = c / (norm(x_real) * norm(x_model) + eps);

[~, idxMax] = max(abs(c));
bestLag = lags(idxMax);

%% ===== ZAROVNÁNÍ SIGNÁLŮ =====
if bestLag > 0
    % model se zpožďuje -> ořízne se začátek reálného signálu
    x_real_aligned  = x_real(1+bestLag:end);
    x_model_aligned = x_model;
elseif bestLag < 0
    % model předbíhá -> ořízne se začátek modelovaného signálu
    shift = abs(bestLag);
    x_real_aligned  = x_real;
    x_model_aligned = x_model(1+shift:end);
else
    x_real_aligned  = x_real;
    x_model_aligned = x_model;
end

%% ===== SJEDNOCENÍ DÉLKY PO ZAROVNÁNÍ =====
N = min(length(x_real_aligned), length(x_model_aligned));
x_real_aligned  = x_real_aligned(1:N);
x_model_aligned = x_model_aligned(1:N);

%% ===== PŘIZPŮSOBENÍ ZISKU =====
% Přizpůsobení hlasitosti modelu reálnému zesilovači
g = (x_real_aligned' * x_model_aligned) / ...
    (x_model_aligned' * x_model_aligned + eps);

x_model_aligned = x_model_aligned * g;

%% ===== ULOŽENÍ GAIN-MATCHED MODELU =====
% Ukládá se stejná verze modelu, která se dále používá
% pro graf a výpočet metrik po zarovnání a přizpůsobení zisku.

x_model_save = x_model_aligned;

% Ochrana WAV souboru proti clippingu
% Pokud signál překročí rozsah [-1; 1], bude WAV soubor dodatečně normalizován.
% Tato operace nemění x_model_aligned pro graf a metriky, pouze uložený soubor.
peak_save = max(abs(x_model_save));

if peak_save > 1
    warning('Gain-matched model exceeds 1.0 FS. WAV will be normalized to avoid clipping.');
    x_model_save = x_model_save / peak_save;
end

audiowrite(outFile_model_gain_matched, x_model_save, fs);

% Dodatečné uložení RAW verze do .mat bez normalizace proti clippingu
save('model_gain_matched_aligned_raw.mat', ...
     'x_model_aligned', 'fs', 'g', 'bestLag');

%% ===== METRIKY =====
lag_ms = bestLag / fs * 1000;

R = corrcoef(x_real_aligned, x_model_aligned);
if numel(R) >= 4
    r = R(1,2);
else
    r = NaN;
end

rmse = sqrt(mean((x_real_aligned - x_model_aligned).^2));

%% ===== ČASOVÁ OSA =====
t = (0:N-1) / fs;

%% ===== SVISLÉ POSUNY =====
offset_real  = -0.9;
offset_model =  0.9;

y_real_plot  = x_real_aligned  + offset_real;
y_model_plot = x_model_aligned + offset_model;

%% ===== GRAF =====
figure('Color','w','Position',[100 100 1500 750]);
hold on;

plot(t, y_real_plot, ...
    'Color', [0.00 0.30 0.75], ...
    'LineWidth', 0.8);

plot(t, y_model_plot, ...
    'Color', [0.90 0.40 0.00], ...
    'LineWidth', 0.8);

yline(offset_real,  ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);
yline(offset_model, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);

grid on;
box on;

xlabel('Time (s)', 'FontSize', 15);
ylabel('Amplitude (FS)', 'FontSize', 15);

title(sprintf('lag %.0f ms  |  r %.3f  |  RMSE %.4f FS  |  gain %.3f', ...
    lag_ms, r, rmse, g), ...
    'FontSize', 16, 'FontWeight', 'normal');

legend({'Zesilovač', 'Model'}, ...
    'Location', 'southoutside', ...
    'Orientation', 'horizontal', ...
    'FontSize', 12);

set(gca, 'FontSize', 13);
xlim([0 t(end)]);
ylim([-2.2 2.2]);

%% ===== VÝPIS VÝSLEDKŮ =====
fprintf('Best lag      : %d samples\n', bestLag);
fprintf('Best lag      : %.3f ms\n', lag_ms);
fprintf('Correlation r : %.6f\n', r);
fprintf('RMSE          : %.6f FS\n', rmse);
fprintf('Gain applied to model: %.6f\n', g);
fprintf('Final length  : %d samples (%.3f s)\n', N, N/fs);
fprintf('Saved WAV     : %s\n', outFile_model_gain_matched);
fprintf('Saved MAT     : model_gain_matched_aligned_raw.mat\n');

%% ===== LOKÁLNÍ FUNKCE =====
function y = trim_silence_edges(x, thr, minLen, fs)
    x = x(:);
    mask = abs(x) > thr;

    if ~any(mask)
        y = x;
        return;
    end

    d = diff([false; mask; false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;

    segLens = ends - starts + 1;
    valid = segLens >= minLen;

    if ~any(valid)
        % Pokud nejsou nalezeny žádné dostatečně dlouhé úseky,
        % vezme se oblast mezi prvním a posledním překročením prahu
        firstIdx = find(mask, 1, 'first');
        lastIdx  = find(mask, 1, 'last');
        y = x(firstIdx:lastIdx);
        return;
    end

    starts = starts(valid);
    ends   = ends(valid);

    firstIdx = starts(1);
    lastIdx  = ends(end);

    % Malá rezerva na okrajích
    pad = round(0.01 * fs); % přibližně 10 ms
    firstIdx = max(1, firstIdx - pad);
    lastIdx  = min(length(x), lastIdx + pad);

    y = x(firstIdx:lastIdx);
end