clear; clc; clear sound;

%% 1) Model
load('FenderJunior_PH_interp_model_M3.mat','model');

%% 2) Vstupní signál
[x, fs] = audioread('Guitar LOOPBACK.wav');
x = x * 1;

if size(x,2) > 1
    x = mean(x,2);
end

assert(fs == model.fs, 'fs nesouhlasi s fs modelu');

%% 3) Zpracování signálu

% První možnost: časově proměnné zpracování.
% Funkce time_varying_ph_process rozděluje vstupní signál na krátké rámce.
% Pro každý rámec zvlášť vypočítá jeho špičkovou hodnotu a podle ní vybere
% nebo interpoluje odpovídající sadu jader Hammersteinového modelu.
% Díky tomu se model může měnit v čase podle aktuální úrovně signálu.
%
[y, levelUsed, g_out, pk_before] = time_varying_ph_process( ...
    x, model, 'OSF', 1, 'NormPeak', []);

 %[y, levelUsed, g_out, pk_before] = time_varying_ph_process2( ...
 %   x, model, 'OSF', 1, 'NormPeak', [], 'MaxOrder', 10)

% Druhá možnost: statické zpracování.
% Funkce ph_process_interp určí úroveň pouze jednou, a to ze špičkové
% hodnoty celého vstupního souboru.
% Podle této jedné hodnoty vybere nebo interpoluje jednu sadu jader,
% která se potom použije pro zpracování celého signálu.
%[y, levelUsed, g_out, pk_before] = ph_process_interp( ...
 %  x, model, 'OSF', 4, 'NormPeak', []);

%% 4) Výpis hodnot
fprintf('Vstupní špička        = %.2f dBFS\n', 20*log10(max(abs(x))+eps));
fprintf('Použita úroveň jader  = %.2f dBFS\n', levelUsed);
fprintf('Výstupní špička       = %.2f dBFS\n', 20*log10(pk_before+eps));

%% 5) Přehrání a uložení výstupu
sound(y, fs);
audiowrite('audio.wav', y, fs);


