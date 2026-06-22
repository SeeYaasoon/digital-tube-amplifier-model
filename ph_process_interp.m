function [y, levelUsed, g_norm, pk_before_norm] = ph_process_interp(x, model, varargin)
%PH_PROCESS_INTERP  Zpracování signálu PH-modelem s interpolací podle úrovně.
%
% [y, levelUsed, g_norm, pk_before_norm] = ph_process_interp(x, model, ...)
%
% Name-Value parametry:
%   'LevelDb'  - pevná úroveň (dBFS). Pokud je prázdné -> odhad podle peaku.
%   'OSF'      - oversampling factor (výchozí 4).
%   'NormPeak' - požadovaný peak na výstupu (0..1). Pokud je [] nebo není zadán,
%                neprovádí se normalizace podle peaku.
%
% g_norm         - koeficient normalizace podle peaku (pokud byla použita).
% pk_before_norm - peak signálu PŘED NormPeak, ale už po kalibračním gainu.

    p = inputParser;
    addParameter(p, 'LevelDb', []);
    addParameter(p, 'OSF', 4);
    addParameter(p, 'NormPeak', []);
    parse(p, varargin{:});
    levelDb_in = p.Results.LevelDb;
    OSF        = p.Results.OSF;
    targetPeak = p.Results.NormPeak;

    x = double(x(:));
    fs = model.fs;
    ord = model.orders;

    % --- odhad úrovně, pokud není zadána explicitně ---
    if isempty(levelDb_in)
        peak = max(abs(x));
        levelDb_in = 20*log10(peak + eps);
    end

    % --- jádra pro tuto úroveň ---
    [H_level, levelUsed] = ph_get_kernels_for_level(model, levelDb_in);
    [Lh, M] = size(H_level);

    % --- oversampling + nelinearita + LTI ---
    [x_os, ~] = os_zeropad_fft(x, fs, OSF);

    y = zeros(size(x));
    for m = 1:M
        z_os = x_os.^ord(m);        % nelinearita na OSF*fs
        z    = resample(z_os, 1, OSF); % anti-alias -> fs
        y    = y + conv(z, H_level(:,m), 'same');
    end

    % --- kalibrace gain(level), pokud je k dispozici ---
    g_cal = 1;
    if isfield(model,'calibLevels') && isfield(model,'gain_vs_level')
        g_cal = interp1(model.calibLevels, model.gain_vs_level, ...
                        levelUsed, 'linear', 'extrap');
        y = g_cal * y;
    elseif isfield(model,'globalGain')
        % fallback, pokud by náhodou zůstal starý globalGain
        y = model.globalGain * y;
    end

    % --- peak před normalizací ---
    pk_before_norm = max(abs(y));

    % --- volitelná NormPeak pro poslech ---
    g_norm = 1;
    if ~isempty(targetPeak)
        [y, g_norm, ~] = norm_peak(y, targetPeak);
    end
end








