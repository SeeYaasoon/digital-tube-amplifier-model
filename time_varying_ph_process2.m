function [y, levelUsedVec, g_norm_vec, pk_before_norm_vec] = time_varying_ph_process2(x, model, varargin)
% TIME_VARYING_PH_PROCESS  Časově proměnné zpracování kytarového signálu pomocí paralelního Hammersteinova modelu.
%
%   [y, levelUsedVec, g_norm_vec, pk_before_norm_vec] = time_varying_ph_process(x, model, ...)
%
%   Parametry typu název-hodnota:
%     'FrameLenMs' - délka rámce v ms (výchozí 128)
%     'HopMs'      - posun rámce v ms (výchozí 32)
%     'LevelDb'    - pokud je zadána skalární hodnota, použije se pevná úroveň;
%                    jinak se úroveň odhaduje pro každý rámec zvlášť
%     'OSF'        - faktor převzorkování (výchozí 4)
%     'NormPeak'   - požadovaná špičková hodnota pro každý rámec nebo [] (výchozí [])
%     'MaxOrder'   - maximální použitý řád Hammersteinova modelu nebo [] pro všechna jádra

    p = inputParser;
    addParameter(p,'FrameLenMs',128); 
    addParameter(p,'HopMs',32);       
    addParameter(p,'LevelDb',[]);
    addParameter(p,'OSF',4);
    addParameter(p,'NormPeak',[]);
    addParameter(p,'MaxOrder',[], @(v) isempty(v) || (isscalar(v) && v >= 1));
    parse(p,varargin{:});

    FrameLenMs = p.Results.FrameLenMs;
    HopMs = p.Results.HopMs;
    levelDb_in = p.Results.LevelDb;
    OSF = p.Results.OSF;
    targetPeak = p.Results.NormPeak;
    maxOrder = p.Results.MaxOrder;

    x = double(x(:));
    fs = model.fs;
    ord = model.orders;
    N = numel(x);

    % rámce
    frameLen = round(FrameLenMs/1000 * fs);
    hop = round(HopMs/1000 * fs);

    if hop <= 0
        hop = max(1, floor(frameLen/2));
    end

    w = hann(frameLen,'periodic');
    nFrames = ceil((N - frameLen)/hop) + 1;

    y = zeros(N,1);
    winSum = zeros(N,1);

    levelUsedVec = zeros(nFrames,1);
    g_norm_vec = ones(nFrames,1);
    pk_before_norm_vec = zeros(nFrames,1);

    for k = 0:(nFrames-1)

        idx = (1:frameLen) + k*hop;
        idx_valid = idx(idx >= 1 & idx <= N);

        frame = zeros(frameLen,1);
        frame(idx_valid - k*hop) = x(idx_valid);

        % odhad úrovně pro tento rámec
        if isempty(levelDb_in)
            peak = max(abs(frame));
            levelDb = 20*log10(peak + eps);
        else
            levelDb = levelDb_in;
        end

        % získání jader pro danou úroveň
        [H_level, levelUsed] = ph_get_kernels_for_level(model, levelDb);

        [Lh, M_available] = size(H_level); %#ok<ASGLU>
        levelUsedVec(k+1) = levelUsed;

        % ------------------------------------------------------------
        % Omezení maximálního použitého řádu / počtu jader
        % ------------------------------------------------------------
        if isempty(maxOrder)
            M = M_available;
        else
            M = min(round(maxOrder), M_available);
        end

        % Bezpečnostní kontrola: nepoužívat více řádů, než existuje v model.orders
        M = min(M, numel(ord));

        % Ponechání pouze povolených jader a řádů
        H_level_used = H_level(:, 1:M);
        ord_used = ord(1:M);
        % ------------------------------------------------------------

        % převzorkování rámce
        [frame_os, ~] = os_zeropad_fft(frame, fs, OSF);

        % zpracování jednotlivých větví
        y_frame = zeros(frameLen,1);

        for m = 1:M
            z_os = frame_os.^ord_used(m);
            z = resample(z_os, 1, OSF); % anti-aliasing při převodu zpět na fs

            % kauzální konvoluce
            c = conv(z, H_level_used(:,m), 'full');
            convRes = c(1:length(z));

            y_frame = y_frame + convRes;
        end

        % kalibrační zesílení
        g_cal = 1;

        if isfield(model,'calibLevels') && isfield(model,'gain_vs_level')
            g_cal = interp1(model.calibLevels, model.gain_vs_level, levelUsed, 'linear', 'extrap');
        elseif isfield(model,'globalGain')
            g_cal = model.globalGain;
        end

        y_frame = g_cal * y_frame;

        % špičková hodnota před normalizací
        pk_before_norm = max(abs(y_frame));
        pk_before_norm_vec(k+1) = pk_before_norm;

        % volitelná normalizace každého rámce
        if ~isempty(targetPeak)
            if pk_before_norm > 0
                g_norm = targetPeak / pk_before_norm;
            else
                g_norm = 1;
            end

            y_frame = g_norm * y_frame;
            g_norm_vec(k+1) = g_norm;
        end

        % overlap-add s oknem
        win = w;

        local_idx = idx_valid - k*hop;

        y(idx_valid) = y(idx_valid) + y_frame(local_idx) .* win(local_idx);
        winSum(idx_valid) = winSum(idx_valid) + win(local_idx);
    end

    % normalizace součtem oken pro rekonstrukci signálu
    nonzero = winSum > 1e-12;
    y(nonzero) = y(nonzero) ./ winSum(nonzero);

    % ochrana proti hodnotám NaN / Inf
    y(~isfinite(y)) = 0;
end