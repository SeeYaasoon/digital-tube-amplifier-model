function model_out = calibrate_ph_gain_vs_level_multi(model, levelsDb, OSF)
%CALIBRATE_PH_GAIN_VS_LEVEL_MULTI  Kalibrace gain(level) podle všech sweepů.
%
%  model_out = calibrate_ph_gain_vs_level_multi(model, levelsDb, OSF)
%
%  levelsDb - vektor úrovní, pro které existují dvojice sweep/odezva, např. -15:-1:-40
%  OSF      - oversampling factor , např. 4.
%
%  Očekávané názvy souborů:
%    inFile  = sprintf('Sweep %ddBFS.wav',  levelDb);
%    outFile = sprintf('MIC M6 %ddBFS VER2.wav', levelDb);

    if nargin < 3 || isempty(OSF), OSF = 4; end

    % před kalibrací odstraníme stará kalibrační pole, pokud existují
    if isfield(model,'calibLevels'),   model = rmfield(model,'calibLevels');   end
    if isfield(model,'gain_vs_level'), model = rmfield(model,'gain_vs_level'); end
    if isfield(model,'globalGain'),    model.globalGain = 1;                   end

    K = numel(levelsDb);
    gains      = zeros(1, K);
    levelsUsed = zeros(1, K);

    for k = 1:K
        L = levelsDb(k);

        % === názvy souborů pro tuto úroveň ===
        inFile  = sprintf('Sweep %ddBFS.wav',  L);
        outFile = sprintf('Response M3 %ddBFS.wav', L);
        fprintf('Calib for level %d dBFS: in="%s", out="%s"\n', L, inFile, outFile);

        % === načtení vstupu a změřeného výstupu ===
        [x, fs1] = audioread(inFile);
        if size(x,2) > 1, x = mean(x,2); end

        [y_meas, fs2] = audioread(outFile);
        if size(y_meas,2) > 1, y_meas = mean(y_meas,2); end

        assert(fs1 == model.fs && fs2 == model.fs, 'fs se neshoduje s model.fs');

        N = min(numel(x), numel(y_meas));
        x      = x(1:N);
        y_meas = y_meas(1:N);

        % === modelový výstup pro tuto úroveň (bez kalibrace) ===
        [y_mod, levelUsed] = ph_process_interp( ...
            x, model, 'LevelDb', L, 'OSF', OSF, 'NormPeak', []);

        y_mod = y_mod(1:N);

        % === RMS gain ===
        rms_meas = sqrt(mean(y_meas.^2));
        rms_mod  = sqrt(mean(y_mod.^2)) + eps;

        gains(k)      = rms_meas / rms_mod;
        levelsUsed(k) = levelUsed;

        fprintf('  used level = %.2f dBFS, gain = %.4f (%.2f dB)\n', ...
            levelUsed, gains(k), 20*log10(gains(k)+eps));
    end

    model_out = model;
    model_out.calibLevels   = levelsUsed;   % může být i jen levelsDb
    model_out.gain_vs_level = gains;
end





