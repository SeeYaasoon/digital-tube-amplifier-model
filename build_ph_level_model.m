function model = build_ph_level_model(levelsDb, kernelFiles, fRefHz, fs_default)
%BUILD_PH_LEVEL_MODEL Sestavení PH-modelu s interpolací podle úrovně.
%
% levelsDb    - vektor úrovní (v dBFS), při kterých se měřil ESS
% kernelFiles - cell pole cest k .mat souborům s jádry
% fRefHz      - referenční frekvence (např. 5000)
% fs_default  - vzorkovací frekvence, pokud není uvedena v souborech
%
% Předpokládá se, že v každém .mat souboru je:
%   HammersteinKernels  (Lh × M)
%   (volitelné) orders  (1×M) – řád pro každý sloupec

    assert(numel(levelsDb) == numel(kernelFiles), ...
        'levelsDb a kernelFiles musí mít stejnou délku');

    G = numel(levelsDb);

    %% 1) Načtení prvního souboru, zjištění rozměrů, fs a ord
    S = load(kernelFiles{1});

    % --- hledáme jádra ---
    Hfield = '';
    candH = {'HammersteinKernels','H','ir'};
    for k = 1:numel(candH)
        if isfield(S, candH{k})
            Hfield = candH{k};
            break;
        end
    end
    if isempty(Hfield)
        error('V souboru %s nebylo nalezeno pole s jádry.', kernelFiles{1});
    end
    H = S.(Hfield);            % Lh × M
    [Lh, M] = size(H);

    % --- hledáme ord (řády) ---
    if isfield(S,'orders') && numel(S.orders)==M
        ord = S.orders(:).';   % 1×M
    else
        ord = 1:M;
    end

    % --- hledáme fs ---
    fs = [];
    candFs = {'fs','Fs','fsamp','Fsamp'};
    for k = 1:numel(candFs)
        if isfield(S, candFs{k})
            fs = S.(candFs{k});
            break;
        end
    end
    if isempty(fs)
        if nargin >= 4 && ~isempty(fs_default)
            fs = fs_default;
        else
            error(['V souboru %s není žádné pole fs/Fs/fsamp/Fsamp. ' ...
                   'Předej fs_default jako čtvrtý argument.'], kernelFiles{1});
        end
    end

    %% 2) Normalizace první sady jader podle sweepu levelsDb(1)
    level1 = levelsDb(1);
    Ameas1 = 10^(level1/20);      % amplituda sweepu při této úrovni

    for m = 1:M
        H(:,m) = H(:,m) / (Ameas1.^ord(m));
    end

    H_all = zeros(Lh, M, G);
    H_all(:,:,1) = H;

    %% 3) Načtení a normalizace ostatních sad
    for g = 2:G
        Sg = load(kernelFiles{g});
        if ~isfield(Sg, Hfield)
            error('V souboru %s chybí pole %s.', kernelFiles{g}, Hfield);
        end
        Hg = Sg.(Hfield);
        if ~isequal(size(Hg), [Lh, M])
            error('Rozměr jader v souboru %s neodpovídá.', kernelFiles{g});
        end

        % pokud má tento soubor vlastní orders – zkontrolujeme shodu
        if isfield(Sg,'orders')
            if ~isequal(reshape(Sg.orders,1,[]), ord)
                error('orders v souboru %s se liší od prvního souboru.', kernelFiles{g});
            end
        end

        levelg = levelsDb(g);           % např. -15, -16, ...
        Ameasg = 10^(levelg/20);

        for m = 1:M
            Hg(:,m) = Hg(:,m) / (Ameasg.^ord(m));
        end

        H_all(:,:,g) = Hg;
    end

    %% 4) FFT a modul na referenční frekvenci fRefHz
    Nfft = 2^nextpow2(Lh * 2);
    if fRefHz >= fs/2
        error('fRefHz=%.1f Hz >= Nyquist=%.1f Hz. Zvol menší fRef.', ...
              fRefHz, fs/2);
    end
    kRef = round(fRefHz / fs * Nfft) + 1;

    mag_vs_level = zeros(M, G);   % [řád × index_úrovně]

    for g = 1:G
        Hf = fft(H_all(:,:,g), Nfft, 1);
        mag_vs_level(:, g) = abs(Hf(kRef, :)).';
    end

    %% 5) Referenční úroveň (střed rozsahu)
    midLevel = mean([min(levelsDb), max(levelsDb)]);
    [~, idxRef] = min(abs(levelsDb - midLevel));

    H_ref   = H_all(:,:,idxRef);
    mag_ref = mag_vs_level(:, idxRef);

    %% 6) Sestavení modelu
    model.fs            = fs;
    model.orders        = ord;
    model.Lh            = Lh;
    model.fRefHz        = fRefHz;

    model.levelsDb      = levelsDb(:).';
    model.kernelFiles   = kernelFiles(:).';

    model.H_ref         = H_ref;
    model.mag_ref       = mag_ref;
    model.mag_vs_level  = mag_vs_level;
    model.H_all         = H_all;   % všechna jádra pro všechny úrovně (časová oblast)
end





