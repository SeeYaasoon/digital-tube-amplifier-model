function [H_level, levelUsed, baseIdx, isExtrap] = ph_get_kernels_for_level(model, levelDb, varargin)
%PH_GET_KERNELS_FOR_LEVEL  Vrací H_level (Lh×M) pro zadanou úroveň levelDb.
%
% [H_level, levelUsed, baseIdx, isExtrap] = ph_get_kernels_for_level(model, levelDb, ...)
% Parametry typu název-hodnota:
%   'InterpMethod' - 'pchip' (výchozí), 'makima', 'spline', 'linear'
%   'AllowExtrap'  - true/false (výchozí false) — povolit extrapolaci úrovně
%
% Model musí obsahovat následující pole:
%   model.levelsDb  (1×G)
%   model.H_all     (Lh×M×G)

    p = inputParser;
    addParameter(p,'InterpMethod','pchip');
    addParameter(p,'AllowExtrap',false);
    parse(p,varargin{:});
    method = validatestring(p.Results.InterpMethod, {'pchip','makima','spline','linear'});
    allowExtrap = p.Results.AllowExtrap;

    % Kontrola, zda jsou přítomna požadovaná pole
    assert(isfield(model,'levelsDb') && isfield(model,'H_all'), 'model: missing required fields levelsDb or H_all');
    levels = double(model.levelsDb(:)).';    % 1×G
    H_all = model.H_all;                     % Lh×M×G
    [Lh, M, G] = size(H_all);
    assert(numel(levels) == G, 'model.levelsDb length must match third dim of H_all');

    % Omezení nebo použití levelUsed podle hodnoty allowExtrap
    if allowExtrap
        levelUsed = double(levelDb);
        isExtrap = (levelUsed < min(levels)) || (levelUsed > max(levels));
    else
        levelUsed = min(max(levelDb, min(levels)), max(levels));
        isExtrap = false;
    end

    % Nalezení nejbližšího měřeného bodu pro baseIdx
    [~, baseIdx] = min(abs(levels - levelUsed));

    % Pokud existuje přesná shoda, použijí se přímo uložená jádra
    if abs(levelUsed - levels(baseIdx)) < 1e-9
        H_level = double(H_all(:,:,baseIdx));
        return;
    end

    % Interpolace celého tenzoru H_all podél osy úrovně
    % změna tvaru matice: G × (Lh*M)
    Hmat = reshape(H_all, Lh*M, G).';  % G × (Lh*M)

    % Volba chování extrapolace pro interp1
    if allowExtrap
        Hinterp = interp1(levels, Hmat, levelUsed, method, 'extrap');  % 1 × (Lh*M)
    else
        Hinterp = interp1(levels, Hmat, levelUsed, method);            % NaN, pokud je hodnota mimo rozsah
    end

    % Pokud se objeví hodnoty NaN, výsledek se nahradí nejbližšími měřenými jádry
    if any(isnan(Hinterp))
        % Záložní řešení: použijí se jádra z nejbližšího měřeného bodu
        warning('ph_get_kernels_for_level:interpNaN', ...
            'Interpolation produced NaN (out of range). Using nearest measured kernels.');
        H_level = double(H_all(:,:,baseIdx));
        return;
    end

    H_level = reshape(Hinterp.', Lh, M);  % Lh × M

    % Dodatečná kontrola: ověření, že se nevyskytují hodnoty NaN nebo Inf
    if any(~isfinite(H_level(:)))
        warning('ph_get_kernels_for_level:nonfinite', 'Resulting H_level contains non-finite values; using nearest measured kernels.');
        H_level = double(H_all(:,:,baseIdx));
    end
end




