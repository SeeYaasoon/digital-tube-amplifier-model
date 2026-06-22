function [y, fs_os] = os_zeropad_fft(x, fs, L)
% Pásmově omezené převzorkování faktorem L pomocí doplnění nul ve spektru.
% x — vektor nebo matice (N×C, sloupce = kanály/jádra)
% Vrací y (N*L × C) a novou vzorkovací frekvenci Fs.

    if nargin<3 || isempty(L), L = 2; end
    x = double(x);
    N = size(x,1); C = size(x,2);

    % FFT podél času, tedy podél první dimenze
    X = fft(x, [], 1);                    % N × C

    % Vytvoření spektra délky N*L: [0..Nyq] + nuly + (Nyq+1..N-1)
    if mod(N,2)==0                        % sudá délka, existuje explicitní Nyquistova frekvence
        kNy = N/2 + 1;                    % index Nyquistovy frekvence
        top = X(1:kNy, :);                % DC..Nyq
        bot = X(kNy+1:end, :);            % zbývající část záporných frekvencí
    else                                  % lichá délka, neexistuje explicitní Nyquistova frekvence
        kNy = (N+1)/2;
        top = X(1:kNy, :);
        bot = X(kNy+1:end, :);
    end

    pad = zeros((L-1)*N, C);
    Xup = [top; pad; bot];                % délka = N*L

    % Inverzní FFT a korekce amplitudy (*L)
    y = real(ifft(Xup, [], 1)) * L;       % N*L × C
    fs_os = fs * L;
end


