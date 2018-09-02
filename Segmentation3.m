clear; 

Fs = 256; % sampling rate
Fst = Fs*60; % total samples

%%%%%%%%%%%%% Create Datasets %%%%%%%%%%%%%
%%%%% Cube %%%%%
% Search directory (normal datasets)
names = dir('datasets\cube'); % Input datasets
names = {names.name};
names = names(3:length(names));

% Pre-allocation for speed
counter = 1;
tablegen = 1;

% Define the wavelet family and the level for decomposition
waveletFunction='db4';
level=5;

% Load normal datasets, partitioning and filtering
for i = 1:length(names)
    raw = load(['datasets\cube\' names{i}]);
    raw = raw(1:Fst, :); % Take datasets from 1-60 seconds
    
    % Filtering
    % Bandpass Filter
    Fn = Fs/2;                         % Nyquist Frequency (Hz)
    Wp = [0.51 64]/Fn;                 % Passband Frequencies (Normalised)
    Ws = [0.5 64.01]/Fn;               % Stopband Frequencies (Normalised)
    Rp = 5;                            % Passband Ripple (dB)
    Rs = 30;                           % Stopband Ripple (dB)
    [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);    % Filter Order
    [z,p,k] = cheby2(n,Rs,Ws);         % Filter Design
    [sosbp,gbp] = zp2sos(z,p,k);       % Convert To Second-Order-Section For Stability
    raw = filtfilt(sosbp,gbp, raw);    % Filter Signal
    
    % Partition data
    EEG{1} = raw(1:20*Fs, :);         % Take datasets from 0-20 seconds
    EEG{2} = raw((20*Fs)+1:40*Fs, :); % Take datasets from 21-40 seconds
    EEG{3} = raw((40*Fs)+1:Fst, :);   % Take datasets from 41-60 seconds
    
    % EEG Sub-band Extraction & Feature Extraction
    for ii = 1:3
        % Discrete Fourier Transform
        EEGfft = fft(EEG{1, ii});
        EEGfftshift = fftshift(EEGfft);
        EEGfftshift = abs(EEGfftshift(length(EEGfftshift)/2:length(EEGfftshift), :));
        
        % Multilevel 1-D wavelet decomposition for both channels
        [C,L]=wavedec(EEG{1, ii}(:, 1),level,waveletFunction);
        [D,L]=wavedec(EEG{1, ii}(:, 2),level,waveletFunction);
        
        % Decompose into EEG Sub-bands
        for iii = 1:2
            if iii == 1
            % Reconstruct single branch from 1-D wavelet coefficients for ch0
                Gamma(:, iii) = wrcoef('d',C,L,waveletFunction,2); %GAMMA, 32 - 64
                Beta(:, iii) = wrcoef('d',C,L,waveletFunction,3); %BETA, 16 - 32
                Alpha(:, iii) = wrcoef('d',C,L,waveletFunction,4); %ALPHA, 8 - 16
                Theta(:, iii) = wrcoef('d',C,L,waveletFunction,5); %THETA, 4 - 8
                Delta(:, iii) = wrcoef('a',C,L,waveletFunction,5); %DELTA, 0 - 4
            elseif iii == 2
            % Reconstruct single branch from 1-D wavelet coefficients for ch1
                Gamma(:, iii) = wrcoef('d',D,L,waveletFunction,2); %GAMMA, 32 - 64
                Beta(:, iii) = wrcoef('d',D,L,waveletFunction,3); %BETA, 16 - 32
                Alpha(:, iii) = wrcoef('d',D,L,waveletFunction,4); %ALPHA, 8 - 16
                Theta(:, iii) = wrcoef('d',D,L,waveletFunction,5); %THETA, 4 - 8
                Delta(:, iii) = wrcoef('a',D,L,waveletFunction,5); %DELTA, 0 - 4
            end
        end
        
        % FFT
        % Perform FFT on each band, then extract features from the FFTs
        Gamma_FFT = fft(Gamma);
        Gammafftshift = fftshift(Gamma_FFT);
        Gammafftshift = abs(Gammafftshift(length(Gammafftshift)/2:length(Gammafftshift), :));
        
        Beta_FFT = fft(Beta);
        Betafftshift = fftshift(Beta_FFT);
        Betafftshift = abs(Betafftshift(length(Betafftshift)/2:length(Betafftshift), :));
        
        Alpha_FFT = fft(Alpha);
        Alphafftshift = fftshift(Alpha_FFT);
        Alphafftshift = abs(Alphafftshift(length(Alphafftshift)/2:length(Alphafftshift), :));
        
        Theta_FFT = fft(Theta);
        Thetafftshift = fftshift(Theta_FFT);
        Thetafftshift = abs(Thetafftshift(length(Thetafftshift)/2:length(Thetafftshift), :));
        
        Delta_FFT = fft(Delta);
        Deltafftshift = fftshift(Delta_FFT);
        Deltafftshift = abs(Deltafftshift(length(Deltafftshift)/2:length(Deltafftshift), :));
        
        % Labelling Classes
        counter = counter + 1;
        action(tablegen, 1) = categorical(cellstr('Cube'));
        
        % Spectral Power
        power_gamma(tablegen, :) = sum(Gammafftshift);
        power_beta(tablegen, :) = sum(Betafftshift);
        power_alpha(tablegen, :) = sum(Alphafftshift);
        power_theta(tablegen, :) = sum(Thetafftshift);
        power_delta(tablegen, :) = sum(Deltafftshift);
        powerdiff_gamma(tablegen, :) = abs(power_gamma(tablegen, 1) - power_gamma(tablegen, 2));
        powerdiff_beta(tablegen, :) = abs(power_beta(tablegen, 1) - power_beta(tablegen, 2));
        powerdiff_alpha(tablegen, :) = abs(power_alpha(tablegen, 1) - power_alpha(tablegen, 2));
        powerdiff_theta(tablegen, :) = abs(power_theta(tablegen, 1) - power_theta(tablegen, 2));
        powerdiff_delta(tablegen, :) = abs(power_delta(tablegen, 1) - power_delta(tablegen, 2));

        % DASDV
        for row = 2:length(EEG{1, 1})
            for iii = 1:2 % Both Channels
                rowcalc_gamma(row, iii) = Gamma(row, iii) - Gamma(row-1, iii);
                rowcalc_beta(row, iii) = Beta(row, iii) - Beta(row-1, iii);
                rowcalc_alpha(row, iii) = Alpha(row, iii) - Alpha(row-1, iii);
                rowcalc_theta(row, iii) = Theta(row, iii) - Theta(row-1, iii);
                rowcalc_delta(row, iii) = Delta(row, iii) - Delta(row-1, iii);
            end
        end
        dasdv_gamma(tablegen, :) = sqrt(abs(sum(rowcalc_beta.^2)/Fst-1));
        dasdv_beta(tablegen, :) = sqrt(abs(sum(rowcalc_beta.^2)/Fst-1));
        dasdv_alpha(tablegen, :) = sqrt(abs(sum(rowcalc_alpha.^2)/Fst-1));
        dasdv_theta(tablegen, :) = sqrt(abs(sum(rowcalc_theta.^2)/Fst-1));
        dasdv_delta(tablegen, :) = sqrt(abs(sum(rowcalc_delta.^2)/Fst-1));
        
        % Other Simple Features
        meanamp(tablegen, :) = mean(EEGfftshift); % Mean
        maxamp(tablegen, :) = max(EEGfftshift); % Max
        varamp(tablegen, :) = var(EEGfftshift); % Variance
        
        % Next row
        tablegen = tablegen + 1;
    end
end

%%%%% Math %%%%%
% Search directory (normal datasets)
names = dir('datasets\math'); % Input datasets
names = {names.name};
names = names(3:length(names));

% Pre-allocation for speed
counter = 1;

% Define the wavelet family and the level for decomposition
waveletFunction='db4';
level=5;

% Load normal datasets, partitioning and filtering
for i = 1:length(names)
    raw = load(['datasets\math\' names{i}]);
    raw = raw(1:Fst, :); % Take datasets from 1-60 seconds
    
    % Filtering
    % Bandpass Filter
    Fn = Fs/2;                         % Nyquist Frequency (Hz)
    Wp = [0.51 64]/Fn;                 % Passband Frequencies (Normalised)
    Ws = [0.5 64.01]/Fn;               % Stopband Frequencies (Normalised)
    Rp = 5;                            % Passband Ripple (dB)
    Rs = 30;                           % Stopband Ripple (dB)
    [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);    % Filter Order
    [z,p,k] = cheby2(n,Rs,Ws);         % Filter Design
    [sosbp,gbp] = zp2sos(z,p,k);       % Convert To Second-Order-Section For Stability
    raw = filtfilt(sosbp,gbp, raw);    % Filter Signal
    
    % Partition data
    EEG{1} = raw(1:20*Fs, :); % Take datasets from 0-20 seconds
    EEG{2} = raw((20*Fs)+1:40*Fs, :); % Take datasets from 21-40 seconds
    EEG{3} = raw((40*Fs)+1:Fst, :); % Take datasets from 41-60 seconds
    
    % EEG Sub-band Extraction & Feature Extraction
    for ii = 1:3
        % Discrete Fourier Transform
        EEGfft = fft(EEG{1, ii});
        EEGfftshift = fftshift(EEGfft);
        EEGfftshift = abs(EEGfftshift(length(EEGfftshift)/2:length(EEGfftshift), :));
        
        % Multilevel 1-D wavelet decomposition for both channels
        [C,L]=wavedec(EEG{1, ii}(:, 1),level,waveletFunction);
        [D,L]=wavedec(EEG{1, ii}(:, 2),level,waveletFunction);
        
        % Decompose into EEG Sub-bands
        for iii = 1:2
            if iii == 1
            % Reconstruct single branch from 1-D wavelet coefficients for ch0
                Gamma(:, iii) = wrcoef('d',C,L,waveletFunction,2); %GAMMA, 32 - 64
                Beta(:, iii) = wrcoef('d',C,L,waveletFunction,3); %BETA, 16 - 32
                Alpha(:, iii) = wrcoef('d',C,L,waveletFunction,4); %ALPHA, 8 - 16
                Theta(:, iii) = wrcoef('d',C,L,waveletFunction,5); %THETA, 4 - 8
                Delta(:, iii) = wrcoef('a',C,L,waveletFunction,5); %DELTA, 0 - 4
            elseif iii == 2
            % Reconstruct single branch from 1-D wavelet coefficients for ch1
                Gamma(:, iii) = wrcoef('d',D,L,waveletFunction,2); %GAMMA, 32 - 64
                Beta(:, iii) = wrcoef('d',D,L,waveletFunction,3); %BETA, 16 - 32
                Alpha(:, iii) = wrcoef('d',D,L,waveletFunction,4); %ALPHA, 8 - 16
                Theta(:, iii) = wrcoef('d',D,L,waveletFunction,5); %THETA, 4 - 8
                Delta(:, iii) = wrcoef('a',D,L,waveletFunction,5); %DELTA, 0 - 4
            end
        end
        
        % FFT
        % Perform FFT on each band, then extract features from the FFTs
        Gamma_FFT = fft(Gamma);
        Gammafftshift = fftshift(Gamma_FFT);
        Gammafftshift = abs(Gammafftshift(length(Gammafftshift)/2:length(Gammafftshift), :));
        
        Beta_FFT = fft(Beta);
        Betafftshift = fftshift(Beta_FFT);
        Betafftshift = abs(Betafftshift(length(Betafftshift)/2:length(Betafftshift), :));
        
        Alpha_FFT = fft(Alpha);
        Alphafftshift = fftshift(Alpha_FFT);
        Alphafftshift = abs(Alphafftshift(length(Alphafftshift)/2:length(Alphafftshift), :));
        
        Theta_FFT = fft(Theta);
        Thetafftshift = fftshift(Theta_FFT);
        Thetafftshift = abs(Thetafftshift(length(Thetafftshift)/2:length(Thetafftshift), :));
        
        Delta_FFT = fft(Delta);
        Deltafftshift = fftshift(Delta_FFT);
        Deltafftshift = abs(Deltafftshift(length(Deltafftshift)/2:length(Deltafftshift), :));
        
        % Labelling Classes
        counter = counter + 1;
        action(tablegen, 1) = categorical(cellstr('Math'));
        
        % Spectral Power
        power_gamma(tablegen, :) = sum(Gammafftshift);
        power_beta(tablegen, :) = sum(Betafftshift);
        power_alpha(tablegen, :) = sum(Alphafftshift);
        power_theta(tablegen, :) = sum(Thetafftshift);
        power_delta(tablegen, :) = sum(Deltafftshift);
        powerdiff_gamma(tablegen, :) = abs(power_gamma(tablegen, 1) - power_gamma(tablegen, 2));
        powerdiff_beta(tablegen, :) = abs(power_beta(tablegen, 1) - power_beta(tablegen, 2));
        powerdiff_alpha(tablegen, :) = abs(power_alpha(tablegen, 1) - power_alpha(tablegen, 2));
        powerdiff_theta(tablegen, :) = abs(power_theta(tablegen, 1) - power_theta(tablegen, 2));
        powerdiff_delta(tablegen, :) = abs(power_delta(tablegen, 1) - power_delta(tablegen, 2));

        % DASDV
        for row = 2:length(EEG{1, 1})
            for iii = 1:2 % Both Channels
                rowcalc_gamma(row, iii) = Gamma(row, iii) - Gamma(row-1, iii);
                rowcalc_beta(row, iii) = Beta(row, iii) - Beta(row-1, iii);
                rowcalc_alpha(row, iii) = Alpha(row, iii) - Alpha(row-1, iii);
                rowcalc_theta(row, iii) = Theta(row, iii) - Theta(row-1, iii);
                rowcalc_delta(row, iii) = Delta(row, iii) - Delta(row-1, iii);
            end
        end
        dasdv_gamma(tablegen, :) = sqrt(abs(sum(rowcalc_beta.^2)/Fst-1));
        dasdv_beta(tablegen, :) = sqrt(abs(sum(rowcalc_beta.^2)/Fst-1));
        dasdv_alpha(tablegen, :) = sqrt(abs(sum(rowcalc_alpha.^2)/Fst-1));
        dasdv_theta(tablegen, :) = sqrt(abs(sum(rowcalc_theta.^2)/Fst-1));
        dasdv_delta(tablegen, :) = sqrt(abs(sum(rowcalc_delta.^2)/Fst-1));
        
        % Other Simple Features
        meanamp(tablegen, :) = mean(EEGfftshift); % Mean
        maxamp(tablegen, :) = max(EEGfftshift); % Max
        varamp(tablegen, :) = var(EEGfftshift); % Variance
        
        % Next row
        tablegen = tablegen + 1;
    end
end

%%%%% Letter %%%%%
% Search directory (normal datasets)
names = dir('datasets\letter'); % Input datasets
names = {names.name};
names = names(3:length(names));

% Pre-allocation for speed
counter = 1;

% Define the wavelet family and the level for decomposition
waveletFunction='db4';
level=5;

% Load normal datasets, partitioning and filtering
for i = 1:length(names)
    raw = load(['datasets\letter\' names{i}]);
    raw = raw(1:Fst, :); % Take datasets from 1-60 seconds
    
    % Filtering
    % Bandpass Filter
    Fn = Fs/2;                         % Nyquist Frequency (Hz)
    Wp = [0.51 64]/Fn;                 % Passband Frequencies (Normalised)
    Ws = [0.5 64.01]/Fn;               % Stopband Frequencies (Normalised)
    Rp = 5;                            % Passband Ripple (dB)
    Rs = 30;                           % Stopband Ripple (dB)
    [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);    % Filter Order
    [z,p,k] = cheby2(n,Rs,Ws);         % Filter Design
    [sosbp,gbp] = zp2sos(z,p,k);       % Convert To Second-Order-Section For Stability
    raw = filtfilt(sosbp,gbp, raw);    % Filter Signal
    
    % Partition data
    EEG{1} = raw(1:20*Fs, :); % Take datasets from 0-20 seconds
    EEG{2} = raw((20*Fs)+1:40*Fs, :); % Take datasets from 21-40 seconds
    EEG{3} = raw((40*Fs)+1:Fst, :); % Take datasets from 41-60 seconds
    
    % EEG Sub-band Extraction & Feature Extraction
    for ii = 1:3
        % Discrete Fourier Transform
        EEGfft = fft(EEG{1, ii});
        EEGfftshift = fftshift(EEGfft);
        EEGfftshift = abs(EEGfftshift(length(EEGfftshift)/2:length(EEGfftshift), :));
        
        % Multilevel 1-D wavelet decomposition for both channels
        [C,L]=wavedec(EEG{1, ii}(:, 1),level,waveletFunction);
        [D,L]=wavedec(EEG{1, ii}(:, 2),level,waveletFunction);
        
        % Decompose into EEG Sub-bands
        for iii = 1:2
            if iii == 1
            % Reconstruct single branch from 1-D wavelet coefficients for ch0
                Gamma(:, iii) = wrcoef('d',C,L,waveletFunction,2); %GAMMA, 32 - 64
                Beta(:, iii) = wrcoef('d',C,L,waveletFunction,3); %BETA, 16 - 32
                Alpha(:, iii) = wrcoef('d',C,L,waveletFunction,4); %ALPHA, 8 - 16
                Theta(:, iii) = wrcoef('d',C,L,waveletFunction,5); %THETA, 4 - 8
                Delta(:, iii) = wrcoef('a',C,L,waveletFunction,5); %DELTA, 0 - 4
            elseif iii == 2
            % Reconstruct single branch from 1-D wavelet coefficients for ch1
                Gamma(:, iii) = wrcoef('d',D,L,waveletFunction,2); %GAMMA, 32 - 64
                Beta(:, iii) = wrcoef('d',D,L,waveletFunction,3); %BETA, 16 - 32
                Alpha(:, iii) = wrcoef('d',D,L,waveletFunction,4); %ALPHA, 8 - 16
                Theta(:, iii) = wrcoef('d',D,L,waveletFunction,5); %THETA, 4 - 8
                Delta(:, iii) = wrcoef('a',D,L,waveletFunction,5); %DELTA, 0 - 4
            end
        end
        
        % FFT
        % Perform FFT on each band, then extract features from the FFTs
        Gamma_FFT = fft(Gamma);
        Gammafftshift = fftshift(Gamma_FFT);
        Gammafftshift = abs(Gammafftshift(length(Gammafftshift)/2:length(Gammafftshift), :));
        
        Beta_FFT = fft(Beta);
        Betafftshift = fftshift(Beta_FFT);
        Betafftshift = abs(Betafftshift(length(Betafftshift)/2:length(Betafftshift), :));
        
        Alpha_FFT = fft(Alpha);
        Alphafftshift = fftshift(Alpha_FFT);
        Alphafftshift = abs(Alphafftshift(length(Alphafftshift)/2:length(Alphafftshift), :));
        
        Theta_FFT = fft(Theta);
        Thetafftshift = fftshift(Theta_FFT);
        Thetafftshift = abs(Thetafftshift(length(Thetafftshift)/2:length(Thetafftshift), :));
        
        Delta_FFT = fft(Delta);
        Deltafftshift = fftshift(Delta_FFT);
        Deltafftshift = abs(Deltafftshift(length(Deltafftshift)/2:length(Deltafftshift), :));
        
        % Labelling Classes
        counter = counter + 1;
        action(tablegen, 1) = categorical(cellstr('Letter'));
        
        % Spectral Power
        power_gamma(tablegen, :) = sum(Gammafftshift);
        power_beta(tablegen, :) = sum(Betafftshift);
        power_alpha(tablegen, :) = sum(Alphafftshift);
        power_theta(tablegen, :) = sum(Thetafftshift);
        power_delta(tablegen, :) = sum(Deltafftshift);
        
        % Spectral Power Difference
        powerdiff_gamma(tablegen, :) = abs(power_gamma(tablegen, 1) - power_gamma(tablegen, 2));
        powerdiff_beta(tablegen, :) = abs(power_beta(tablegen, 1) - power_beta(tablegen, 2));
        powerdiff_alpha(tablegen, :) = abs(power_alpha(tablegen, 1) - power_alpha(tablegen, 2));
        powerdiff_theta(tablegen, :) = abs(power_theta(tablegen, 1) - power_theta(tablegen, 2));
        powerdiff_delta(tablegen, :) = abs(power_delta(tablegen, 1) - power_delta(tablegen, 2));

        % DASDV
        for row = 2:length(EEG{1, 1})
            for iii = 1:2 % Both Channels
                rowcalc_gamma(row, iii) = Gamma(row, iii) - Gamma(row-1, iii);
                rowcalc_beta(row, iii) = Beta(row, iii) - Beta(row-1, iii);
                rowcalc_alpha(row, iii) = Alpha(row, iii) - Alpha(row-1, iii);
                rowcalc_theta(row, iii) = Theta(row, iii) - Theta(row-1, iii);
                rowcalc_delta(row, iii) = Delta(row, iii) - Delta(row-1, iii);
            end
        end
        dasdv_gamma(tablegen, :) = sqrt(abs(sum(rowcalc_beta.^2)/Fst-1));
        dasdv_beta(tablegen, :) = sqrt(abs(sum(rowcalc_beta.^2)/Fst-1));
        dasdv_alpha(tablegen, :) = sqrt(abs(sum(rowcalc_alpha.^2)/Fst-1));
        dasdv_theta(tablegen, :) = sqrt(abs(sum(rowcalc_theta.^2)/Fst-1));
        dasdv_delta(tablegen, :) = sqrt(abs(sum(rowcalc_delta.^2)/Fst-1));
        
        % Other Simple Features
        meanamp(tablegen, :) = mean(EEGfftshift); % Mean
        maxamp(tablegen, :) = max(EEGfftshift); % Max
        varamp(tablegen, :) = var(EEGfftshift); % Variance
        
        % Next row
        tablegen = tablegen + 1;
    end
end

%%%%% Nothing %%%%%
% Search directory (normal datasets)
names = dir('datasets\nothing'); % Input datasets
names = {names.name};
names = names(3:length(names));

% Pre-allocation for speed
counter = 1;

% Define the wavelet family and the level for decomposition
waveletFunction='db4';
level=5;

% Load normal datasets, partitioning and filtering
for i = 1:length(names)
    raw = load(['datasets\nothing\' names{i}]);
    raw = raw(1:Fst, :); % Take datasets from 1-60 seconds
    
    % Filtering
    % Bandpass Filter
    Fn = Fs/2;                         % Nyquist Frequency (Hz)
    Wp = [0.51 64]/Fn;                 % Passband Frequencies (Normalised)
    Ws = [0.5 64.01]/Fn;               % Stopband Frequencies (Normalised)
    Rp = 5;                            % Passband Ripple (dB)
    Rs = 30;                           % Stopband Ripple (dB)
    [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);    % Filter Order
    [z,p,k] = cheby2(n,Rs,Ws);         % Filter Design
    [sosbp,gbp] = zp2sos(z,p,k);       % Convert To Second-Order-Section For Stability
    raw = filtfilt(sosbp,gbp, raw);    % Filter Signal
    
    % Partition data
    EEG{1} = raw(1:20*Fs, :);          % Take datasets from 0-20 seconds
    EEG{2} = raw((20*Fs)+1:40*Fs, :);  % Take datasets from 21-40 seconds
    EEG{3} = raw((40*Fs)+1:Fst, :);    % Take datasets from 41-60 seconds
    
    % EEG Sub-band Extraction & Feature Extraction
    for ii = 1:3
        % Discrete Fourier Transform
        EEGfft = fft(EEG{1, ii});
        EEGfftshift = fftshift(EEGfft);
        EEGfftshift = abs(EEGfftshift(length(EEGfftshift)/2:length(EEGfftshift), :));
        
        % Multilevel 1-D wavelet decomposition for both channels
        [C,L]=wavedec(EEG{1, ii}(:, 1),level,waveletFunction);
        [D,L]=wavedec(EEG{1, ii}(:, 2),level,waveletFunction);
        
        % Decompose into EEG Sub-bands
        for iii = 1:2
            if iii == 1
            % Reconstruct single branch from 1-D wavelet coefficients for ch0
                Gamma(:, iii) = wrcoef('d',C,L,waveletFunction,2); %GAMMA, 32 - 64
                Beta(:, iii) = wrcoef('d',C,L,waveletFunction,3); %BETA, 16 - 32
                Alpha(:, iii) = wrcoef('d',C,L,waveletFunction,4); %ALPHA, 8 - 16
                Theta(:, iii) = wrcoef('d',C,L,waveletFunction,5); %THETA, 4 - 8
                Delta(:, iii) = wrcoef('a',C,L,waveletFunction,5); %DELTA, 0 - 4
            elseif iii == 2
            % Reconstruct single branch from 1-D wavelet coefficients for ch1
                Gamma(:, iii) = wrcoef('d',D,L,waveletFunction,2); %GAMMA, 32 - 64
                Beta(:, iii) = wrcoef('d',D,L,waveletFunction,3); %BETA, 16 - 32
                Alpha(:, iii) = wrcoef('d',D,L,waveletFunction,4); %ALPHA, 8 - 16
                Theta(:, iii) = wrcoef('d',D,L,waveletFunction,5); %THETA, 4 - 8
                Delta(:, iii) = wrcoef('a',D,L,waveletFunction,5); %DELTA, 0 - 4
            end
        end
        
        % FFT
        % Perform FFT on each band, then extract features from the FFTs
        Gamma_FFT = fft(Gamma);
        Gammafftshift = fftshift(Gamma_FFT);
        Gammafftshift = abs(Gammafftshift(length(Gammafftshift)/2:length(Gammafftshift), :));
        
        Beta_FFT = fft(Beta);
        Betafftshift = fftshift(Beta_FFT);
        Betafftshift = abs(Betafftshift(length(Betafftshift)/2:length(Betafftshift), :));
        
        Alpha_FFT = fft(Alpha);
        Alphafftshift = fftshift(Alpha_FFT);
        Alphafftshift = abs(Alphafftshift(length(Alphafftshift)/2:length(Alphafftshift), :));
        
        Theta_FFT = fft(Theta);
        Thetafftshift = fftshift(Theta_FFT);
        Thetafftshift = abs(Thetafftshift(length(Thetafftshift)/2:length(Thetafftshift), :));
        
        Delta_FFT = fft(Delta);
        Deltafftshift = fftshift(Delta_FFT);
        Deltafftshift = abs(Deltafftshift(length(Deltafftshift)/2:length(Deltafftshift), :));
        
        % Labelling Classes
        counter = counter + 1;
        action(tablegen, 1) = categorical(cellstr('Nothing'));
        
        % Spectral Power
        power_gamma(tablegen, :) = sum(Gammafftshift);
        power_beta(tablegen, :) = sum(Betafftshift);
        power_alpha(tablegen, :) = sum(Alphafftshift);
        power_theta(tablegen, :) = sum(Thetafftshift);
        power_delta(tablegen, :) = sum(Deltafftshift);
        powerdiff_gamma(tablegen, :) = abs(power_gamma(tablegen, 1) - power_gamma(tablegen, 2));
        powerdiff_beta(tablegen, :) = abs(power_beta(tablegen, 1) - power_beta(tablegen, 2));
        powerdiff_alpha(tablegen, :) = abs(power_alpha(tablegen, 1) - power_alpha(tablegen, 2));
        powerdiff_theta(tablegen, :) = abs(power_theta(tablegen, 1) - power_theta(tablegen, 2));
        powerdiff_delta(tablegen, :) = abs(power_delta(tablegen, 1) - power_delta(tablegen, 2));

        % DASDV
        for row = 2:length(EEG{1, 1})
            for iii = 1:2 % Both Channels
                rowcalc_gamma(row, iii) = Gamma(row, iii) - Gamma(row-1, iii);
                rowcalc_beta(row, iii) = Beta(row, iii) - Beta(row-1, iii);
                rowcalc_alpha(row, iii) = Alpha(row, iii) - Alpha(row-1, iii);
                rowcalc_theta(row, iii) = Theta(row, iii) - Theta(row-1, iii);
                rowcalc_delta(row, iii) = Delta(row, iii) - Delta(row-1, iii);
            end
        end
        dasdv_gamma(tablegen, :) = sqrt(abs(sum(rowcalc_beta.^2)/Fst-1));
        dasdv_beta(tablegen, :) = sqrt(abs(sum(rowcalc_beta.^2)/Fst-1));
        dasdv_alpha(tablegen, :) = sqrt(abs(sum(rowcalc_alpha.^2)/Fst-1));
        dasdv_theta(tablegen, :) = sqrt(abs(sum(rowcalc_theta.^2)/Fst-1));
        dasdv_delta(tablegen, :) = sqrt(abs(sum(rowcalc_delta.^2)/Fst-1));
        
        % Other Simple Features
        meanamp(tablegen, :) = mean(EEGfftshift); % Mean
        maxamp(tablegen, :) = max(EEGfftshift); % Max
        varamp(tablegen, :) = var(EEGfftshift); % Variance
        
        % Next row
        tablegen = tablegen + 1;
    end
end

% Table Generation
ClassName = table(action, power_gamma(:, 1), power_gamma(:, 2), power_beta(:, 1), power_beta(:, 2), power_alpha(:, 1), power_alpha(:, 2), power_theta(:, 1), power_theta(:, 2), power_delta(:, 1), power_delta(:, 2), dasdv_gamma(:, 1), dasdv_gamma(:, 2), dasdv_beta(:, 1), dasdv_beta(:, 2), dasdv_alpha(:, 1), dasdv_alpha(:, 2), dasdv_theta(:, 1), dasdv_theta(:, 2), dasdv_delta(:, 1), dasdv_delta(:, 2), powerdiff_gamma(:, 1), powerdiff_beta(:, 1), powerdiff_alpha(:, 1), powerdiff_theta(:, 1), powerdiff_delta(:, 1), meanamp(:, 1), meanamp(:, 2), maxamp(:, 1), maxamp(:, 2), varamp(:, 1), varamp(:, 2),...
        'VariableNames',{'Action', 'Gamma_powerA', 'Gamma_powerB', 'Beta_powerA', 'Beta_powerB', 'Alpha_powerA', 'Alpha_powerB', 'Theta_powerA', 'Theta_powerB', 'Delta_powerA', 'Delta_powerB', 'Gamma_DASDVA', 'Gamma_DASDVB', 'Beta_DASDVA', 'Beta_DASDVB', 'Alpha_DASDVA', 'Alpha_DASDVB', 'Theta_DASDVA', 'Theta_DASDVB', 'Delta_DASDVA', 'Delta_DASDVB', 'PowerDiff_Gamma', 'PowerDiff_Beta', 'PowerDiff_Alpha', 'PowerDiff_Theta', 'PowerDiff_Delta', 'MeanA', 'MeanB', 'MaxA', 'MaxB', 'VarianceA', 'VarianceB'});
