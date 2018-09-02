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
    Rs = 50;                           % Stopband Ripple (dB)
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
        
        for bin = 1:32 % 32 frequencies
            fftbin{bin} = EEGfftshift(((25*(bin-1))+1):(25*bin), :); % 1 frequency = 1 bin
            for iii = 1:2 % Both Channels
                if iii == 1
                    % FFT Feature Extraction
                    temp(bin, iii) = mean(fftbin{bin}(:, iii)); % 20 features/columns each for 2 rows/channels
                    FFT(:, bin) = horzcat(temp(bin, iii)); % Partition into Channel A
                elseif iii == 2
                    % FFT Feature Extraction
                    temp(bin, iii) = mean(fftbin{bin}(:, iii)); % 20 features/columns each for 2 rows/channels
                    FFT(:, bin+32) = horzcat(temp(bin, iii)); % Partition into Channel B
                end
            end
        end
        cube(counter, :) = FFT(1, :)'; % 64 features/rows, 90 datasets/columns
        counter = counter + 1;
        action(tablegen, 1) = categorical(cellstr('Cube')); % Labelling

        % Root-Mean-Square
        rms_gamma(tablegen, :) = sqrt(mean(Gamma.^2));
        rms_beta(tablegen, :) = sqrt(mean(Beta.^2));
        rms_alpha(tablegen, :) = sqrt(mean(Alpha.^2));
        rms_theta(tablegen, :) = sqrt(mean(Theta.^2));
        rms_delta(tablegen, :) = sqrt(mean(Delta.^2));

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

        tablegen = tablegen + 1; % Next row
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
    Rs = 50;                           % Stopband Ripple (dB)
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
        
        EEGfft = fft(EEG{1, ii});
        EEGfftshift = fftshift(EEGfft);
        EEGfftshift = abs(EEGfftshift(length(EEGfftshift)/2:length(EEGfftshift), :));
        for bin = 1:32 % 32 frequencies
            fftbin{bin} = EEGfftshift(((25*(bin-1))+1):(25*bin), :); % 1 frequency = 1 bin
            for iii = 1:2 % Both Channels
                if iii == 1
                    % FFT Feature Extraction
                    temp(bin, iii) = mean(fftbin{bin}(:, iii)); % 20 features/columns each for 2 rows/channels
                    FFT(:, bin) = horzcat(temp(bin, iii)); % Partition into Channel A
                elseif iii == 2
                    % FFT Feature Extraction
                    temp(bin, iii) = mean(fftbin{bin}(:, iii)); % 20 features/columns each for 2 rows/channels
                    FFT(:, bin+32) = horzcat(temp(bin, iii)); % Partition into Channel B
                end
            end
        end
        math(counter, :) = FFT(1, :)'; % 64 features/rows, 90 datasets/columns
        counter = counter + 1;
        action(tablegen, 1) = categorical(cellstr('Math')); % Labelling

        % Root-Mean-Square
        rms_gamma(tablegen, :) = sqrt(mean(Gamma.^2));
        rms_beta(tablegen, :) = sqrt(mean(Beta.^2));
        rms_alpha(tablegen, :) = sqrt(mean(Alpha.^2));
        rms_theta(tablegen, :) = sqrt(mean(Theta.^2));
        rms_delta(tablegen, :) = sqrt(mean(Delta.^2));

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

        tablegen = tablegen + 1; % Next row
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
    Rs = 50;                           % Stopband Ripple (dB)
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
        
        EEGfft = fft(EEG{1, ii});
        EEGfftshift = fftshift(EEGfft);
        EEGfftshift = abs(EEGfftshift(length(EEGfftshift)/2:length(EEGfftshift), :));
        for bin = 1:32 % 32 frequencies
            fftbin{bin} = EEGfftshift(((25*(bin-1))+1):(25*bin), :); % 1 frequency = 1 bin
            for iii = 1:2 % Both Channels
                if iii == 1
                    % FFT Feature Extraction
                    temp(bin, iii) = mean(fftbin{bin}(:, iii)); % 20 features/columns each for 2 rows/channels
                    FFT(:, bin) = horzcat(temp(bin, iii)); % Partition into Channel A
                elseif iii == 2
                    % FFT Feature Extraction
                    temp(bin, iii) = mean(fftbin{bin}(:, iii)); % 20 features/columns each for 2 rows/channels
                    FFT(:, bin+32) = horzcat(temp(bin, iii)); % Partition into Channel B
                end
            end
        end
        letter(counter, :) = FFT(1, :)'; % 64 features/rows, 90 datasets/columns
        counter = counter + 1;
        action(tablegen, 1) = categorical(cellstr('Letter')); % Labelling

        % Root-Mean-Square
        rms_gamma(tablegen, :) = sqrt(mean(Gamma.^2));
        rms_beta(tablegen, :) = sqrt(mean(Beta.^2));
        rms_alpha(tablegen, :) = sqrt(mean(Alpha.^2));
        rms_theta(tablegen, :) = sqrt(mean(Theta.^2));
        rms_delta(tablegen, :) = sqrt(mean(Delta.^2));

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

        tablegen = tablegen + 1; % Next row
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
    Rs = 50;                           % Stopband Ripple (dB)
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
        
        EEGfft = fft(EEG{1, ii});
        EEGfftshift = fftshift(EEGfft);
        EEGfftshift = abs(EEGfftshift(length(EEGfftshift)/2:length(EEGfftshift), :));
        for bin = 1:32 % 32 frequencies
            fftbin{bin} = EEGfftshift(((25*(bin-1))+1):(25*bin), :); % 1 frequency = 1 bin
            for iii = 1:2 % Both Channels
                if iii == 1
                    % FFT Feature Extraction
                    temp(bin, iii) = mean(fftbin{bin}(:, iii)); % 20 features/columns each for 2 rows/channels
                    FFT(:, bin) = horzcat(temp(bin, iii)); % Partition into Channel A
                elseif iii == 2
                    % FFT Feature Extraction
                    temp(bin, iii) = mean(fftbin{bin}(:, iii)); % 20 features/columns each for 2 rows/channels
                    FFT(:, bin+32) = horzcat(temp(bin, iii)); % Partition into Channel B
                end
            end
        end
        nothing(counter, :) = FFT(1, :)'; % 64 features/rows, 90 datasets/columns
        counter = counter + 1;
        action(tablegen, 1) = categorical(cellstr('Nothing')); % Labelling

        % Root-Mean-Square
        rms_gamma(tablegen, :) = sqrt(mean(Gamma.^2));
        rms_beta(tablegen, :) = sqrt(mean(Beta.^2));
        rms_alpha(tablegen, :) = sqrt(mean(Alpha.^2));
        rms_theta(tablegen, :) = sqrt(mean(Theta.^2));
        rms_delta(tablegen, :) = sqrt(mean(Delta.^2));

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

        tablegen = tablegen + 1; % Next row
    end
end

% Table Generation
FFTBin = vertcat(cube, math, letter, nothing);
FFTTable = array2table(FFTBin);
ClassName = table(action, rms_gamma(:, 1), rms_gamma(:, 2), rms_beta(:, 1), rms_beta(:, 2), rms_alpha(:, 1), rms_alpha(:, 2), rms_theta(:, 1), rms_theta(:, 2), rms_delta(:, 1), rms_delta(:, 2), dasdv_gamma(:, 1), dasdv_gamma(:, 2), dasdv_beta(:, 1), dasdv_beta(:, 2), dasdv_alpha(:, 1), dasdv_alpha(:, 2), dasdv_theta(:, 1), dasdv_theta(:, 2), dasdv_delta(:, 1), dasdv_delta(:, 2),...
        'VariableNames',{'Action', 'Gamma_RMSA', 'Gamma_RMSB', 'Beta_RMSA', 'Beta_RMSB', 'Alpha_RMSA', 'Alpha_RMSB', 'Theta_RMSA', 'Theta_RMSB', 'Delta_RMSA', 'Delta_RMSB', 'Gamma_DASDVA', 'Gamma_DASDVB', 'Beta_DASDVA', 'Beta_DASDVB', 'Alpha_DASDVA', 'Alpha_DASDVB', 'Theta_DASDVA', 'Theta_DASDVB', 'Delta_DASDVA', 'Delta_DASDVB'});

AFinal = [ClassName FFTTable]
