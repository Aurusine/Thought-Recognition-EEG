clear; 

Fs = 256; % sampling rate
Fst = Fs*60; % total samples

%%%%%%%%%%%%% Create Datasets %%%%%%%%%%%%%

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
    Wp = [0.51 32]/Fn;                 % Passband Frequencies (Normalised)
    Ws = [0.5 32.01]/Fn;               % Stopband Frequencies (Normalised)
    Rp = 5;                            % Passband Ripple (dB)
    Rs = 50;                           % Stopband Ripple (dB)
    [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);    % Filter Order
    [z,p,k] = cheby2(n,Rs,Ws);         % Filter Design
    [sosbp,gbp] = zp2sos(z,p,k);       % Convert To Second-Order-Section For Stability
    raw = filtfilt(sosbp,gbp, raw);    % Filter Signal
    
    % Partition data
    EEG{1} = raw(Fs+1:20*Fs, :); % Take datasets from 1-20 seconds
    EEG{2} = raw((20*Fs)+1:40*Fs, :); % Take datasets from 21-40 seconds
    EEG{3} = raw((40*Fs)+1:Fst, :); % Take datasets from 41-60 seconds
    
    % EEG Sub-band Extraction & Feature Extraction
    for ii = 1:3
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
        cube(counter, :) = FFT(1, :)'; % 64 features/rows, 90 datasets/columns
        counter = counter + 1;
        action(tablegen, 1) = categorical(cellstr('Cube')); % Labelling

        % Root-Mean-Square
        rms(tablegen, :) = sqrt(mean(EEG{1, ii}.^2));
        
        meanamp(tablegen, :) = mean(EEG{1, ii}); % Mean
        maxamp(tablegen, :) = max(EEG{1, ii}); % Max
        minamp(tablegen, :) = min(EEG{1, ii}); % Min
        varamp(tablegen, :) = var(EEG{1, ii}); % Variance

        % DASDV
        for row = 2:length(EEG{1, 1})
            for iii = 1:2 % Both Channels
                rowcalc(row, iii) = EEG{1, ii}(row, iii) - EEG{1, ii}(row-1, iii);
            end
        end
        dasdv(tablegen, :) = sqrt(abs(sum(rowcalc.^2)/Fst-1));

        tablegen = tablegen + 1; % Next row
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Search directory (normal datasets)
names = dir('datasets\math'); % Input datasets
names = {names.name};
names = names(3:length(names));

% Pre-allocation for speed
counter = 1;

% Load normal datasets, partitioning and filtering
for i = 1:length(names)
    raw = load(['datasets\math\' names{i}]);
    raw = raw(1:Fst, :); % Take datasets from 1-60 seconds
    
    % Filtering
    % Bandpass Filter
    Fn = Fs/2;                         % Nyquist Frequency (Hz)
    Wp = [0.51 32]/Fn;                 % Passband Frequencies (Normalised)
    Ws = [0.5 32.01]/Fn;               % Stopband Frequencies (Normalised)
    Rp = 5;                            % Passband Ripple (dB)
    Rs = 50;                           % Stopband Ripple (dB)
    [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);    % Filter Order
    [z,p,k] = cheby2(n,Rs,Ws);         % Filter Design
    [sosbp,gbp] = zp2sos(z,p,k);       % Convert To Second-Order-Section For Stability
    raw = filtfilt(sosbp,gbp, raw);    % Filter Signal
    
    % Partition data
    EEG{1} = raw(Fs+1:20*Fs, :); % Take datasets from 1-20 seconds
    EEG{2} = raw((20*Fs)+1:40*Fs, :); % Take datasets from 21-40 seconds
    EEG{3} = raw((40*Fs)+1:Fst, :); % Take datasets from 41-60 seconds
    
    % FFT
    for ii = 1:3
        EEGfft = fft(EEG{1, ii});
        EEGfftshift = fftshift(EEGfft);
        EEGfftshift = abs(EEGfftshift(length(EEGfftshift)/2:length(EEGfftshift), :));
        for bin = 1:32 % 32 frequencies
            fftbin{bin} = EEGfftshift(((25*(bin-1))+1):(25*bin), :); % 1 frequency = 1 bin
            for iii = 1:2 % Both Channels
                if iii == 1
                    temp(bin, iii) = mean(fftbin{bin}(:, iii)); % 20 features/columns each for 2 rows/channels
                    FFT(:, bin) = horzcat(temp(bin, iii)); % Partition into Channel A
                elseif iii == 2
                    temp(bin, iii) = mean(fftbin{bin}(:, iii)); % 20 features/columns each for 2 rows/channels
                    FFT(:, bin+32) = horzcat(temp(bin, iii)); % Partition into Channel B
                end
            end
        end
        math(counter, :) = FFT(1, :)'; % 64 features/rows, 90 datasets/columns
        counter = counter + 1;
        action(tablegen, 1) = categorical(cellstr('Math')); % Labelling
    
        % Root-Mean-Square
        rms(tablegen, :) = sqrt(mean(EEG{1, ii}.^2));
        meanamp(tablegen, :) = mean(EEG{1, ii}); % Mean
        maxamp(tablegen, :) = max(EEG{1, ii}); % Max
        minamp(tablegen, :) = min(EEG{1, ii}); % Min
        varamp(tablegen, :) = var(EEG{1, ii}); % Variance

        % DASDV
        for row = 2:length(EEG{1, 1})
            for iii = 1:2 % Both Channels
                rowcalc(row, iii) = EEG{1, ii}(row, iii) - EEG{1, ii}(row-1, iii);
            end
        end
        dasdv(tablegen, :) = sqrt(abs(sum(rowcalc.^2)/Fst-1));

        tablegen = tablegen + 1; % Next row
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Search directory (normal datasets)
names = dir('datasets\letter'); % Input datasets
names = {names.name};
names = names(3:length(names));

% Pre-allocation for speed
counter = 1;

% Load normal datasets, partitioning and filtering
for i = 1:length(names)
    raw = load(['datasets\letter\' names{i}]);
    raw = raw(1:Fst, :); % Take datasets from 1-60 seconds
    
    % Filtering
    % Bandpass Filter
    Fn = Fs/2;                         % Nyquist Frequency (Hz)
    Wp = [0.51 32]/Fn;                 % Passband Frequencies (Normalised)
    Ws = [0.5 32.01]/Fn;               % Stopband Frequencies (Normalised)
    Rp = 5;                            % Passband Ripple (dB)
    Rs = 50;                           % Stopband Ripple (dB)
    [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);    % Filter Order
    [z,p,k] = cheby2(n,Rs,Ws);         % Filter Design
    [sosbp,gbp] = zp2sos(z,p,k);       % Convert To Second-Order-Section For Stability
    raw = filtfilt(sosbp,gbp, raw);    % Filter Signal

    % Partition data
    EEG{1} = raw(Fs+1:20*Fs, :); % Take datasets from 1-20 seconds
    EEG{2} = raw((20*Fs)+1:40*Fs, :); % Take datasets from 21-40 seconds
    EEG{3} = raw((40*Fs)+1:Fst, :); % Take datasets from 41-60 seconds
    
    % FFT
    for ii = 1:3
        EEGfft = fft(EEG{1, ii});
        EEGfftshift = fftshift(EEGfft);
        EEGfftshift = abs(EEGfftshift(length(EEGfftshift)/2:length(EEGfftshift), :));
        for bin = 1:32 % 32 frequencies
            fftbin{bin} = EEGfftshift(((25*(bin-1))+1):(25*bin), :); % 1 frequency = 1 bin
            for iii = 1:2 % Both Channels
                if iii == 1
                    temp(bin, iii) = mean(fftbin{bin}(:, iii)); % 20 features/columns each for 2 rows/channels
                    FFT(:, bin) = horzcat(temp(bin, iii)); % Partition into Channel A
                elseif iii == 2
                    temp(bin, iii) = mean(fftbin{bin}(:, iii)); % 20 features/columns each for 2 rows/channels
                    FFT(:, bin+32) = horzcat(temp(bin, iii)); % Partition into Channel B
                end
            end
        end
    letter(counter, :) = FFT(1, :)'; % 64 features/rows, 90 datasets/columns
    counter = counter + 1;
    action(tablegen, 1) = categorical(cellstr('Letter')); % Labelling
    
        % Root-Mean-Square
        rms(tablegen, :) = sqrt(mean(EEG{1, ii}.^2));
        
        meanamp(tablegen, :) = mean(EEG{1, ii}); % Mean
        maxamp(tablegen, :) = max(EEG{1, ii}); % Max
        minamp(tablegen, :) = min(EEG{1, ii}); % Min
        varamp(tablegen, :) = var(EEG{1, ii}); % Variance

        % DASDV
        for row = 2:length(EEG{1, 1})
            for iii = 1:2 % Both Channels
                rowcalc(row, iii) = EEG{1, ii}(row, iii) - EEG{1, ii}(row-1, iii);
            end
        end
    dasdv(tablegen, :) = sqrt(abs(sum(rowcalc.^2)/Fst-1));
    
    tablegen = tablegen + 1; % Next row
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Search directory (normal datasets)
names = dir('datasets\nothing'); % Input datasets
names = {names.name};
names = names(3:length(names));

% Pre-allocation for speed
counter = 1;

% Load normal datasets, partitioning and filtering
for i = 1:length(names)
    raw = load(['datasets\nothing\' names{i}]);
    raw = raw(1:Fst, :); % Take datasets from 1-60 seconds
    
    % Filtering
    % Bandpass Filter
    Fn = Fs/2;                         % Nyquist Frequency (Hz)
    Wp = [0.51 32]/Fn;                 % Passband Frequencies (Normalised)
    Ws = [0.5 32.01]/Fn;               % Stopband Frequencies (Normalised)
    Rp = 5;                            % Passband Ripple (dB)
    Rs = 50;                           % Stopband Ripple (dB)
    [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);    % Filter Order
    [z,p,k] = cheby2(n,Rs,Ws);         % Filter Design
    [sosbp,gbp] = zp2sos(z,p,k);       % Convert To Second-Order-Section For Stability
    raw = filtfilt(sosbp,gbp, raw);    % Filter Signal
    
    % Partition data
    EEG{1} = raw(Fs+1:20*Fs, :); % Take datasets from 1-20 seconds
    EEG{2} = raw((20*Fs)+1:40*Fs, :); % Take datasets from 21-40 seconds
    EEG{3} = raw((40*Fs)+1:Fst, :); % Take datasets from 41-60 seconds
    
    % FFT
    for ii = 1:3
        EEGfft = fft(EEG{1, ii});
        EEGfftshift = fftshift(EEGfft);
        EEGfftshift = abs(EEGfftshift(length(EEGfftshift)/2:length(EEGfftshift), :));
        for bin = 1:32 % 32 frequencies
            fftbin{bin} = EEGfftshift(((25*(bin-1))+1):(25*bin), :); % 1 frequency = 1 bin
            for iii = 1:2 % Both Channels
                if iii == 1
                    temp(bin, iii) = mean(fftbin{bin}(:, iii)); % 20 features/columns each for 2 rows/channels
                    FFT(:, bin) = horzcat(temp(bin, iii)); % Partition into Channel A
                elseif iii == 2
                    temp(bin, iii) = mean(fftbin{bin}(:, iii)); % 20 features/columns each for 2 rows/channels
                    FFT(:, bin+32) = horzcat(temp(bin, iii)); % Partition into Channel B
                end
            end
        end
    nothing(counter, :) = FFT(1, :)'; % 64 features/rows, 90 datasets/columns
    counter = counter + 1;
    action(tablegen, 1) = categorical(cellstr('Nothing')); % Labelling
    
    % Root-Mean-Square
    rms(counter, :) = sqrt(mean(EEGfftshift{1, ii}.^2));
    meanamp(counter, :) = mean(EEGfftshift{1, ii}); % Mean
    maxamp(counter, :) = max(EEGfftshift{1, ii}); % Max
    minamp(counter, :) = min(EEGfftshift{1, ii}); % Min
    varamp(counter, :) = var(EEGfftshift{1, ii}); % Variance

        % DASDV
        for row = 2:length(EEG{1, 1})
            for iii = 1:2 % Both Channels
                rowcalc(row, iii) = EEG{1, ii}(row, iii) - EEG{1, ii}(row-1, iii);
            end
        end
    dasdv(tablegen, :) = sqrt(abs(sum(rowcalc.^2)/Fst-1));
    
    tablegen = tablegen + 1; % Next row
    end
end

% Table Generation
FFTBin = vertcat(cube, math, letter, nothing);
FFTTable = array2table(FFTBin);
ClassName = table(action, rms(:, 1), rms(:, 2), meanamp(:, 1), meanamp(:, 2), maxamp(:, 1), maxamp(:, 2), minamp(:, 1), minamp(:, 2), varamp(:, 1), varamp(:, 2), dasdv(:, 1), dasdv(:, 2),...
        'VariableNames',{'Action', 'RMSA', 'RMSB', 'MeanA', 'MeanB', 'MaxA', 'MaxB', 'MinA', 'MinB', 'VarianceA', 'VarianceB', 'DASDVA', 'DASDVB'});
    
Final = [ClassName FFTTable]