clear; 

fs = 2048; % sampling rate
fst = fs*10; % total samples

%%%%%%%%%%%%% Cylindrical Grip %%%%%%%%%%%%%

% Load files (normal datasets)
names = dir('Cylin** e1 e2.txt');
names = {names.name};

% Pre-allocation for speed
rmsc_raw = zeros(2, length(names));
mavc = zeros(2, length(names));
actionc(1:length(names), :) = categorical(cellstr('Cylindrical')); % Labelling
rowcalc = zeros(fst, 2);
dasdvc = zeros(length(names), 2);

% Segmentation and data modification
for ii = 1:length(names)
    raw = load([names{ii}]);
    raw = raw(1:fst, :);
    rmsc_raw(:, ii) = sqrt(mean(raw.^2)); % Absolute RMS values
    actionc(ii, :) = categorical(cellstr('Cylinder')); % Labelling
end

% Load files (maximum grip force datasets)
maxnames = dir('Max*Cylin e1 e2.txt');
maxnames = {maxnames.name};

% Pre-allocation for speed
rmsc_max = zeros(2, length(maxnames));
GFPC = zeros(1, length(maxnames));

% Calculation of reference GFP
for ii = 1:length(maxnames)
    maxraw = load([maxnames{ii}]);
    maxraw = maxraw(1:fst, :);
    rmsc_max(:, ii) = sqrt(mean(maxraw.^2)); % Row 1 = J, Row 2 = M, Row 3 = T
    % Grip-Force Percentage calculation
        if strcmp(maxnames{ii}(1:9), 'MaxJCylin');
            GFPC(:, ii) = mean(rmsc_raw(:, ii)./rmsc_max(:, 1)); % GFPC(:, 1) = J
        elseif strcmp(maxnames{ii}(1:9), 'MaxMCylin');
            GFPC(:, ii) = mean(rmsc_raw(:, ii)./rmsc_max(:, 2)); % GFPC(:, 2) = M
        elseif strcmp(maxnames{ii}(1:9), 'MaxTCylin');
            GFPC(:, ii) = mean(rmsc_raw(:, ii)./rmsc_max(:, 3)); % GFPC(:, 3) = T
        end
end

% Pre-allocation for speed
raw_normalized = zeros(fst, 2);
rmsc = zeros(2, length(names));
ffttemp = zeros(2, 20);
meanfreqc = zeros(length(names), 2);

% Normalization of other subject's values based on reference GFP
for ii = 1:length(names)
    raw = load([names{ii}]);
    raw = raw(1:fst, :);
    for i = 1:fst
        if strcmp(names{ii}(1:6), 'CylinJ');
            raw_normalized(i, :) = raw(i, :)*(GFPC(:, 2)/GFPC(:, 1));
        elseif strcmp(names{ii}(1:6), 'CylinM');
            raw_normalized(i, :) = raw(i, :); % Using M as reference
        elseif strcmp(names{ii}(1:6), 'CylinT');
            raw_normalized(i, :) = raw(i, :)*(GFPC(:, 2)/GFPC(:, 3));
        end
    end
    
    % Feature extraction
    mavc(:, ii) = abs(sum(raw_normalized)/fst); % Mean Absolute Value
    rmsc(:, ii) = sqrt(mean(raw_normalized.^2)); % Root-Mean-Square
    % DASDV
    for row = 2:fst
        for i = 1:2 % Both Channels
            rowcalc(row, i) = raw_normalized(row, i) - raw_normalized(row-1, i);
        end
    end
    dasdvc(ii, :) = sqrt(abs(sum(rowcalc.^2)/fst-1));
    % FFT
    rawfft = fft(raw_normalized);
    rawfftshift = fftshift(rawfft);
    rawfftshift = abs(rawfftshift(fst-5000:fst, :));
    for i = 1:2 % Both Channels
        meanfreqc(ii, i) = mean(rawfftshift(:, i));
    end
end

%%%%%%%%%%%%% Key Grip %%%%%%%%%%%%%

% Load files (normal datasets)
names = dir('Key** e1 e2.txt');
names = {names.name};

% Pre-allocation for speed
rmsk_raw = zeros(2, length(names));
mavk = zeros(2, length(names));
actionk(1:length(names), :) = categorical(cellstr('Key')); % Labelling
rowcalc = zeros(fst, 2);
dasdvk = zeros(length(names), 2);

% Segmentation and data modification
for ii = 1:length(names)
    raw = load([names{ii}]);
    raw = raw(1:fst, :);
    rmsk_raw(:, ii) = sqrt(mean(raw.^2)); % Absolute RMS values
    actionk(ii, :) = categorical(cellstr('Key')); % Labelling
end

% Load files (maximum grip force datasets)
maxnames = dir('Max*Key e1 e2.txt');
maxnames = {maxnames.name};

% Pre-allocation for speed
rmsk_max = zeros(2, length(maxnames));
GFPK = zeros(1, length(maxnames));

% Calculation of reference GFP
for ii = 1:length(maxnames)
    maxraw = load([maxnames{ii}]);
    maxraw = maxraw(1:fst, :);
    rmsk_max(:, ii) = sqrt(mean(maxraw.^2)); % Row 1 = J, Row 2 = M, Row 3 = T
    % Grip-Force Percentage calculation
        if strcmp(maxnames{ii}(1:7), 'MaxJKey');
            GFPK(:, ii) = mean(rmsk_raw(:, ii)./rmsk_max(:, 1)); % GFPK(:, 1) = J
        elseif strcmp(maxnames{ii}(1:7), 'MaxMKey');
            GFPK(:, ii) = mean(rmsk_raw(:, ii)./rmsk_max(:, 2)); % GFPK(:, 2) = M
        elseif strcmp(maxnames{ii}(1:7), 'MaxTKey');
            GFPK(:, ii) = mean(rmsk_raw(:, ii)./rmsk_max(:, 3)); % GFPK(:, 3) = T
        end
end

% Pre-allocation for speed
raw_normalized = zeros(fst, 2);
rmsk = zeros(2, length(names));
meanfreqk = zeros(length(names), 2);

% Normalization of other subject's values based on reference GFP
for ii = 1:length(names)
    raw = load([names{ii}]);
    raw = raw(1:fst, :);
    for i = 1:fst
        if strcmp(names{ii}(1:4), 'KeyJ');
            raw_normalized(i, :) = raw(i, :)*(GFPK(:, 2)/GFPK(:, 1));
        elseif strcmp(names{ii}(1:4), 'KeyM');
            raw_normalized(i, :) = raw(i, :); % Using M as reference
        elseif strcmp(names{ii}(1:4), 'KeyT');
            raw_normalized(i, :) = raw(i, :)*(GFPK(:, 2)/GFPK(:, 3));
        end
    end
    
    % Feature extraction
    mavk(:, ii) = abs(sum(raw_normalized)/fst); % Mean Absolute Value
    rmsk(:, ii) = sqrt(mean(raw_normalized.^2)); % Root-Mean-Square
    % DASDV
    for row = 2:fst
        for i = 1:2 % Both Channels
            rowcalc(row, i) = raw_normalized(row, i) - raw_normalized(row-1, i);
        end
    end
    dasdvk(ii, :) = sqrt(abs(sum(rowcalc.^2)/fst-1));
    % FFT
    rawfft = fft(raw_normalized);
    rawfftshift = fftshift(rawfft);
    rawfftshift = abs(rawfftshift(fst-5000:fst, :));
    for i = 1:2 % Both Channels
        meanfreqk(ii, i) = mean(rawfftshift(:, i));
    end
end

%%%%%%%%%%%%% Pinch Grip %%%%%%%%%%%%%

% Load files (normal datasets)
names = dir('Pinch** e1 e2.txt');
names = {names.name};

% Pre-allocation for speed
rmsp_raw = zeros(2, length(names));
mavp = zeros(2, length(names));
actionp(1:length(names), :) = categorical(cellstr('Pinch')); % Labelling
rowcalc = zeros(fst, 2);
dasdvp = zeros(length(names), 2);

% Segmentation and data modification
for ii = 1:length(names)
    raw = load([names{ii}]);
    raw = raw(1:fst, :);
    rmsp_raw(:, ii) = sqrt(mean(raw.^2)); % Absolute RMS values
end

% Load files (maximum grip force datasets)
maxnames = dir('Max*Pinch e1 e2.txt');
maxnames = {maxnames.name};

% Pre-allocation for speed
rmsp_max = zeros(2, length(maxnames));
GFPP = zeros(1, length(maxnames));

% Calculation of reference GFP
for ii = 1:length(maxnames)
    maxraw = load([maxnames{ii}]);
    maxraw = maxraw(1:fst, :);
    rmsp_max(:, ii) = sqrt(mean(maxraw.^2)); % Row 1 = J, Row 2 = M, Row 3 = T
    % Grip-Force Percentage calculation
        if strcmp(maxnames{ii}(1:9), 'MaxJPinch');
            GFPP(:, ii) = mean(rmsp_raw(:, ii)./rmsp_max(:, 1)); % GFPP(:, 1) = J
        elseif strcmp(maxnames{ii}(1:9), 'MaxMPinch');
            GFPP(:, ii) = mean(rmsp_raw(:, ii)./rmsp_max(:, 2)); % GFPP(:, 2) = M
        elseif strcmp(maxnames{ii}(1:9), 'MaxTPinch');
            GFPP(:, ii) = mean(rmsp_raw(:, ii)./rmsp_max(:, 3)); % GFPP(:, 3) = T
        end
end

% Pre-allocation for speed
raw_normalized = zeros(fst, 2);
rmsp = zeros(2, length(names));
meanfreqp = zeros(length(names), 2);

% Normalization of other subject's values based on reference GFP
for ii = 1:length(names)
    raw = load([names{ii}]);
    raw = raw(1:fst, :);
    for i = 1:fst
        if strcmp(names{ii}(1:6), 'PinchJ');
            raw_normalized(i, :) = raw(i, :)*(GFPP(:, 2)/GFPP(:, 1));
        elseif strcmp(names{ii}(1:6), 'PinchM');
            raw_normalized(i, :) = raw(i, :); % Using M as reference
        elseif strcmp(names{ii}(1:6), 'PinchT');
            raw_normalized(i, :) = raw(i, :)*(GFPP(:, 2)/GFPP(:, 3));
        end
    end
    
    % Feature extraction
    mavp(:, ii) = abs(sum(raw_normalized)/fst); % Mean Absolute Value
    rmsp(:, ii) = sqrt(mean(raw_normalized.^2)); % Root-Mean-Square
    % DASDV
    for row = 2:fst
        for i = 1:2 % Both Channels
            rowcalc(row, i) = raw_normalized(row, i) - raw_normalized(row-1, i);
        end
    end
    dasdvp(ii, :) = sqrt(abs(sum(rowcalc.^2)/fst-1));
    % FFT
    rawfft = fft(raw_normalized);
    rawfftshift = fftshift(rawfft);
    rawfftshift = abs(rawfftshift(fst-5000:fst, :));
    for i = 1:2 % Both Channels
        meanfreqp(ii, i) = mean(rawfftshift(:, i));
    end
end


%%%%%%%%%%%%% Open %%%%%%%%%%%%%

% Load files (normal datasets)
names = dir('Open** e1 e2.txt');
names = {names.name};

% Pre-allocation for speed
rmso = zeros(2, length(names));
mavo = zeros(2, length(names));
actiono(1:length(names), :) = categorical(cellstr('Open')); % Labelling
rowcalc = zeros(fst, 2);
dasdvo = zeros(length(names), 2);
meanfreqo = zeros(length(names), 2);

% Segmentation and data modification
for ii = 1:length(names)
    raw = load([names{ii}]);
    raw = raw(1:fst, :);
    actiono(ii, :) = categorical(cellstr('Open')); % Labelling
    % Feature extraction (max GFP not possible for Open as it is not a grip)
    mavo(:, ii) = abs(sum(raw)/fst); % Mean Absolute Value
    rmso(:, ii) = sqrt(mean(raw.^2)); % Root-Mean-Square
    % DASDV
    for row = 2:fst
        for i = 1:2 % Both Channels
            rowcalc(row, i) = raw(row, i) - raw(row-1, i);
        end
    end
    dasdvo(ii, :) = sqrt(abs(sum(rowcalc.^2)/fst-1));
    % FFT
    rawfft = fft(raw);
    rawfftshift = fftshift(rawfft);
    rawfftshift = abs(rawfftshift(fst-5000:fst, :));
    for i = 1:2 % Both Channels
        meanfreqo(ii, i) = mean(rawfftshift(:, i));
    end
end

% Table Creation
rms = vertcat(rmsc', rmsk', rmsp', rmso');
mav = vertcat(mavc', mavk', mavp', mavo');
action = vertcat(actionc, actionk, actionp, actiono);
dasdv = vertcat(dasdvc, dasdvk, dasdvp, dasdvo);
meanfreq = vertcat(meanfreqc, meanfreqk, meanfreqp, meanfreqo);
Final = table(action, rms(:, 1), rms(:, 2), mav(:, 1), mav(:, 2), dasdv(:, 1), dasdv(:, 2), meanfreq(:, 1), meanfreq(:, 2), ...
        'VariableNames',{'Action' 'RMSA' 'RMSB', 'MAVA', 'MAVB', 'DASDVA', 'DASDVB', 'MeanFreqA', 'MeanFreqB'});
    
Final