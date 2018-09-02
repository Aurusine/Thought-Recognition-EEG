clear; 

%%%%% To change from e#1 e#3 to e#2 e#3 or etc, ctrl-f 'e#1 e#3' to 'e#2 e#3'.

fs = 256; % sampling rate
fst = fs*20; % total samples

% Pre-allocation for speed (static variables)
rowcalc = zeros(fst, 2);
meanfreqtemp = zeros(20, 2);
counter = 1;
fftbin = {zeros(1, 20)};
fftmean = zeros(44, 40); % 4 tasks, 11 datasets over 3 people
raw_normalized = zeros((fst-fs), 2); % 9 seconds

%%%%%%%%%%%%% Cylindrical Grip %%%%%%%%%%%%%

% Search directory (normal datasets)
names = dir('Cylin** e1 e2.txt');
names = {names.name};

% Search directory (maximum grip force datasets)
maxnames = dir('Max*Cylin e1 e2.txt');
maxnames = {maxnames.name};

% Pre-allocation for speed
rmsc_raw = zeros(2, length(names));
rmsc = zeros(2, length(names));
actionc(1:length(names), :) = categorical(cellstr('Cylindrical')); % Labelling
dasdvc = zeros(length(names), 2);
rmsc_max = zeros(2, length(maxnames));
GFPC = zeros(2, length(names));
GFPCUser = zeros(2, length(maxnames));
meanfreqc = zeros(length(names), 40);

% Load normal datasets, segmentation and data modification
for i = 1:length(names)
    raw = load([names{i}]);
    raw = raw(fs+1:fst, :); % Take datasets from 1-9 seconds
    rmsc_raw(:, i) = sqrt(mean(raw.^2)); % Absolute RMS values
    
    % Load maximum grip force datasets
    for ii = 1:length(maxnames)
        maxraw = load([maxnames{ii}]);
        maxraw = maxraw(1:fst, :);
        rmsc_max(:, ii) = sqrt(mean(maxraw.^2)); % Row 1 = J, Row 2 = M, Row 3 = T
        % Grip-Force Percentage calculation
        if strcmp(maxnames{ii}(1:9), 'MaxJCylin') && strcmp(names{i}(1:6), 'CylinJ');
            GFPC(:, i) = rmsc_raw(:, i)./rmsc_max(:, 1); % GFPC(:, 1) = J
        elseif strcmp(maxnames{ii}(1:9), 'MaxMCylin') && strcmp(names{i}(1:6), 'CylinM');
            GFPC(:, i) = rmsc_raw(:, i)./rmsc_max(:, 2); % GFPC(:, 2) = M
        elseif strcmp(maxnames{ii}(1:9), 'MaxTCylin') && strcmp(names{i}(1:6), 'CylinT');
            GFPC(:, i) = rmsc_raw(:, i)./rmsc_max(:, 3); % GFPC(:, 3) = T
        end
    end
end

% Calculation of reference grip force
for ii = 1:3
    for i = 1:2
        if ii == 1;
            GFPCUser(i, ii) = mean(GFPC(i, 1:4));
        elseif ii == 2;
            GFPCUser(i, ii) = mean(GFPC(i, 5:7));
        elseif ii == 3;
            GFPCUser(i, ii) = mean(GFPC(i, 8:11));
        end
    end
end

% Normalization of other subject's values based on reference GFP
for ii = 1:length(names)
    raw = load([names{ii}]);
    raw = raw(fs+1:fst, :); % Take datasets from 1-9 seconds
    for i = 1:(fst-fs)
        for iii = 1:2
            if strcmp(names{ii}(1:6), 'CylinJ');
                raw_normalized(i, iii) = raw(i, iii)*(GFPCUser(iii, 3)/GFPCUser(iii, 1));
            elseif strcmp(names{ii}(1:6), 'CylinM');
                raw_normalized(i, iii) = raw(i, iii)*(GFPCUser(iii, 3)/GFPCUser(iii, 2));
            elseif strcmp(names{ii}(1:6), 'CylinT');
                raw_normalized(i, iii) = raw(i, iii); % Using T as reference
            end
        end
    end
    
    % Feature extraction
    % Root-Mean-Square
    rmsc(:, ii) = sqrt(mean(raw_normalized.^2));

    % DASDV
    for row = 2:(fst-fs)
        for i = 1:2 % Both Channels
            rowcalc(row, i) = raw_normalized(row, i) - raw_normalized(row-1, i);
        end
    end
    dasdvc(ii, :) = sqrt(abs(sum(rowcalc.^2)/fst-1));

    % FFT
    rawfft = fft(raw_normalized);
    rawfftshift = fftshift(rawfft);
    rawfftshift = abs(rawfftshift(length(rawfftshift)/2:length(rawfftshift), :));
    for bin = 1:20 % 20 bins
        fftbin{bin} = rawfftshift(((250*(bin-1))+1):(250*bin), :); % 25 frequencies = 1 bin
        for i = 1:2 % Both Channels
            if i == 1
                meanfreqtemp(bin, i) = mean(fftbin{bin}(:, i)); % 20 features/columns each for 2 rows/channels
                meanfreqc(:, bin) = horzcat(meanfreqtemp(bin, i)); % Segment into Channel A
            elseif i == 2
                meanfreqtemp(bin, i) = mean(fftbin{bin}(:, i)); % 20 features/columns each for 2 rows/channels
                meanfreqc(:, bin+20) = horzcat(meanfreqtemp(bin, i)); % Segment into Channel B
            end
        end
    end
    fftmean(counter, :) = meanfreqc(i, :); % 40 features/columns, 11 datasets/rows
    counter = counter + 1;
end

%%%%%%%%%%%%% Key Grip %%%%%%%%%%%%%

% Search directory (normal datasets)
names = dir('Key** e1 e2.txt');
names = {names.name};

% Search directory (maximum grip force datasets)
maxnames = dir('Max*Key e1 e2.txt');
maxnames = {maxnames.name};

% Pre-allocation for speed
rmsk_raw = zeros(2, length(names));
rmsk = zeros(2, length(names));
actionk(1:length(names), :) = categorical(cellstr('Key')); % Labelling
dasdvk = zeros(length(names), 2);
rmsk_max = zeros(2, length(maxnames));
GFPK = zeros(2, length(names));
GFPKUser = zeros(2, length(maxnames));
meanfreqk = zeros(length(names), 40);

% Load normal datasets, segmentation and data modification
for i = 1:length(names)
    raw = load([names{i}]);
    raw = raw(fs+1:fst, :); % Take datasets from 1-9 seconds
    rmsk_raw(:, i) = sqrt(mean(raw.^2)); % Absolute RMS values
    
    % Load maximum grip force datasets
    for ii = 1:length(maxnames)
        maxraw = load([maxnames{ii}]);
        maxraw = maxraw(1:fst, :);
        rmsk_max(:, ii) = sqrt(mean(maxraw.^2)); % Row 1 = J, Row 2 = M, Row 3 = T
        % Grip-Force Percentage calculation
        if strcmp(maxnames{ii}(1:7), 'MaxJKey') && strcmp(names{i}(1:4), 'KeyJ');
            GFPK(:, i) = rmsk_raw(:, i)./rmsk_max(:, 1); % GFPK(:, 1) = J
        elseif strcmp(maxnames{ii}(1:7), 'MaxMKey') && strcmp(names{i}(1:4), 'KeyM');
            GFPK(:, i) = rmsk_raw(:, i)./rmsk_max(:, 2); % GFPK(:, 2) = M
        elseif strcmp(maxnames{ii}(1:7), 'MaxTKey') && strcmp(names{i}(1:4), 'KeyT');
            GFPK(:, i) = rmsk_raw(:, i)./rmsk_max(:, 3); % GFPK(:, 3) = T
        end
    end
end

% Calculation of reference grip force
for ii = 1:3
    for i = 1:2
        if ii == 1;
            GFPKUser(i, ii) = mean(GFPK(i, 1:4));
        elseif ii == 2;
            GFPKUser(i, ii) = mean(GFPK(i, 5:7));
        elseif ii == 3;
            GFPKUser(i, ii) = mean(GFPK(i, 8:11));
        end
    end
end

% Normalization of other subject's values based on reference GFP
for ii = 1:length(names)
    raw = load([names{ii}]);
    raw = raw(fs+1:fst, :); % Take datasets from 1-9 seconds
    for i = 1:(fst-fs)
        for iii = 1:2
            if strcmp(names{ii}(1:4), 'KeyJ');
                raw_normalized(i, iii) = raw(i, iii)*(GFPKUser(iii, 3)/GFPKUser(iii, 1));
            elseif strcmp(names{ii}(1:4), 'KeyM');
                raw_normalized(i, iii) = raw(i, iii)*(GFPKUser(iii, 3)/GFPKUser(iii, 2));
            elseif strcmp(names{ii}(1:4), 'KeyT');
                raw_normalized(i, iii) = raw(i, iii); % Using T as reference
            end
        end
    end
    
    % Feature extraction
    % Root-Mean-Square
    rmsk(:, ii) = sqrt(mean(raw_normalized.^2));

    % DASDV
    for row = 2:(fst-fs)
        for i = 1:2 % Both Channels
            rowcalc(row, i) = raw_normalized(row, i) - raw_normalized(row-1, i);
        end
    end
    dasdvk(ii, :) = sqrt(abs(sum(rowcalc.^2)/fst-1));

    % FFT
    rawfft = fft(raw_normalized);
    rawfftshift = fftshift(rawfft);
    rawfftshift = abs(rawfftshift(length(rawfftshift)/2:length(rawfftshift), :));
    for bin = 1:20 % 20 bins
        fftbin{bin} = rawfftshift(((250*(bin-1))+1):(250*bin), :); % 25 frequencies = 1 bin
        for i = 1:2 % Both Channels
            if i == 1
                meanfreqtemp(bin, i) = mean(fftbin{bin}(:, i)); % 20 features/columns each for 2 rows/channels
                meanfreqk(:, bin) = horzcat(meanfreqtemp(bin, i)); % Segment into Channel A
            elseif i == 2
                meanfreqtemp(bin, i) = mean(fftbin{bin}(:, i)); % 20 features/columns each for 2 rows/channels
                meanfreqk(:, bin+20) = horzcat(meanfreqtemp(bin, i)); % Segment into Channel B
            end
        end
    end
    fftmean(counter, :) = meanfreqk(i, :); % 40 features/columns, 11 datasets/rows
    counter = counter + 1;
end

%%%%%%%%%%%%% Pinch Grip %%%%%%%%%%%%%

% Search directory (normal datasets)
names = dir('Pinch** e1 e2.txt');
names = {names.name};

% Search directory (maximum grip force datasets)
maxnames = dir('Max*Pinch e1 e2.txt');
maxnames = {maxnames.name};

% Pre-allocation for speed
rmsp_raw = zeros(2, length(names));
rmsp = zeros(2, length(names));
actionp(1:length(names), :) = categorical(cellstr('Pinch')); % Labelling
dasdvp = zeros(length(names), 2);
rmsp_max = zeros(2, length(maxnames));
GFPP = zeros(2, length(names));
GFPPUser = zeros(2, length(maxnames));
meanfreqp = zeros(length(names), 40);

% Load normal datasets, segmentation and data modification
for i = 1:length(names)
    raw = load([names{i}]);
    raw = raw(fs+1:fst, :); % Take datasets from 1-9 seconds
    rmsp_raw(:, i) = sqrt(mean(raw.^2)); % Absolute RMS values
    
    % Load maximum grip force datasets
    for ii = 1:length(maxnames)
        maxraw = load([maxnames{ii}]);
        maxraw = maxraw(1:fst, :);
        rmsp_max(:, ii) = sqrt(mean(maxraw.^2)); % Row 1 = J, Row 2 = M, Row 3 = T
        % Grip-Force Percentage calculation
        if strcmp(maxnames{ii}(1:9), 'MaxJPinch') && strcmp(names{i}(1:6), 'PinchJ');
            GFPP(:, i) = rmsp_raw(:, i)./rmsp_max(:, 1); % GFPP(:, 1) = J
        elseif strcmp(maxnames{ii}(1:9), 'MaxMPinch') && strcmp(names{i}(1:6), 'PinchM');
            GFPP(:, i) = rmsp_raw(:, i)./rmsp_max(:, 2); % GFPP(:, 2) = M
        elseif strcmp(maxnames{ii}(1:9), 'MaxTPinch') && strcmp(names{i}(1:6), 'PinchT');
            GFPP(:, i) = rmsp_raw(:, i)./rmsp_max(:, 3); % GFPP(:, 3) = T
        end
    end
end

% Calculation of reference grip force
for ii = 1:3
    for i = 1:2
        if ii == 1;
            GFPPUser(i, ii) = mean(GFPP(i, 1:4));
        elseif ii == 2;
            GFPPUser(i, ii) = mean(GFPP(i, 5:7));
        elseif ii == 3;
            GFPPUser(i, ii) = mean(GFPP(i, 8:11));
        end
    end
end

% Normalization of other subject's values based on reference GFP
for ii = 1:length(names)
    raw = load([names{ii}]);
    raw = raw(fs+1:fst, :); % Take datasets from 1-9 seconds
    for i = 1:(fst-fs)
        for iii = 1:2
            if strcmp(names{ii}(1:6), 'PinchJ');
                raw_normalized(i, iii) = raw(i, iii)*(GFPPUser(iii, 3)/GFPPUser(iii, 1));
            elseif strcmp(names{ii}(1:6), 'PinchM');
                raw_normalized(i, iii) = raw(i, iii)*(GFPPUser(iii, 3)/GFPPUser(iii, 2));
            elseif strcmp(names{ii}(1:6), 'PinchT');
                raw_normalized(i, iii) = raw(i, iii); % Using T as reference
            end
        end
    end
    
    % Feature extraction
    % Root-Mean-Square
    rmsp(:, ii) = sqrt(mean(raw_normalized.^2));

    % DASDV
    for row = 2:(fst-fs)
        for i = 1:2 % Both Channels
            rowcalc(row, i) = raw_normalized(row, i) - raw_normalized(row-1, i);
        end
    end
    dasdvp(ii, :) = sqrt(abs(sum(rowcalc.^2)/fst-1));

    % FFT
    rawfft = fft(raw_normalized);
    rawfftshift = fftshift(rawfft);
    rawfftshift = abs(rawfftshift(length(rawfftshift)/2:length(rawfftshift), :));
    for bin = 1:20 % 20 bins
        fftbin{bin} = rawfftshift(((250*(bin-1))+1):(250*bin), :); % 25 frequencies = 1 bin
        for i = 1:2 % Both Channels
            if i == 1
                meanfreqtemp(bin, i) = mean(fftbin{bin}(:, i)); % 20 features/columns each for 2 rows/channels
                meanfreqp(:, bin) = horzcat(meanfreqtemp(bin, i)); % Segment into Channel A
            elseif i == 2
                meanfreqtemp(bin, i) = mean(fftbin{bin}(:, i)); % 20 features/columns each for 2 rows/channels
                meanfreqp(:, bin+20) = horzcat(meanfreqtemp(bin, i)); % Segment into Channel B
            end
        end
    end
    fftmean(counter, :) = meanfreqp(i, :); % 40 features/columns, 11 datasets/rows
    counter = counter + 1;
end

%%%%%%%%%%%%% Open %%%%%%%%%%%%%

% Search directory (normal datasets)
names = dir('Open** e1 e2.txt');
names = {names.name};

% Pre-allocation for speed
rmso = zeros(2, length(names));
actiono(1:length(names), :) = categorical(cellstr('Open')); % Labelling
dasdvo = zeros(length(names), 2);
meanfreqo = zeros(length(names), 40);

% Load normal datasets, segmentation and data modification
for i = 1:length(names)
    raw = load([names{i}]);
    raw = raw(1:fst, :);
    actiono(i, :) = categorical(cellstr('Open')); % Labelling
    
    % Feature extraction (max GFP not possible for Open as it is not a grip)
    % Root-Mean-Square
    rmso(:, i) = sqrt(mean(raw.^2));
    
    % DASDV
    for row = 2:fst
        for ii = 1:2 % Both Channels
            rowcalc(row, ii) = raw(row, ii) - raw(row-1, ii);
        end
    end
    dasdvo(i, :) = sqrt(abs(sum(rowcalc.^2)/fst-1));
    
    % FFT
    rawfft = fft(raw);
    rawfftshift = fftshift(rawfft);
    rawfftshift = abs(rawfftshift(length(rawfftshift)/2:length(rawfftshift), :));
    for bin = 1:20 % 20 bins
        fftbin{bin} = rawfftshift(((250*(bin-1))+1):(250*bin), :); % 25 frequencies = 1 bin
        for ii = 1:2 % Both Channels
            if ii == 1
                meanfreqtemp(bin, ii) = mean(fftbin{bin}(:, ii)); % 20 features/columns each for 2 rows/channels
                meanfreqo(:, bin) = horzcat(meanfreqtemp(bin, ii)); % Segment into Channel A
            elseif ii == 2
                meanfreqtemp(bin, ii) = mean(fftbin{bin}(:, ii)); % 20 features/columns each for 2 rows/channels
                meanfreqo(:, bin+20) = horzcat(meanfreqtemp(bin, ii)); % Segment into Channel B
            end
        end
    end
    fftmean(counter, :) = meanfreqo(i, :); % 40 features/columns, 11 datasets/rows
    counter = counter + 1;
end

% Table Generation
rms = vertcat(rmsc', rmsk', rmsp', rmso');
action = vertcat(actionc, actionk, actionp, actiono);
dasdv = vertcat(dasdvc, dasdvk, dasdvp, dasdvo);
MAVFFT = array2table(fftmean);
OtherFeatures = table(action, rms(:, 1), rms(:, 2), dasdv(:, 1), dasdv(:, 2),...
        'VariableNames',{'Action' 'RMSA' 'RMSB', 'DASDVA', 'DASDVB'});
    
Final = [OtherFeatures MAVFFT]

% Plotting FFT
% plot(rawfftshift(:, 1)) %, xlim([0 500])