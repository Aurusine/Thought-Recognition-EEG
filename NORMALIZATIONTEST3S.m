clear; 

fs = 2048; % sampling rate
fst = fs*10; % total samples

%%%%%%%%%%%%% Cylindrical Grip %%%%%%%%%%%%%

% Search directory (normal datasets)
names = dir('Cylin** e1 e2.txt');
names = {names.name};

% Search directory (maximum grip force datasets)
maxnames = dir('Max*Cylin e1 e2.txt');
maxnames = {maxnames.name};

% Pre-allocation for speed
rmsc_raw = zeros(2, length(names));
mavc = zeros(2, length(names));
actionc(1:length(names), :) = categorical(cellstr('Cylindrical')); % Labelling
rowcalc = zeros((fst-fs), 2);
dasdvc = zeros(length(names), 2);
rmsc_max = zeros(2, length(maxnames));
GFPC = zeros(2, length(names));
GFPCUser = zeros(2, length(maxnames));

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

% Pre-allocation for speed
raw_normalized = zeros((fst-fs), 2);
rmsc = zeros(2, length(names));

% Normalization of other subject's values based on reference GFP
for ii = 1:length(names)
    raw = load([names{ii}]);
    raw = raw(fs+1:fst, :); % Take datasets from 1-9 seconds
    for i = 1:(fst-fs)
        for iii = 1:2
            if strcmp(names{ii}(1:6), 'CylinJ');
                raw_normalized(i, iii) = raw(i, iii)*(GFPCUser(iii, 3)./GFPCUser(iii, 1));
            elseif strcmp(names{ii}(1:6), 'CylinM');
                raw_normalized(i, iii) = raw(i, iii)*(GFPCUser(iii, 3)./GFPCUser(iii, 2));
            elseif strcmp(names{ii}(1:6), 'CylinT');
                raw_normalized(i, iii) = raw(i, iii); % Using T as reference
            end
        end
    end
    
    % Feature extraction
    mavc(:, ii) = abs(sum(raw_normalized)/(fst-fs)); % Mean Absolute Value
    rmsc(:, ii) = sqrt(mean(raw_normalized.^2)); % Root-Mean-Square
    % DASDV
    for row = 2:(fst-fs)
        for i = 1:2 % Both Channels
            rowcalc(row, i) = raw_normalized(row, i) - raw_normalized(row-1, i);
        end
    end
    dasdvc(ii, :) = sqrt(abs(sum(rowcalc.^2)/fst-1));
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
mavk = zeros(2, length(names));
actionk(1:length(names), :) = categorical(cellstr('Key')); % Labelling
dasdvk = zeros(length(names), 2);
rmsk_max = zeros(2, length(maxnames));
GFPK = zeros(2, length(names));
GFPKUser = zeros(2, length(maxnames));
rmsk = zeros(2, length(names));

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
                raw_normalized(i, iii) = raw(i, iii)*(GFPKUser(iii, 3)./GFPKUser(iii, 1));
            elseif strcmp(names{ii}(1:4), 'KeyM');
                raw_normalized(i, iii) = raw(i, iii)*(GFPKUser(iii, 3)./GFPKUser(iii, 2));
            elseif strcmp(names{ii}(1:4), 'KeyT');
                raw_normalized(i, iii) = raw(i, iii); % Using T as reference
            end
        end
    end
    
    % Feature extraction
    mavk(:, ii) = abs(sum(raw_normalized)/(fst-fs)); % Mean Absolute Value
    rmsk(:, ii) = sqrt(mean(raw_normalized.^2)); % Root-Mean-Square
    % DASDV
    for row = 2:(fst-fs)
        for i = 1:2 % Both Channels
            rowcalc(row, i) = raw_normalized(row, i) - raw_normalized(row-1, i);
        end
    end
    dasdvk(ii, :) = sqrt(abs(sum(rowcalc.^2)/fst-1));
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
mavp = zeros(2, length(names));
actionp(1:length(names), :) = categorical(cellstr('Pinch')); % Labelling
dasdvp = zeros(length(names), 2);
rmsp_max = zeros(2, length(maxnames));
GFPP = zeros(2, length(names));
GFPPUser = zeros(2, length(maxnames));
rmsp = zeros(2, length(names));

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
    mavp(:, ii) = abs(sum(raw_normalized)/(fst-fs)); % Mean Absolute Value
    rmsp(:, ii) = sqrt(mean(raw_normalized.^2)); % Root-Mean-Square
    % DASDV
    for row = 2:(fst-fs)
        for i = 1:2 % Both Channels
            rowcalc(row, i) = raw_normalized(row, i) - raw_normalized(row-1, i);
        end
    end
    dasdvp(ii, :) = sqrt(abs(sum(rowcalc.^2)/fst-1));
end


%%%%%%%%%%%%% Open %%%%%%%%%%%%%

% Search directory (normal datasets)
names = dir('Open** e1 e2.txt');
names = {names.name};

% Pre-allocation for speed
rmso = zeros(2, length(names));
mavo = zeros(2, length(names));
actiono(1:length(names), :) = categorical(cellstr('Open')); % Labelling
dasdvo = zeros(length(names), 2);

% Load normal datasets, segmentation and data modification
for ii = 1:length(names)
    raw = load([names{ii}]);
    raw = raw(fs+1:fst, :); % Take datasets from 1-9 seconds
    actiono(ii, :) = categorical(cellstr('Open')); % Labelling
    % Feature extraction (max GFP not possible for Open as it is not a grip)
    mavo(:, ii) = abs(sum(raw)/(fst-fs)); % Mean Absolute Value
    rmso(:, ii) = sqrt(mean(raw.^2)); % Root-Mean-Square
    % DASDV
    for row = 2:(fst-fs)
        for i = 1:2 % Both Channels
            rowcalc(row, i) = raw(row, i) - raw(row-1, i);
        end
    end
    dasdvo(ii, :) = sqrt(abs(sum(rowcalc.^2)/fst-1));
end

% Table Creation
rms = vertcat(rmsc', rmsk', rmsp', rmso');
mav = vertcat(mavc', mavk', mavp', mavo');
action = vertcat(actionc, actionk, actionp, actiono);
dasdv = vertcat(dasdvc, dasdvk, dasdvp, dasdvo);
Final = table(action, rms(:, 1), rms(:, 2), mav(:, 1), mav(:, 2), dasdv(:, 1), dasdv(:, 2),...
        'VariableNames',{'Action' 'RMSA' 'RMSB', 'MAVA', 'MAVB', 'DASDVA', 'DASDVB'});

    
Final