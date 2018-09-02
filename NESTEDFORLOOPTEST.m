clear; 

fs = 2048; % sampling rate
fst = fs*10; % total samples

%%%%%%%%%%%%% Cylindrical %%%%%%%%%%%%%

% Load files (normal datasets)
names = dir('Cylin** e1 e2.txt');
names = {names.name};

% Pre-allocation for speed
rmsc = zeros(2, length(names));
mavc = zeros(2, length(names));
actionc(1:length(names), :) = categorical(cellstr('Cylindrical')); % Labelling
rowcalc = zeros(fst, 2);
dasdvc = zeros(length(names), 2);

% Segmentation and data modification
for ii = 1:length(names)
    raw = load([names{ii}]);
    raw = raw(1:fst, :);
    actionc(ii, :) = categorical(cellstr('Cylindrical')); % Labelling
    % Feature extraction
    mavc(:, ii) = abs(sum(raw)/fst); % Mean Absolute Value
    rmsc(:, ii) = sqrt(mean(raw.^2)); % Root-Mean-Square
    % DASDV
    for row = 2:fst
        for i = 1:2 % Both Channels
            rowcalc(row, i) = raw(row, i) - raw(row-1, i);
        end
    end
    dasdvc(ii, :) = sqrt(abs(sum(rowcalc.^2)/fst-1));
end

%%%%%%%%%%%%% Key %%%%%%%%%%%%%

% Load files (normal datasets)
names = dir('Key** e1 e2.txt');
names = {names.name};

% Pre-allocation for speed
rmsk = zeros(2, length(names));
mavk = zeros(2, length(names));
actionk(1:length(names), :) = categorical(cellstr('Key')); % Labelling
rowcalc = zeros(fst, 2);
dasdvk = zeros(length(names), 2);

% Segmentation and data modification
for ii = 1:length(names)
    raw = load([names{ii}]);
    raw = raw(1:fst, :);
    actionk(ii, :) = categorical(cellstr('Key')); % Labelling
    % Feature extraction
    mavk(:, ii) = abs(sum(raw)/fst); % Mean Absolute Value
    rmsk(:, ii) = sqrt(mean(raw.^2)); % Root-Mean-Square
    % DASDV
    for row = 2:fst
        for i = 1:2 % Both Channels
            rowcalc(row, i) = raw(row, i) - raw(row-1, i);
        end
    end
    dasdvk(ii, :) = sqrt(abs(sum(rowcalc.^2)/fst-1));
end

%%%%%%%%%%%%% Pinch %%%%%%%%%%%%%

% Load files (normal datasets)
names = dir('Pinch** e1 e2.txt');
names = {names.name};

% Pre-allocation for speed
rmsp = zeros(2, length(names));
mavp = zeros(2, length(names));
actionp(1:length(names), :) = categorical(cellstr('Pinch')); % Labelling
rowcalc = zeros(fst, 2);
dasdvp = zeros(length(names), 2);

% Segmentation and data modification
for ii = 1:length(names)
    raw = load([names{ii}]);
    raw = raw(1:fst, :);
    actionp(ii, :) = categorical(cellstr('Pinch')); % Labelling
    % Feature extraction
    mavp(:, ii) = abs(sum(raw)/fst); % Mean Absolute Value
    rmsp(:, ii) = sqrt(mean(raw.^2)); % Root-Mean-Square
    % DASDV
    for row = 2:fst
        for i = 1:2 % Both Channels
            rowcalc(row, i) = raw(row, i) - raw(row-1, i);
        end
    end
    dasdvp(ii, :) = sqrt(abs(sum(rowcalc.^2)/fst-1));
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

% Segmentation and data modification
for ii = 1:length(names)
    raw = load([names{ii}]);
    raw = raw(1:fst, :);
    actiono(ii, :) = categorical(cellstr('Open')); % Labelling
    % Feature extraction
    mavo(:, ii) = abs(sum(raw)/fst); % Mean Absolute Value
    rmso(:, ii) = sqrt(mean(raw.^2)); % Root-Mean-Square
    % DASDV
    for row = 2:fst
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