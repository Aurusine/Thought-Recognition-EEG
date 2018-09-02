clear; 

fs = 2048; % sampling rate
fst = fs*10; % total samples

% Pre-allocation for speed (static variables)
rowcalc = zeros(fst, 2);
meanfreqtemp = zeros(20, 2);
counter = 1;
fftbin = {zeros(1, 20)};
fftmean = zeros(44, 40); % 4 tasks, 11 datasets over 3 people

%%%%%%%%%%%%% Pinch %%%%%%%%%%%%%

% Load files (normal datasets)
names = dir('Pinch** e1 e2.txt');
names = {names.name};

% Pre-allocation for speed
rmso = zeros(2, length(names));
mavo = zeros(2, length(names));
actiono(1:length(names), :) = categorical(cellstr('Pinch')); % Labelling
dasdvo = zeros(length(names), 2);
meanfreqo = zeros(length(names), 40);

% Segmentation and data modification
for ii = 1:length(names)
    raw = load([names{ii}]);
    raw = raw(1:fst, :);
    actiono(ii, :) = categorical(cellstr('Pinch')); % Labelling
    % Feature extraction (max GFP not possible for Pinch as it is not a grip)
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
    rawfftshift = abs(rawfftshift(length(rawfftshift)/2:length(rawfftshift), :));
    for bin = 1:20 % 20 bins
        fftbin{bin} = rawfftshift(((250*(bin-1))+1):(250*bin), :); % 25 frequencies = 1 bin
        for i = 1:2 % Both Channels
            if i == 1
                meanfreqtemp(bin, i) = mean(fftbin{bin}(:, i)); % 20 features/columns each for 2 rows/channels
                meanfreqo(:, bin) = horzcat(meanfreqtemp(bin, i)); % Segment into Channel A
            elseif i == 2
                meanfreqtemp(bin, i) = mean(fftbin{bin}(:, i)); % 20 features/columns each for 2 rows/channels
                meanfreqo(:, bin+20) = horzcat(meanfreqtemp(bin, i)); % Segment into Channel B
            end
        end
    end
    fftmean(counter, :) = meanfreqo(ii, :); % 40 features/columns, 11 datasets/rows
    counter = counter + 1;
end

% Plotting FFT
fshift = (1:10241)*(fs/fst);
plot(fshift, rawfftshift(:, 1)) , xlim([0 500])