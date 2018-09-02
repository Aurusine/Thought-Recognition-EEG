clear; 

fs = 2048; % sampling rate
fst = fs*10; % total samples

%%%%%%%%%%%%% Cylindrical Grip %%%%%%%%%%%%%

% Search directory (normal datasets)
names = dir('Cylin** e1 e2.txt');
names = {names.name};

rmscounter = 1;

for i = 1:length(names)
    raw = load([names{i}]);
    raw = raw(1:fst, :); % Segment to 10 seconds
    for x = 1:3
        rawtemp{x} = raw((fs)+(3*fs*(x-1)+1):(fs+(3*fs*(x))), :); % Take 3 datasets from 1-9 seconds
    end
    for x = 1:3
        rmsc_raw(:, rmscounter) = rawtemp{1, x};
        rmscounter = rmscounter + 1;
    end
    rmsc_raw(:, i) = sqrt(mean(raw.^2)); % Absolute RMS values
end

length(rawtemp)