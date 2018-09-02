clear;
cubel=load('datasets\open, close1.txt');
cubel=cubel(1:2560,:);

% Sampling Frequency
Fs=256;
dalt=1/Fs; % Partition 'Fs' samples to 1 second

% Bandpass Filter
Fn = Fs/2;                                              % Nyquist Frequency (Hz)
Wp = [0.51 64]/Fn;                                      % Passband Frequencies (Normalised)
Ws = [0.5 64.01]/Fn;                                    % Stopband Frequencies (Normalised)
Rp = 5;                                                 % Passband Ripple (dB)
Rs = 30;                                                % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                         % Filter Order
[z,p,k] = cheby2(n,Rs,Ws);                              % Filter Design
[sosbp,gbp] = zp2sos(z,p,k);                            % Convert To Second-Order-Section For Stability
s_filt = filtfilt(sosbp,gbp, cubel);                    % Filter Signal

% Select the channels from the signal
for k = 1:2
     v = genvarname('Ch', who); % Define unique variable names
     eval([v '=s_filt(:,k)']); % Assign loaded data to variables
end

t=(0:length(Ch)-1)*dalt;
figure(1);
subplot(211);
plot(s_filt);
title('Filtered Signal');
legend('Ch1', 'Ch2')

fs=fft(Ch,256);        % Fast Fourier Transform
pp=fs.*conj(fs)/256;   % Calculation of power spectrum
ff=(0:127)/256/dalt;   % Calculate the frequency value of each point
subplot(212);
plot(ff,pp(1:128));
ylabel('Power spectral density');
xlabel('Frequency'); xlim([0 128])
title('Signal Power Spectrum');
% Colour graph
noisecol=[64 64 128 128];
gammacol=[32 32 64 64];
betacol=[16 16 32 32];
alphacol=[8 8 16 16];
thetacol=[4 4 8 8];
deltacol=[0 0 4 4];
totalcol=[0 1000 1000 0];
patch(noisecol,totalcol,[0 1 0],'FaceAlpha',0.2,'EdgeColor','none')
patch(gammacol,totalcol,[1 0 0],'FaceAlpha',0.2,'EdgeColor','none')
patch(betacol,totalcol,[1 1 0],'FaceAlpha',0.2,'EdgeColor','none')
patch(alphacol,totalcol,[0 0 1],'FaceAlpha',0.2,'EdgeColor','none')
patch(thetacol,totalcol,[0 1 1],'FaceAlpha',0.2,'EdgeColor','none')
patch(deltacol,totalcol,[1 0 1],'FaceAlpha',0.2,'EdgeColor','none')

% Define the wavelet family and the level for decomposition
 waveletFunction='db4';
 level=5;

% Multilevel 1-D wavelet decomposition
 [C,L]=wavedec(Ch,level,waveletFunction);
 [D,L]=wavedec(Ch1,level,waveletFunction);

for i = 1:2
    if i == 1
    % Reconstruct single branch from 1-D wavelet coefficients for ch0
        Gamma(:, i) = wrcoef('d',C,L,waveletFunction,2); %GAMMA, 32 - 64
        Beta(:, i) = wrcoef('d',C,L,waveletFunction,3); %BETA, 16 - 32
        Alpha(:, i) = wrcoef('d',C,L,waveletFunction,4); %ALPHA, 8 - 16
        Theta(:, i) = wrcoef('d',C,L,waveletFunction,5); %THETA, 4 - 8
        Delta(:, i) = wrcoef('a',C,L,waveletFunction,5); %DELTA, 0 - 4
    elseif i == 2
    % Reconstruct single branch from 1-D wavelet coefficients for ch1
        Gamma(:, i) = wrcoef('d',D,L,waveletFunction,2); %GAMMA, 32 - 64
        Beta(:, i) = wrcoef('d',D,L,waveletFunction,3); %BETA, 16 - 32
        Alpha(:, i) = wrcoef('d',D,L,waveletFunction,4); %ALPHA, 8 - 16
        Theta(:, i) = wrcoef('d',D,L,waveletFunction,5); %THETA, 4 - 8
        Delta(:, i) = wrcoef('a',D,L,waveletFunction,5); %DELTA, 0 - 4
    end
end

figure('Name', 'Channel 1');
subplot(5,1,1), plot(Gamma(:, 1)), ylabel('gamma'), xlim([0 2560]), ylim([-40 40])
subplot(5,1,2), plot(Beta(:, 1)), ylabel('beta'), xlim([0 2560]), ylim([-40 40])
subplot(5,1,3), plot(Alpha(:, 1)), ylabel('alpha'), xlim([0 2560]), ylim([-40 40])
subplot(5,1,4), plot(Theta(:, 1)), ylabel('theta'), xlim([0 2560]), ylim([-20 20])
subplot(5,1,5), plot(Delta(:, 1)), ylabel('delta'), xlim([0 2560]), ylim([-20 20])

figure('Name', 'Channel 2');
subplot(5,1,1), plot(Gamma(:, 2)), ylabel('gamma'), xlim([0 2560]), ylim([-40 40])
subplot(5,1,2), plot(Beta(:, 2)), ylabel('beta'), xlim([0 2560]), ylim([-40 40])
subplot(5,1,3), plot(Alpha(:, 2)), ylabel('alpha'), xlim([0 2560]), ylim([-40 40])
subplot(5,1,4), plot(Theta(:, 2)), ylabel('theta'), xlim([0 2560]), ylim([-20 20])
subplot(5,1,5), plot(Delta(:, 2)), ylabel('delta'), xlim([0 2560]), ylim([-20 20])

EEGfft = fft(Theta(:, 1));
EEGfftshift = fftshift(EEGfft);
EEGfftshift = abs(EEGfftshift(length(EEGfftshift)/2:length(EEGfftshift), :));
figure(4);
plot(EEGfftshift(:, 1)), xlim([0 500])