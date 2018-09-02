clear;
cube1=load('C:\Users\Lilac\Documents\University\Biomedical Instrumentation\123\open, close2.txt');

% Sampling Frequency
fs=256;
dalt=1/fs; % Partition 'Fs' samples to 1 second

% Select both channels from the signal
for k = 1:2;
     v = genvarname('Ch', who); % Define unique variable names
     eval([v '=cube1(:,k)']); % Assign loaded data to variables
end;

t=[0:length(Ch)-1]*dalt;
figure(1);
subplot(311);
plot(Ch);
title('Original Signal');


fs=fft(Ch,256);        % Fast Fourier Transform
pp=fs.*conj(fs)/256;   % Calculation of power spectrum
ff=(0:127)/256/dalt;   % Calculate the frequency value that each point corresponds to
subplot(312);
plot(ff,pp(1:128));
ylabel('Power spectral density');
xlabel('Frequency');
title('Signal Power Spectrum');



l=wmaxlev(Ch,'db4');     % Take db4-scale decomposition of the variable s, the maximum decomposition scale is smaller than the actual
sd=wden(Ch,'minimaxi','s','mln',real(l),'db4'); % with addition noise dbN symmetric minimax principle minimaxi Threshold 's' for the soft threshold 'mln' estimated noise level in the different layers and adjust the threshold
subplot(313);plot(sd);xlabel('Filter out noise');

[LoD,HiD,LoR,HiR]=wfilters('db4');
figure(2);
subplot(221);
stem(LoD);
title('LoD Low-pass decomposition filter');
grid;
subplot(222);
stem(HiD);
title('HiD High-pass decomposition filter');
grid;
subplot(223);
stem(LoR);
title('LoR Low-pass reconstruction filter');
grid;
subplot(224);
stem(HiR);title('HiR High-pass reconstruction filter');
grid

% Define the wavelet family and the level for decomposition
 waveletFunction='db4';
 level=6;
 
% Multilevel 1-D wavelet decomposition for ch0
[C,L]=wavedec(Ch,level,waveletFunction);

% 1-D detail coefficients for ch0
C6=appcoef(C,L,waveletFunction,6);%Scale 64
D6=detcoef(C,L,6);%Scale 64
D5=detcoef(C,L,5);%Scale 32
D4=detcoef(C,L,4);%Scale 16
D3=detcoef(C,L,3);%Scale 8
D2=detcoef(C,L,2);%Scale 4
D1=detcoef(C,L,1);%Scale 2

figure(3);
subplot(611);plot(D6);ylabel('D6');
subplot(612);plot(D5);ylabel('D5');
subplot(613);plot(D4);ylabel('D4');
subplot(614);plot(D3);ylabel('D3');
subplot(615);plot(D2);ylabel('D2');
subplot(616);plot(D1);ylabel('D1');
title('Detail coefficient');

[C,L]=wavedec(Ch,6,waveletFunction);%Scale 64
C6=appcoef(C,L,waveletFunction,6);
[C,L]=wavedec(Ch,5,waveletFunction);%Scale 32
C5=appcoef(C,L,waveletFunction,5);
[C,L]=wavedec(Ch,4,waveletFunction);%Scale 16
C4=appcoef(C,L,waveletFunction,4);
[C,L]=wavedec(Ch,3,waveletFunction);%Scale 8
C3=appcoef(C,L,waveletFunction,3);
[C,L]=wavedec(Ch,2,waveletFunction);%Scale 4
C2=appcoef(C,L,waveletFunction,2);
[C,L]=wavedec(Ch,1,waveletFunction);%Scale 2
C1=appcoef(C,L,waveletFunction,1);

figure(4);
subplot(611);plot(C6);ylabel('C6');
subplot(612);plot(C5);ylabel('C5');
subplot(613);plot(C4);ylabel('C4');
subplot(614);plot(C3);ylabel('C3');
subplot(615);plot(C2);ylabel('C2');
subplot(616);plot(C1);ylabel('C1');
title('Approximate coefficient');

%fs=256 Hz
%?-wave(0-3Hz)?-wave(4-7Hz)?-wave(8~15Hz);?-wave(16~31Hz)?-wave(32-64Hz);
%******************************
% Multilevel 1-D wavelet decomposition for ch0
[C,L]=wavedec(Ch,6,waveletFunction); %Scale 256~64

% 1-D detail coefficients for ch0
D1=detcoef(C,L,1);      %Scale 2, 128 - 256
D2=detcoef(C,L,2);      %Scale 4, 64 - 128
D3=detcoef(C,L,3);      %Scale 8, 32 - 64
D4=detcoef(C,L,4);      %Scale 16, 16 - 32
D5=detcoef(C,L,5);      %Scale 32, 8 - 16
D6=detcoef(C,L,6);      %Scale 64, 4 - 8
C6=appcoef(C,L,waveletFunction,6); %Scale 256, 0 - 4

% Reconstruct single branch from 1-D wavelet coefficients for ch0
SRC1=wrcoef('a',C,L,waveletFunction,1);
SRC2=wrcoef('a',C,L,waveletFunction,2);
SRC3=wrcoef('a',C,L,waveletFunction,3);
SRC4=wrcoef('a',C,L,waveletFunction,4);
SRC5=wrcoef('a',C,L,waveletFunction,5);
SRC6=wrcoef('a',C,L,waveletFunction,6);    % reconstruction

SRD1=wrcoef('d',C,L,waveletFunction,1);
SRD2=wrcoef('d',C,L,waveletFunction,2);
SRD3=wrcoef('d',C,L,waveletFunction,3);
SRD4=wrcoef('d',C,L,waveletFunction,4);
SRD5=wrcoef('d',C,L,waveletFunction,5);
SRD6=wrcoef('d',C,L,waveletFunction,6);

gamma=[SRD3];       %?-wave(32~64Hz)
beta=[SRD4];        %?-wave(16~31Hz)
alpha=[SRD5];       %?-wave(8~15Hz)
theta=[SRD6];       %?-wave(4~7Hz)
delta=[SRC6];       %?-wave(0~3Hz)

%delta
fs1=fft(delta,128);         %Fast fourier transform
pp1=fs1.*conj(fs1)/128;     %Calculation of power spectrum
ff1=((0:50)/128)*256;      %Calculate the frequency value of each point corresponds to

%theta
fs2=fft(theta,128);         %Fast fourier transform
pp2=fs2.*conj(fs2)/128;     %Calculation of power spectrum
ff2=((0:50)/128)*256;      %Calculate the frequency value of each point corresponds to

%alpha
fs3=fft(alpha,128);         %Fast fourier transform
pp3=fs3.*conj(fs3)/128;     %Calculation of power spectrum
ff3=((0:50)/128)*256;      %Calculate the frequency value of each point corresponds to

%beta
fs4=fft(beta,128);          %Fast fourier transform
pp4=fs4.*conj(fs4)/128;     %Calculation of power spectrum
ff4=((0:50)/128)*256;      %Calculate the frequency value of each point corresponds to

%gamma
fs5=fft(gamma,128);         %Fast fourier transform
pp5=fs5.*conj(fs5)/128;     %Calculation of power spectrum
ff5=((0:50)/128)*256;      %Calculate the frequency value of each point corresponds to

%Graph of PSD alpha and beta wave is plotted
figure (5);
plot(ff1,pp1(1:51),ff2,pp2(1:51),ff3,pp3(1:51),ff4,pp4(1:51),ff5,pp5(1:51));
xlabel('Frequency(Hz)'),
ylabel('Power spectral density µV ^2/Hz')
legend('Delta','theta','alpha','beta','gamma')

figure('Name', 'Channel 1');
subplot(5,1,1), plot(gamma), ylabel('gamma')
subplot(5,1,2), plot(beta), ylabel('beta')
subplot(5,1,3), plot(alpha), ylabel('alpha')
subplot(5,1,4), plot(theta), ylabel('theta')
subplot(5,1,5), plot(delta), ylabel('delta')