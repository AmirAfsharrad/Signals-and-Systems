clc
clear
warning off

%% Q1 - Part 1
% SoundMaker.m

%% Q1 - Part 2
load('G_N_I.mat');
Fs = 200000;
m1 = SoundMaker(gain, note, interval, Fs);

%% Q1 - Part 3
sound(m1, Fs);

%% Q1 - Part 4
filename = 'music.wav';
audiowrite(filename, m1, Fs);

%% Q1 - Part 5
sound(m1, 80000);
% It still has an acceptable quality.

%% Q1 - Part 6 - 20 kHz
sound(m1, 20000);

%% Q1 - Part 6 - 100 kHz
sound(m1, 100000);

%% Q1 - Part 6 - 400 kHz
% sound(m1, 400000);
% According to the doc, Fs should be in this interval : 1000 Hz <= Fs <= 384000 Hz

%% Q1 - Part 7 - Adding DC offset to Notes

noisynote = SoundMaker(gain, note + 1, interval, Fs);
sound(noisynote, Fs);

%% Q1 - Part 7 - Adding DC offset to Intervals

noisynote = SoundMaker(gain, note , interval + 0.5, Fs);
sound(noisynote, Fs);

%% Q1 - Part 8 - Noisy Notes
s = size(gain);
noise = (randn(s(2), 1))';
noisynote = SoundMaker(gain, note + noise, interval, Fs);
sound(noisynote, Fs);

%% Q1 - Part 8 - Noisy Intervals
noisyinterval = SoundMaker(gain, note, interval + noise, Fs);
sound(noisyinterval, Fs);

%% Q1 - Part 8 - Noisy Gains
noisygain = SoundMaker(gain + noise, note, interval, Fs);
sound(noisygain, Fs);

%% Q1 - Part 9 - Noisy Signal
s = size(m1);
noise = (randn(s(2),1))'; % Gaussian noise with sigma = 1 and mean = 0
noisySignal = m1 + noise;
sound(noisySignal, Fs);

%% Q1 - Part 10 - Maximum Noise Power
s = size(m1);
noise = 20 .*(randn(s(2),1))'; 
maxnoisySignal = m1 + noise;
sound(maxnoisySignal, Fs);
P1_1 = sum(noise.*noise)./s(2)

 %% Q1 - Part 11 - P2 = P1/5
s = size(m1);
noise = 20./sqrt(5) .*(randn(s(2),1))'; 
m2 = m1 + noise;
sound(m2, Fs);
newPower = sum(noise.*noise)./s(2)
filename = 'noisyMusic.wav';
audiowrite(filename, m2, Fs);

%% Q1 - Part 12 - FFT of m1 and m2
Fs = 200000;
L = length(m1);
Y1 = fft(m1);
Y2 = fft(m2);

P2_1 = abs(Y1/L);
P1_1 = P2_1(1:L/2+1);
P1_1(2:end-1) = 2*P1_1(2:end-1);

f = Fs*(0:(L/2))/L;

figure
subplot(2,1,1);

plot(f, P1_1)
title('Single-Sided Amplitude Spectrum of m1(t)')
xlabel('f (Hz)')
ylabel('|P1_{m1}(f)|')
xlim([100,600]);


P2_2 = abs(Y2/L);
P1_2 = P2_2(1:L/2+1);
P1_2(2:end-1) = 2*P1_2(2:end-1);

subplot(2,1,2);
plot(f, P1_2)
title('Single-Sided Amplitude Spectrum of m2(t)')
xlabel('f (Hz)')
ylabel('|P1_{m2}(f)|')
xlim([100,600]);

%% Q1 - Part 16 - Low Pass Filter

M = length(m2);                   % signal length
Fs = 200000;                      % sampling rate

x = m2;

% Filter parameters:
L = M;    % filter length
fc = 445;   % cutoff frequency

hsupp = (-(L-1)/2:(L-1)/2);

hideal = (2*fc/Fs)*sinc(2*fc*hsupp/Fs);
h =  hideal; % h is our filter

% in courtesy of DSP Stack Exchange users

% Choose the next power of 2 greater than L+M-1
Nfft = 2^(ceil(log2(L+M-1))); % or 2^nextpow2(L+M-1)

x_withzeros = [ x zeros(1,Nfft-M) ];
h_withzeros = [ h zeros(1,Nfft-L) ];

figure
plotFFT(fft(h), Fs, 0, 600, 'Fourier Transform of h[n]', '|H(jw)|')

% This procedure is the same as y = conv(x,h) but takes ""far"" less time
% We have also used the conv function (See line 145)

X = fft(x_withzeros); % signal
H = fft(h_withzeros); % filter
Y = X .* H;
y = ifft(Y);

figure
plotFFT(fft(y), Fs, 0 , 600, 'Fourier Transform of the filtered version of m_2[n]', '|M_2 (jw)|');
% The same operation in time domain: (Run Time  ~  1 hr)
% y = conv(h,x)

filtered = y(1.2*10^6: 5*10^6);
sound(filtered, Fs)
filename = 'm2_f.wav';
audiowrite(filename, filtered, Fs);

%% Q1 - Part 16 - Band Pass Hamming Filter

hsupp = (-(L-1)/2:(L-1)/2);
fc = 430;
h1 = (2*fc/Fs)*sinc(2*fc*hsupp/Fs).*hamming(L)';
fc = 210;
h2 = (2*fc/Fs)*sinc(2*fc*hsupp/Fs).*hamming(L)';
h = h1 -h2;

% Choose the next power of 2 greater than L+M-1
Nfft = 2^(ceil(log2(L+M-1))); % or 2^nextpow2(L+M-1)

x_withzeros = [ x zeros(1,Nfft-M) ];
h_withzeros = [ h zeros(1,Nfft-L) ];


% This procedure is the same as y = conv(x,h) but takes ""far"" less time
% We have also used the conv function (See line 145)

X = fft(x_withzeros); % signal
H = fft(h_withzeros); % filter
Y = X .* H;
y = ifft(Y);

figure
plotFFT(Y, Fs, 0, 600, 'Fourier Transform of m2_{best}[n]', '|H(jw)|')

figure
plotFFT(H, Fs, 0, 600, 'Fourier Transform of h_{best}[n]', '|H(jw)|')

filtered_best = y(1.2*10^6: 5*10^6);
sound(filtered_best, Fs)
filename = 'm2_best.wav';
audiowrite(filename, filtered_best, Fs);

%% Q2_Part02_a
clc
n = 1 : 20;
x1 = [ones(1,6) zeros(1,5)];        % signal defenition
h1 = exp(-n(1:11)).*x1;             % signal defenition
y1_conv = conv(x1,h1);              % convolution calulation using conv
y1_MyConv = MyConv(x1,h1);          % convolution calulation using MyConv
subplot(2,2,1)                      % signal plotting
stem(x1,'x')
title('x1[n]');
subplot(2,2,2)
stem(h1,'x')
title('h1[n]');
subplot(2,2,3:4)
stem(y1_conv,'x')
hold on
stem(y1_MyConv,'x')
legend('conv result','MyConv result');
title('x1*h1[n]');
%% Q2_Part02_b
clc
clear all
syms t
x2_t = @(t) -heaviside(t+1) + 3*heaviside(t) - 2*heaviside(t-1);                % Continuous-Time signal definition
h2_t = @(t) (heaviside(t) - heaviside(t-1)*heaviside(t) - heaviside(t-1));      % Continuous-Time signal definition
Fs = 10;                                                                        % Sampling Frequency
x2 = double(subs(x2_t,t,-1.5 : 1/Fs : 3));                                      % Discrete-Time signal definition (sampling)
h2 = double(subs(h2_t,t,-1.5 : 1/Fs : 3));                                      % Discrete-Time signal definition (sampling)
y2_conv = conv(x2,h2);                                                          % convolution calulation using conv
y2_MyConv = MyConv(x2,h2);                                                      % convolution calulation using MyConv
subplot(2,2,1)                                                                  % signal plotting
stem(x2,'x')
title('x2[n]');
subplot(2,2,2)
stem(h2,'x')
title('h2[n]');
subplot(2,2,3:4)
stem(y2_conv,'x')
hold on
stem(y2_MyConv,'x')
legend('conv result','MyConv result');
title('x2*h2[n]');
%% Q2_Part02_c
clc
clear all
close all
syms t
f1 = 1;
f2 = 100;
x3_t = @(t) sin(2*pi*f1*t)+sin(2*pi*f2*t);                              % Continuous-Time signal definition
h3_t = @(t) sinc(2*f1*t)*(heaviside(t+10/f1)-heaviside(t-10/f1));       % Continuous-Time signal definition
L=1.5;                                                                  % interval margin
Fs=201;                                                                 % Sampling Frequency
x3 = double(subs(x3_t,t,-L : 1/Fs : L));                                % Discrete-Time signal definition (sampling)
h3 = double(subs(h3_t,t,-L : 1/Fs : L));                                % Discrete-Time signal definition (sampling)
y3_conv = conv(x3,h3);                                                  % convolution calulation using conv
y3_MyConv = MyConv(x3,h3);                                              % convolution calulation using MyConv
subplot(2,2,1)                                                          % signal plotting
stem(x3,'x')
title('x3[n]');
subplot(2,2,2)
stem(h3,'x')
title('h3[n]');
subplot(2,2,3:4)
stem(y3_conv,'x')
hold on
stem(y3_MyConv,'x')
legend('conv result','MyConv result');
title('x3*h3[n]');
%% Q2_Part02_d
n = -20 : 1 : 20;
x4 = zeros(1,21);               % signal definition
x4(7) = -1;                     % signal definition
x4(11) = 1;                     % signal definition
x4(15) = 2;                     % signal definition
h4 = zeros(1,21);               % signal definition
h4(11) = 1;                     % signal definition
y4_conv = conv(x4,h4);          % convolution calulation using conv
y4_MyConv = MyConv(x4,h4);      % convolution calulation using MyConv
y4_MyConv(2)=2;
subplot(2,2,1)                  % signal plotting
stem(n(11:31),x4,'x')
title('x4[n]');
subplot(2,2,2)
stem(n(11:31),h4,'x')
title('h4[n]');
subplot(2,2,3:4)
stem(n,y4_conv,'x')
hold on
stem(n,y4_MyConv,'x')
legend('conv result','MyConv result');
title('x4*h4[n]');