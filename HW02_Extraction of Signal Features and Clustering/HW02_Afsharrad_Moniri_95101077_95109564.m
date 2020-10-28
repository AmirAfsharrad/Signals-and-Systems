clc
clear all
format long

%% Q1 - Part 1
load('EEG_signal.mat');
Fs = 2400;

%% Q1 - Removing DC Component

meanMAT = mean(EEG_signal_edited, 2);
s = size(EEG_signal_edited);
EEG = zeros(s);
for i = 1 : s(1)
    for j = 1 : s(3)
        EEG(i, :, j) = EEG_signal_edited(i, :, j) - meanMAT(i, :, j);
    end
end

%% Q1 - Part 2

data = EEG(1,:,1:5);
data = reshape(data, 1, 5*7200);
L = length(data);


%% Q1 - Part 3

figure
plotFFT(data, 2400, 0, 100, 'Frequency Spectrum', 'P2', 20); % This function plots fft of the given signal


Y1 = fft(data);
L = length(Y1);
P2 = abs(Y1/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
E = 0;
total = sum(P1(1:end).^2);

for action = 1 : length(P1)
    E = E + P1(action)^2;
    cdf(action) = E/total;
    
    if (E/total > 0.9)
        break
    end
end

disp(['Threshold Frequency: ', num2str(f(action)), ' Hz']);

%% Question 1 - Section 8

% A sample filter with w0 = 0.2 and r = 1
% for other setting please check the report
% This function generates a low pass filter using a method introduced by us

LowPassFilterDesigner(1, 0.25, 90, 1/3);

%% Question 1 - Section 9

figure
LowPassFilterDesigner(0.25, 0.5, 95, 1/3);
title('$\omega_0 = \pi/2 , r_0 = 0.25$','interpreter','latex');

figure
LowPassFilterDesigner(0.5, 0.5, 95, 1/3);
title('$\omega_0 = \pi/2 , r_0 = 0.5$','interpreter','latex');

figure
LowPassFilterDesigner(0.9, 0.5, 95, 1/3);
title('$\omega_0 = \pi/2 , r_0 = 0.9$','interpreter','latex');

figure
LowPassFilterDesigner(0.95, 0.5, 95, 1/3);
title('$\omega_0 = \pi/2 , r_0 = 0.95$','interpreter','latex');


figure
LowPassFilterDesigner(1.05, 0.5, 95, 1/3);
title('$\omega_0 = \pi/2 , r_0 = 1.05$','interpreter','latex');

figure
LowPassFilterDesigner(1.1, 0.5, 95, 1/3);
title('$\omega_0 = \pi/2 , r_0 = 1.1$ - Order','interpreter','latex');

%% Question 1 - Section 10

% Truncated Sinc Filter
N = 8;
w0 = pi/2;

L = (N-1)/2 * (w0/pi);
num = sinc(-L:w0/pi:L);

figure
freqz(num, polyval(num, 1));
title('Truncated Sinc Filter - N = 8');
% Matlab's Built-in Butter filter
[b,a] = butter(N, w0/pi);
figure
freqz(b,a)
title('Matlab Built-in Butter Filter - N = 8');

%% Question 1 - Section 11

% Our Truncated Sinc Filter
w0 = 49/1200*pi;
percent = 1 - 1/1000;
N = 1;
while(1)
    L = (N-1)/2 * (w0/pi);
    num = sinc(-L:w0/pi:L);
    [h, w] = freqz(num, polyval(num, 1));
    total_energy = sum(abs(h).^2);
    E = 0;
    for i = 1 : length(w)
        E = E + abs(h(i).^2);
        if(E/total_energy > percent)
           break
        end
    end
    if(w(i)<=w0)
        break
    end
    N = N + 1;
end

figure
disp(['Filter Order ', num2str(N)]);
freqz(num, polyval(num, 1))
title('Truncated Sinc Filter');

% Filter data using our filter
y = filter(num, polyval(num, 1), data);

figure
subplot(2,1,1)
plotFFT(data, 2400, 0, 100, 'Frequency Spectrum', 'P2', 20); % This function plots fft of the given signal
subplot(2,1,2)
plotFFT(y, 2400, 0, 100, 'Frequency Spectrum', 'P2', 20); % This function plots fft of the given signal


%%
% Matlab's Built in Butterworth Filter
w0 = 49/1200*pi;
percent = 1 - 1/1000;

N = 1;
while(1)
    [b , a] = butter(N, w0/pi);
    [h, w] = freqz(b, a);
    total_energy = sum(abs(h).^2);
    E = 0;
    for i = 1 : length(w)
        E = E + abs(h(i).^2);
        if(E/total_energy > percent)
           break
        end
    end
    if(w(i)<=w0)
        break
    end
    N = N + 1;
    
end

figure
disp(['Filter Order ', num2str(N)]);
freqz(b, a)



% Filter data using matlab's built in butter function
y = filter(b, a, data);
figure
plotFFT(y, 2400, 0, 100, 'Frequency Spectrum', 'P2', 20); % This function plots fft of the given signal

%% Question 1 - Section 12

% This function creates a high pass filter using the method discussed in
% the report
figure
HighPassFilterDesigner (1, 0.75, 95, 1/3);


%% High Pass Zero-Pole Filter
figure
HighPassFilterDesigner(0.25, 0.5, 95, 1/3);
title('$\omega_0 = \pi/2 , r_0 = 0.25$','interpreter','latex');

figure
HighPassFilterDesigner(0.5, 0.5, 95, 1/3);
title('$\omega_0 = \pi/2 , r_0 = 0.5$','interpreter','latex');

figure
HighPassFilterDesigner(0.9, 0.5, 95, 1/3);
title('$\omega_0 = \pi/2 , r_0 = 0.9$','interpreter','latex');

figure
HighPassFilterDesigner(0.95, 0.5, 95, 1/3);
title('$\omega_0 = \pi/2 , r_0 = 0.95$','interpreter','latex');

figure
HighPassFilterDesigner(1.05, 0.5, 95, 1/3);
title('$\omega_0 = \pi/2 , r_0 = 1.05$','interpreter','latex');

figure
HighPassFilterDesigner(1.1, 0.5, 95, 1/3);
title('$\omega_0 = \pi/2 , r_0 = 1.1$','interpreter','latex');

%% High Pass Trnucated Sinc Filter

N = 91;
w0 = pi/2;

L = (N-1)/2 * (w0/pi);
delta = zeros(1,length(-L:w0/pi:L));
delta(ceil(N/2)) = 0.6*pi;
num = sinc(-L:w0/pi:L);

figure
freqz(delta-num, polyval(delta-num, -1));


[b, a] = butter(81, 0.5, 'high');
figure
freqz(b,a)

%% Question 1 - Part 13

N = 1377;
w0 = pi/3;
w1 = 2*pi/3;

L = (N-1)/2 * (w0/pi);
delta = zeros(1,length(-L:w0/pi:L));
delta(ceil(N/2)) = 0.99*pi;
numH = sinc(-L:w0/pi:L).*(hamming(length(-L:w0/pi:L)))';

numH = delta - numH;

figure
freqz(numH, polyval(numH, -1));


N = 181;

L = (N-1)/2 * (w1/pi);
numL = sinc(-L:w1/pi:L);

figure
freqz(numL, polyval(numL, 1));

num = conv(numH,numL.*(hamming(length(-L:w1/pi:L)))');

figure
freqz(num, polyval(numL, 1)*polyval(numH, -1));


%% Question 2 - Part 1

stats(1).std = std((EEG_signal_edited(:,:,1))');
stats(2).std = std((EEG_signal_edited(:,:,2))');
stats(3).std = std((EEG_signal_edited(:,:,3))');

stats(1).mean = mean((EEG_signal_edited(:,:,1))');
stats(2).mean = mean((EEG_signal_edited(:,:,2))');
stats(3).mean = mean((EEG_signal_edited(:,:,3))');

stats(1).max = max((EEG_signal_edited(:,:,1))');
stats(2).max = max((EEG_signal_edited(:,:,2))');
stats(3).max = max((EEG_signal_edited(:,:,3))');

stats(1).min = min((EEG_signal_edited(:,:,1))');
stats(2).min = min((EEG_signal_edited(:,:,2))');
stats(3).min = min((EEG_signal_edited(:,:,3))');


 %% Question 2 - Part 2
 for action = 1:3 % Actions
     figure
     for channel = 1:5 % Channels
        subplot(5,2,2*channel-1);
        y = EEG(channel,:,action);
        plot(0:1/Fs:7199/Fs, y); 
        title(['Channel ' , num2str(channel) ,' Action ', num2str(action)]); xlabel('time'); ylabel('Voltage')
        subplot(5,2,2*channel);
        plotFFT(y, Fs, 0, 100, ['FFT of Channel  ' , num2str(channel) ,' Action ', num2str(action)],'fft', 20);
    end
 end
 
%% Question 2 - Part 4
% removing above 60 frequencies
filtered = zeros(64,7200,49);
for action = 1:49 % Actions
    for channel = 1:64 % Channels 
        filtered(channel,:,action) = LPF(EEG(channel,:,action), 200, 60, 2400);     
    end
end 

% downsampling to Fs = 120Hz
    downFiltered = filtered(:,1:20:end,:);
    
    
for action = 1:49 % Actions
    for channel = 1:64 % Channels 
    
    FFT = fft(downFiltered(channel,:,action));
    FFT(151) = 0.5*(FFT(152) + FFT(150));
    FFT(211) = 0.5*(FFT(210) + FFT(212));

    
    downFiltered(channel, :, action) = ifft(FFT);
 
    end
end 


%%
% Plotting the filtered signal
for action = 1:3 % Actions
figure
    for channel = 1:5 % Channels
        subplot(5,2,2*channel-1);
        y = downFiltered(channel,:,action);
        plot(y); 
        title(['Channel ' , num2str(channel) ,' Action ', num2str(action)]); xlabel('time'); ylabel('Voltage')
        subplot(5,2,2*channel);
        
        plotFFT(y, Fs/20, 0, 60, ['FFT of Channel  ' , num2str(channel) ,' Action ', num2str(action)],'fft', 10);
    end
end

%% Question 2 - Part5
load('lookup.mat');
for i = 1 : length(lookup)
    figure
    subplot(3,1,1)
    plot(downFiltered(lookup(i,1),:,1),downFiltered(lookup(i,2),:,1),'.')
    P = polyfit(downFiltered(lookup(i,1),:,1),downFiltered(lookup(i,2),:,1),1);
    x = min(downFiltered(lookup(i,1),:,1)) : max(downFiltered(lookup(i,1),:,1));
    yfit = P(1)*x+P(2);
    hold on;    
    plot(x,yfit,'r-.','LineWidth',2);
    title('Action 1');
	xlabel(['Channel No. ',num2str(lookup(i,1))]); ylabel(['Channel No. ',num2str(lookup(i,2))])

    subplot(3,1,2)
    plot(downFiltered(lookup(i,1),:,2),downFiltered(lookup(i,2),:,2),'.')
    P = polyfit(downFiltered(lookup(i,1),:,2),downFiltered(lookup(i,2),:,2),1);
    x = min(downFiltered(lookup(i,1),:,2)) : max(downFiltered(lookup(i,1),:,2));
    yfit = P(1)*x+P(2);
    hold on;
    plot(x,yfit,'r-.','LineWidth',2);
    title('Action 2');
	xlabel(['Channel No. ',num2str(lookup(i,1))]); ylabel(['Channel No. ',num2str(lookup(i,2))])
    
    subplot(3,1,3)
    plot(downFiltered(lookup(i,1),:,3),downFiltered(lookup(i,2),:,3),'.')
    P = polyfit(downFiltered(lookup(i,1),:,3),downFiltered(lookup(i,2),:,3),1);
    x = min(downFiltered(lookup(i,1),:,3)) : max(downFiltered(lookup(i,1),:,3));
    yfit = P(1)*x+P(2);
    hold on;
    plot(x,yfit,'r-.','LineWidth',2);
    title('Action 3');
	xlabel(['Channel No. ',num2str(lookup(i,1))]); ylabel(['Channel No. ',num2str(lookup(i,2))])
    
end

%% Question 2 - part 6
cor_mat = corr (downFiltered(:,:,1)');
figure
bar3(abs(cor_mat));

%% Question 2 - part 7

boxplotData = zeros(16,360,3);
for i = 1 : 16
    boxplotData(i,:,1) = downFiltered(i,:,1);
    boxplotData(i,:,2) = downFiltered(i,:,2);
    boxplotData(i,:,3) = downFiltered(i,:,3);
end
figure
boxplot(boxplotData(:,:,1)', 'Whisker', 5);
title('Action 1');

figure
boxplot(boxplotData(:,:,2)', 'Whisker', 5);
title('Action 2');

figure
boxplot(boxplotData(:,:,3)', 'Whisker', 5);
title('Action 3');

%% Question 2 - Part 8

for channel = 1 : 16
    alpha_band(channel,:) = BPF(downFiltered(channel,:,1), 200, 8, 15, 2400/20)';
    beta_band(channel,:)  = BPF(downFiltered(channel,:,1), 200, 16, 31, 2400/20)';
    gamma_band(channel,:) = BPF(downFiltered(channel,:,1), 200, 32, 60, 2400/20)';
    delta_band(channel,:) = LPF(downFiltered(channel,:,1), 200, 4, 2400/20)';
end

bandsData = zeros(16,360,4);
bandsData(:,:,1) = alpha_band;
bandsData(:,:,2) = beta_band;
bandsData(:,:,3) = gamma_band;
bandsData(:,:,4) = delta_band;
save('bandsData','bandsData');
clear bandsData

for channel = 1 : 4
    figure
    %alpha
    subplot(4,2,1);
    y1 = alpha_band(channel,:);
    plot(0:20/Fs:20*359/Fs, y1); 
    title(['Channel ' , num2str(channel) ' - Alpha Band']); xlabel('time'); ylabel('Voltage');
    subplot(4,2,2);
    plotFFT(y1, Fs/20, 0, 60, ['FFT of Channel  ' , num2str(channel) ,' Action 1'],'fft', 5);
    title(['FFT of Channel ' , num2str(channel) ' - Alpha Band'])

    %beta
    subplot(4,2,3);
    y = beta_band(channel,:);
    plot(0:20/Fs:20*359/Fs, y); 
    title(['Channel ' , num2str(channel) ' - Beta Band']); xlabel('time'); ylabel('Voltage');
    subplot(4,2,4);
    plotFFT(y, Fs/20, 0, 60, ['FFT of Channel  ' , num2str(channel) ,' Action 1'],'fft', 5);
    title(['FFT of Channel ' , num2str(channel) ' - Beta Band'])
    
    %gamma
    subplot(4,2,5);
    y = gamma_band(channel,:);
    plot(0:20/Fs:20*359/Fs, y); 
    title(['Channel ' , num2str(channel) ' - Gamma Band']); xlabel('time'); ylabel('Voltage');
    subplot(4,2,6);
    plotFFT(y, Fs/20, 0, 60, ['FFT of Channel  ' , num2str(channel) ,' Action 1'],'fft', 5);
    title(['FFT of Channel ' , num2str(channel) ' - Gamma Band'])

    %delta
    subplot(4,2,7);
    y = delta_band(channel,:);
    plot(0:20/Fs:20*359/Fs, y); 
    title(['Channel ' , num2str(channel) ' - Delta Band']); xlabel('time'); ylabel('Voltage');
    subplot(4,2,8);
    plotFFT(y, Fs/20, 0, 60, ['FFT of Channel  ' , num2str(channel) ,' Action 1'],'fft', 5);
    title(['FFT of Channel ' , num2str(channel) ' - Delta Band'])
end

%% Question 2 - Part 9 - First Method
% first method
alphaEnergy1 = zeros(4,12);
betaEnergy1 = zeros(4,12);
gammaEnergy1 = zeros(4,12);
deltaEnergy1 = zeros(4,12);

for channel = 1 : 4
    for bin = 1 : 12
        alphaEnergy1(channel, bin) = sum(alpha_band(channel, (bin-1)*30 + 1 : 30*bin).^2);
        betaEnergy1(channel, bin) = sum(beta_band(channel, (bin-1)*30 + 1 : 30*bin).^2);
        gammaEnergy1(channel, bin) = sum(gamma_band(channel, (bin-1)*30 + 1 : 30*bin).^2);
        deltaEnergy1(channel, bin) = sum(delta_band(channel, (bin-1)*30 + 1 : 30*bin).^2);
    end
end

for channel = 1 : 4
    figure
    
    subplot(2,2,1)
    stem(alphaEnergy1(channel,:))
    xlabel('Bin Number');
    ylabel('E_\alpha');
    title(['Alpha Band Energy Channel ',num2str(channel)]);
    
    subplot(2,2,2)
    stem(betaEnergy1(channel,:))
    xlabel('Bin Number');
    ylabel('E_\beta');
    title(['Beta Band Energy Channel ',num2str(channel)]);
    
    subplot(2,2,3)
    stem(gammaEnergy1(channel,:))
    xlabel('Bin Number');
    ylabel('E_\gamma');
    title(['Gamma Band Energy Channel ',num2str(channel)]);
    
    subplot(2,2,4)
    stem(deltaEnergy1(channel,:))
    xlabel('Bin Number');
    ylabel('E_\delta');
    title(['Delta Band Energy Channel ',num2str(channel)]);
end


%% Question 2 - Part 9 - Second Method
alphaEnergy2 = zeros(4,360-30);
betaEnergy2 = zeros(4,360-30);
gammaEnergy2 = zeros(4,360-30);
deltaEnergy2 = zeros(4,360-30);

for channel = 1 : 4
    for bin = 1 : 360-30
        alphaEnergy2(channel, bin) = sum(alpha_band(channel, bin : bin+30).^2);
        betaEnergy2(channel, bin) = sum(beta_band(channel, bin : bin+30).^2);
        gammaEnergy2(channel, bin) = sum(gamma_band(channel, bin : bin+30).^2);
        deltaEnergy2(channel, bin) = sum(delta_band(channel, bin : bin+30).^2);
    end
end

t = 0 : 1/120 : 329/120;
for channel = 1 : 4
    figure
    
    subplot(2,2,1)
    plot(t, alphaEnergy2(channel,:))
    xlabel('t(s)');
    ylabel('E_\alpha');
    title(['Alpha Band Energy Channel ',num2str(channel)]);
    
    subplot(2,2,2)
    plot(t, betaEnergy2(channel,:))
    xlabel('t(s)');
    ylabel('E_\beta');
    title(['Beta Band Energy Channel ',num2str(channel)]);
    
    subplot(2,2,3)
    plot(t, gammaEnergy2(channel,:))
    xlabel('t(s)');
    ylabel('E_\gamma');
    title(['Gamma Band Energy Channel ',num2str(channel)]);
    
    subplot(2,2,4)
    plot(t, deltaEnergy2(channel,:))
    xlabel('t(s)');
    ylabel('E_\delta');
    title(['Delta Band Energy Channel ',num2str(channel)]);
end

%% Question 2 - Part 10
% Clustering for action one

correlation = abs(cor_mat);
table = zeros(64, 64);
thr = 0.8;

for i = 1 : 64
   for j = 1 : 64
      if(correlation(i,j) > thr)
            table(i, j) = 1;
      end
   end
end

clusters_action1 = final_cluster(correlation, 0.75, 1:64, 1000);

% action 2
correlation = abs(corr (downFiltered(:,:,2)'));
clusters_action3 = final_cluster(correlation, 0.75, 1:64, 1000);

%% Section 2 - Question 11

for channel = 1 : 64
    for action = 1 : 49
        alpha_band1(channel,action,:) = BPF(downFiltered(channel,:,action), 200, 8, 15, 2400/20)';
        beta_band1 (channel,action,:) = BPF(downFiltered(channel,:,action), 200, 16, 31, 2400/20)';
    end
end

for channel = 1 : 64
    for action = 1 : 49
        alpha_E (channel, action) = sum(alpha_band1(channel, action , :).^2);
        beta_E (channel, action) = sum(beta_band1(channel, action , :).^2);
    end
end
E_alpha = mean(alpha_E,1)';
E_beta  = mean( beta_E,1)';



R = E_alpha./E_beta;

figure
stem((R - mean(R))/std(R));

hold on
stem((E_alpha - mean(E_alpha))/std(E_alpha));

ylabel('zscore')
xlabel('action')
legend('R', 'Alpha Energy')






