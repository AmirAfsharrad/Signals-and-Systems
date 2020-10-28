%% Question 1
clear
clc
img1 = imread('img03.jpg');
%% down-sampling (the oringinal image is too high-resolution)
% img1 = img1(1:3:end,1:3:end,:);
img01 = im2double(img1);

%% Original Image
figure
imshow(img01)
title('Original Image')

%% Additive White Gaussian Noise
img01_AWGN = imnoise(img01,'gaussian',0,0.8);

figure
imshow(img01_AWGN)
title('Additive White Gaussian Noise')

%% Salt&Pepper (Impulse) Noise
img01_salt_pepper = imnoise(img01,'salt & pepper',.2);

figure
imshow(img01_salt_pepper)
title('Salt&Pepper Noise')

%% Shot (Poisson) Noise
scale = 1e10;
img01_shot = scale*imnoise(img01/scale,'poisson');

figure
imshow(img01_shot)
title('Shot (Poisson) Noise')

%% Speckle (Multiplicative) Noise
img01_speckle = imnoise(img01,'speckle',0.2);

figure
imshow(img01_speckle)
title('Speckle (Multiplicative) Noise')

%% Uniform Additive Noise
UAN = random('uniform',-0.5,0.5,[size(img01,1),size(img01,2),3]);
img01_unif_additive = img01 + UAN;

figure
imshow(img01_unif_additive)
title('Uniform Additve Noise')

%% mirroring noisy images for better filtering
img01_AWGN = mirror_edge(img01_AWGN, 20);
img01_salt_pepper = mirror_edge(img01_salt_pepper, 20);
img01_shot = mirror_edge(img01_shot, 20);
img01_speckle = mirror_edge(img01_speckle, 20);
img01_unif_additive = mirror_edge(img01_unif_additive, 20);
%% linear smoothig filter 1 - moving average
h = fspecial('average',[4 4]);

% Additive White Gaussian Noise 
array = img01_AWGN;
img_filtered(:,:,1) = imfilter(array(:,:,1),h);
img_filtered(:,:,2) = imfilter(array(:,:,2),h);
img_filtered(:,:,3) = imfilter(array(:,:,3),h);
img_filtered = remove_edge(img_filtered, 20);
figure
imshow(img_filtered)
title('Additive White Gaussian Noise Filtered with moving average')
clear img_filtered

% Salt&Pepper Noise
array = img01_salt_pepper;
img_filtered(:,:,1) = imfilter(array(:,:,1),h);
img_filtered(:,:,2) = imfilter(array(:,:,2),h);
img_filtered(:,:,3) = imfilter(array(:,:,3),h);
img_filtered = remove_edge(img_filtered, 20);
figure
imshow(img_filtered)
title('Salt&Pepper Noise Filtered with moving average')
clear img_filtered

% Shot (Poisson) Noise
array = img01_shot;
img_filtered(:,:,1) = imfilter(array(:,:,1),h);
img_filtered(:,:,2) = imfilter(array(:,:,2),h);
img_filtered(:,:,3) = imfilter(array(:,:,3),h);
img_filtered = remove_edge(img_filtered, 20);
figure
imshow(img_filtered)
title('Shot (Poisson) Noise Filtered with moving average')
clear img_filtered

% Speckle (Multiplicative) Noise
array = img01_speckle;
img_filtered(:,:,1) = imfilter(array(:,:,1),h);
img_filtered(:,:,2) = imfilter(array(:,:,2),h);
img_filtered(:,:,3) = imfilter(array(:,:,3),h);
img_filtered = remove_edge(img_filtered, 20);
figure
imshow(img_filtered)
title('Speckle (Multiplicative) Noise Filtered with moving average')
clear img_filtered
%%
% Uniform Additve Noise
array = img01_unif_additive;
img_filtered(:,:,1) = imfilter(array(:,:,1),h);
img_filtered(:,:,2) = imfilter(array(:,:,2),h);
img_filtered(:,:,3) = imfilter(array(:,:,3),h);
img_filtered = remove_edge(img_filtered, 20);
figure
imshow(img_filtered)
title('Uniform Additve Noise Filtered with moving average')
clear img_filtered

%% linear smoothig filter 2 - Gaussian Filter
sigma = 2;
% Additive White Gaussian Noise 
array = img01_AWGN;
img_filtered(:,:,1) = imgaussfilt(array(:,:,1),sigma);
img_filtered(:,:,2) = imgaussfilt(array(:,:,2),sigma);
img_filtered(:,:,3) = imgaussfilt(array(:,:,3),sigma);
img_filtered = remove_edge(img_filtered, 20);
figure
imshow(img_filtered)
title('Additive White Gaussian Noise Filtered with Gaussian filter')
clear img_filtered

% Salt&Pepper Noise
array = img01_salt_pepper;
img_filtered(:,:,1) = imgaussfilt(array(:,:,1),sigma);
img_filtered(:,:,2) = imgaussfilt(array(:,:,2),sigma);
img_filtered(:,:,3) = imgaussfilt(array(:,:,3),sigma);
img_filtered = remove_edge(img_filtered, 20);

figure
imshow(img_filtered)
title('Salt&Pepper Noise Filtered with Gaussian filter')
clear img_filtered

% Shot (Poisson) Noise
array = img01_shot;
img_filtered(:,:,1) = imgaussfilt(array(:,:,1),sigma);
img_filtered(:,:,2) = imgaussfilt(array(:,:,2),sigma);
img_filtered(:,:,3) = imgaussfilt(array(:,:,3),sigma);
img_filtered = remove_edge(img_filtered, 20);
figure
imshow(img_filtered)
title('Shot (Poisson) Noise Filtered with Gaussian filter')
clear img_filtered

% Speckle (Multiplicative) Noise
array = img01_speckle;
img_filtered(:,:,1) = imgaussfilt(array(:,:,1),sigma);
img_filtered(:,:,2) = imgaussfilt(array(:,:,2),sigma);
img_filtered(:,:,3) = imgaussfilt(array(:,:,3),sigma);
img_filtered = remove_edge(img_filtered, 20);
figure
imshow(img_filtered)
title('Speckle (Multiplicative) Noise Filtered with Gaussian filter')
clear img_filtered

% Uniform Additve Noise
array = img01_unif_additive;
img_filtered(:,:,1) = imfilter(array(:,:,1),sigma);
img_filtered(:,:,2) = imfilter(array(:,:,2),sigma);
img_filtered(:,:,3) = imfilter(array(:,:,3),sigma);
img_filtered = remove_edge(img_filtered, 20);
figure
imshow(img_filtered)
title('Uniform Additve Noise Filtered with Gaussian filter')
clear img_filtered

%% Wiener Filtering
k = 4;
% Additive White Gaussian Noise 
array = img01_AWGN;
img_filtered(:,:,1) = wiener2(array(:,:,1),[k k]);
img_filtered(:,:,2) = wiener2(array(:,:,2),[k k]);
img_filtered(:,:,3) = wiener2(array(:,:,3),[k k]);
img_filtered = remove_edge(img_filtered, 20);
figure
imshow(img_filtered)
title('Additive White Gaussian Noise Filtered with Wiener filter')
clear img_filtered

% Salt&Pepper Noise
array = img01_salt_pepper;
img_filtered(:,:,1) = wiener2(array(:,:,1),[k k]);
img_filtered(:,:,2) = wiener2(array(:,:,2),[k k]);
img_filtered(:,:,3) = wiener2(array(:,:,3),[k k]);
img_filtered = remove_edge(img_filtered, 20);

figure
imshow(img_filtered)
title('Salt&Pepper Noise Filtered with Wiener filter')
clear img_filtered

% Shot (Poisson) Noise
array = img01_shot;
img_filtered(:,:,1) = wiener2(array(:,:,1),[k k]);
img_filtered(:,:,2) = wiener2(array(:,:,2),[k k]);
img_filtered(:,:,3) = wiener2(array(:,:,3),[k k]);
img_filtered = remove_edge(img_filtered, 20);
figure
imshow(img_filtered)
title('Shot (Poisson) Noise Filtered with Wiener filter')
clear img_filtered

% Speckle (Multiplicative) Noise
array = img01_speckle;
img_filtered(:,:,1) = wiener2(array(:,:,1),[k k]);
img_filtered(:,:,2) = wiener2(array(:,:,2),[k k]);
img_filtered(:,:,3) = wiener2(array(:,:,3),[k k]);
img_filtered = remove_edge(img_filtered, 20);
figure
imshow(img_filtered)
title('Speckle (Multiplicative) Noise Filtered with Wiener filter')
clear img_filtered

% Uniform Additve Noise
array = img01_unif_additive;
img_filtered(:,:,1) = wiener2(array(:,:,1),[k k]);
img_filtered(:,:,2) = wiener2(array(:,:,2),[k k]);
img_filtered(:,:,3) = wiener2(array(:,:,3),[k k]);
img_filtered = remove_edge(img_filtered, 20);
figure
imshow(img_filtered)
title('Uniform Additve Noise Filtered with Wiener filter')
clear img_filtered

%% Median Filter
n = 2;
% Additive White Gaussian Noise 
array = img01_AWGN;
img_filtered(:,:,1) = mymedfilt(array(:,:,1),n);
img_filtered(:,:,2) = mymedfilt(array(:,:,2),n);
img_filtered(:,:,3) = mymedfilt(array(:,:,3),n);
img_filtered = remove_edge(img_filtered, 20);
figure
imshow(img_filtered)
title('Additive White Gaussian Noise Filtered with Median filter')
clear img_filtered

% Salt&Pepper Noise
array = img01_salt_pepper;
img_filtered(:,:,1) = mymedfilt(array(:,:,1),n);
img_filtered(:,:,2) = mymedfilt(array(:,:,2),n);
img_filtered(:,:,3) = mymedfilt(array(:,:,3),n);
img_filtered = remove_edge(img_filtered, 20);

figure
imshow(img_filtered)
title('Salt&Pepper Noise Filtered with Median filter')
clear img_filtered

% Shot (Poisson) Noise
array = img01_shot;
img_filtered(:,:,1) = mymedfilt(array(:,:,1),n);
img_filtered(:,:,2) = mymedfilt(array(:,:,2),n);
img_filtered(:,:,3) = mymedfilt(array(:,:,3),n);
img_filtered = remove_edge(img_filtered, 20);
figure
imshow(img_filtered)
title('Shot (Poisson) Noise Filtered with Median filter')
clear img_filtered

% Speckle (Multiplicative) Noise
array = img01_speckle;
img_filtered(:,:,1) = mymedfilt(array(:,:,1),n);
img_filtered(:,:,2) = mymedfilt(array(:,:,2),n);
img_filtered(:,:,3) = mymedfilt(array(:,:,3),n);
img_filtered = remove_edge(img_filtered, 20);
figure
imshow(img_filtered)
title('Speckle (Multiplicative) Noise Filtered with Median filter')
clear img_filtered

% Uniform Additve Noise
array = img01_unif_additive;
img_filtered(:,:,1) = mymedfilt(array(:,:,1),n);
img_filtered(:,:,2) = mymedfilt(array(:,:,2),n);
img_filtered(:,:,3) = mymedfilt(array(:,:,3),n);
img_filtered = remove_edge(img_filtered, 20);
figure
imshow(img_filtered)
title('Uniform Additve Noise Filtered with Median filter')
clear img_filtered

%% Question 2 - Part 1 (Image Processing)
clear all
close all


img1 = imread('01.jpg');
img2 = imread('02.jpg');

img1_gray = rgb2gray(img1);
img2_gray = rgb2gray(img2);

Fourier_img1 = fftshift(fft(img1_gray));
Fourier_img2 = fftshift(fft(img2_gray));


figure
subplot(231)
imshow(img1_gray);
title('Photo');
subplot(232)
imshow(abs(Fourier_img1), [min(min(abs(Fourier_img1))), max(max(abs(Fourier_img1)))]);
title('Abs');
subplot(233)
imshow(angle(Fourier_img1));
title('Phase');

subplot(234)
imshow(img2_gray)
title('Photo');
subplot(235)
imshow(abs(Fourier_img2), [min(min(abs(Fourier_img2))), max(max(abs(Fourier_img2)))]);
title('Abs');
subplot(236)
imshow(angle(Fourier_img2));
title('Phase');

%%

fft1 = fft(img1_gray);
fft2 = fft(img2_gray);
fft3 = abs(fft1).*exp(1j.*angle(fft2));
fft4 = abs(fft2).*exp(1j.*angle(fft1));
img4 = real(ifft(fft4));
img3 = real(ifft(fft3));


figure

subplot(221)
imshow(img1_gray);
title('Img 1');

subplot(222)
imshow(img2_gray);
title('Img 2');

subplot(223)
imshow(img3);
title('Abs1 & Phase2');

subplot(224)
imshow(img4);
title('Abs2 & Phase1');


figure

subplot(221)
imshow(img1_gray);
title('Img 1');

subplot(222)
imshow(img2_gray);
title('Img 2');

subplot(223)
imshow(img3, [min(min(img3)), max(max(img3))]);
title('Abs1 & Phase2');

subplot(224)
imshow(img4, [min(min(img4)), max(max(img4))]);
title('Abs2 & Phase1');

%% Sound Processing
clc
clear
[sounds, fs_sounds] = audioread('sound.wav');
[singing, fs_sing] = audioread('singing.wav');

%% Zero Phase
fft_sing = fft(singing);
newfft = abs(fft_sing);
newsound = ifft(newfft);

sound(newsound, fs_sing);
audiowrite('NoPhaseSound.wav', newsound, fs_sing);

%% Sound with the Same Sampling Frequency

t_sing = 0:1/fs_sing:(length(singing)-1)/fs_sing;
t_sounds = 0:1/fs_sounds:(length(sounds)-1)/fs_sounds;
T_sing = max(t_sing);
T_sounds = max(t_sounds);

f = lcm(fs_sounds, fs_sing);
t = 0 :  1/f : min(T_sing, T_sounds);

sing_spline = spline(t_sing, singing, t);
sound_spline = spline(t_sounds, sounds, t);

sing_new = sing_spline(1:90:end);
sounds_new = sound_spline(1:90:end);
f = f/90;

%%

fft_sounds = fft(sounds_new);
fft_sing = fft(sing_new);

fft1 = abs(fft_sounds).*exp(1j*angle(fft_sing));
fft2 = abs(fft_sing).*exp(1j*angle(fft_sounds));

sound1 = real(ifft(fft1));
sound2 = real(ifft(fft2));

audiowrite('AbsSoundPhaseSing.wav', sound1, f);
audiowrite('AbsSingPhaseSound.wav', sound2, f);

%% Question 3
clc
clear all
close all

% Load Image
img = imread('coins3.jpg');
img_gray  = rgb2gray(img);


% Apply Gaussian Filter
G = fspecial('gaussian', [5,5], 2);
smooth = imfilter(img_gray,G,'same');


% Convert into a binary photo
bin = smooth < 45;

% Edge detection
Lx = [-1 0 1; -2 0 2; -1 0 1];
Ly = [1 2 1; 0 0 0; -1 -2 -1];

filtx = conv2(Lx, double(bin));
filty = conv2(Ly, double(bin));
final = sqrt(filtx.^2 + filty.^2);
edges = imcomplement(final);
binedges = (edges ~= 1);


% Detect circle using Hough transform
x0 = [];
y0 = [];
radius = [];
for r = 1:50
    [Potentialcircle(r).y0, Potentialcircle(r).x0, Potentialcircle(r).Accumulator(:,:)] = houghcircle(binedges, r, 150);
    x0 = [x0; Potentialcircle(r).x0];
    y0 = [y0; Potentialcircle(r).y0];
    radius = [radius, repmat(r, 1, length(Potentialcircle(r).y0))];

end

circles = [x0, y0, radius'];

% Remove Redundant Circles
[C ia ic] = unique(circles(:,1:2), 'rows');
circles = circles(ia, :);

[C ia ic] = unique(circles(:,2:3), 'rows');
circles = circles(ia, :);

[C ia ic] = unique(circles(:,[1,3]), 'rows');
circles = circles(ia, :);

s = size(circles);
disp('Number of Circles = ');
disp(s(1))

Circles = array2table(circles);
Circles.Properties.VariableNames{1} = 'x0';
Circles.Properties.VariableNames{2} = 'y0';
Circles.Properties.VariableNames{3} = 'R';
Circles

%% Final Plots
figure;
surf(Potentialcircle(48).Accumulator)
title('Accumulator Matrix for r = 48px')
xlabel('a')
ylabel('b')

figure
hold on
imshow(img)
viscircles(circles(:,1:2),circles(:,3) ,'LineStyle','--');


figure
subplot(231)
imshow(img)
title('Photo')

subplot(232)
imshow(img_gray)
title('Grayscale Photo')

subplot(233)
imshow(smooth)
title('Gaussian Filtered')

subplot(234)
imshow(bin)
title('Binary Photo')

subplot(2,3,[5,6])
imshow(binedges)
title('Edges')



