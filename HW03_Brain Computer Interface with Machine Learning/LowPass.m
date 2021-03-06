function y = LowPass(x)
%DOFILTER Filters input x and returns output y.

% MATLAB Code
% Generated by MATLAB(R) 9.0 and the Signal Processing Toolbox 7.2.
% Generated on: 24-Mar-2018 23:09:17

persistent Hd;

if isempty(Hd)
    
    Fpass = 60;    % Passband Frequency
    Fstop = 65;    % Stopband Frequency
    Apass = 1;     % Passband Ripple (dB)
    Astop = 60;    % Stopband Attenuation (dB)
    Fs    = 2400;  % Sampling Frequency
    
    h = fdesign.lowpass('fp,fst,ap,ast', Fpass, Fstop, Apass, Astop, Fs);
    
    Hd = design(h, 'equiripple', ...
        'MinOrder', 'any', ...
        'StopbandShape', 'flat');
    
    
    
    set(Hd,'PersistentMemory',true);
    
end

y = filter(Hd,x);
y = filter(Hd,x);

