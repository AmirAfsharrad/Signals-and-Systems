function y = BPF(x, N, fl, fh, Fs)


% Truncated Sinc Filter
wh = 2*((fh)/Fs)*pi;
wl = 2*((fl)/Fs)*pi;

Ll = (N-1)/2 * (wl/pi);
numl = sinc(-Ll:wl/pi:Ll).*(hamming(length(-Ll:wl/pi:Ll))');



Lh = (N-1)/2 * (wh/pi);
numh = sinc(-Lh:wh/pi:Lh).*(hamming(length(-Lh:wh/pi:Lh))');

num = numh/polyval(numh, 1) - numl/polyval(numl, 1);
% c= 0.5*(fl+fh);
% OMEGA = 2*c*pi/Fs;

% yh = filter(numh,polyval(numh, exp(1j*OMEGA)), x);
% yl = filter(numl,polyval(numl, exp(1j*OMEGA)), x);
y = filter(num, 1, x);
end
% y = yh - yl;


%
% Truncated Sinc Filter
% wl = 2*(fl/Fs)*pi;
% wh = 2*(fh/Fs)*pi;
% 
% Ll = (N-1)/2 * (wl/pi);
% Lh = (N-1)/2 * (wh/pi);
% 
% numl = sinc(-Ll:wl/pi:Ll).*(hamming(length(-Ll:wl/pi:Ll))');
% numh = sinc(-Lh:wh/pi:Lh).*(hamming(length(-Lh:wh/pi:Lh))');
% 
% c = (wl+wh)/2;
% num = numh - numl;
% OMEGA = 2*c/Fs;
% freqz(num, polyval(num, exp(1j*OMEGA)))
% 
% %freqz(num,polyval(num,1))
% y = filter(num,polyval(num, exp(1j*OMEGA)), x);
% 
% %y = filter(Hd,x);