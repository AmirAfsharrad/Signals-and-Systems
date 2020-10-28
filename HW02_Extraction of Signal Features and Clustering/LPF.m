function y = LPF(x, N, fc, Fs)

% Truncated Sinc Filter
w0 = 2*(fc/Fs)*pi;

L = (N-1)/2 * (w0/pi);
num = sinc(-L:w0/pi:L);
%freqz(num,polyval(num,1))
y = filter(num,polyval(num, 1), x);
%y = filter(Hd,x);