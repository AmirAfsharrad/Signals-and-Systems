function [] = plotFFT( Y1, Fs, inf, sup, titlestr, ylabelstr )
L = length(Y1);
P2_1 = abs(Y1/L);
P1_1 = P2_1(1:L/2+1);
P1_1(2:end-1) = 2*P1_1(2:end-1);

f = Fs*(0:(L/2))/L;

plot(f, P1_1)
title(titlestr)
xlabel('f (Hz)')
ylabel(ylabelstr)
xlim([inf,sup]);
end

