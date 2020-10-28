function y = doFilt(x, h)
    %{
        Inputs:
            h: Impulse response of the filter
            x: Matrix of inputs (time,signals)
        Outut:
            y: Matrix of outputs (time,signals)
    %}
    y = ifft(fft(x).*fft(h)');
end