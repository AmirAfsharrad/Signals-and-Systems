function [ Y ] = SoundMaker( gain, note, interval, Fs)
    a = length(gain);
    Y = 0;
    for i = 1 : a
        f = 440.*2.^(note(i)./12);
        omega = 2 .* pi .* f;
        t = 0 : 1/Fs : interval(i);
        Y = [Y , gain(i) .* sin(omega.*t)];
    end
end