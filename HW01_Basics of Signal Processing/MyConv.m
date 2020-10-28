function Y = MyConv(u,v)
    c = [u zeros(1,length(v)-1)];
    r = [u(1) zeros(1,length(v)-1)];
    ConvolutionMatrix = toeplitz(c,r);
    Y = ConvolutionMatrix*v';
    Y = Y';
end

