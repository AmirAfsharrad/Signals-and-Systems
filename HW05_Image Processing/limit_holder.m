function y = limit_holder (x,n)
    % this function operates as follows:
    % y = x     if 1 < x < n
    % y = 1     if x <= 1
    % y = n     if x >= n
    y = (x-1).*heaviside(x-1)+1;
    y(y>n) = n;
end
