function [b, a] = LowPassFilterDesigner (r, w0, percent, AM_coef)
n = 2;
percent = percent/100;
while(1)
    b = zeros(1,n);
    w1 = w0 + (AM_coef)*(1-w0);
    if(mod(n,2) == 0)
       b(1:n/2) = linspace(w1,1-2*(1-w1)/n,n/2); 
       b(n/2 + 1 : end) = -b(1:n/2);
       b = -r*exp(1j*b*pi);
    else
       m = n-1;
       b(1:m/2) = linspace(w1,1-2*(1-w1)/m,m/2); 
       b(m/2 + 1 : end-1) = -b(1:m/2);
       b(n) = 1;
       b = -r*exp(1j*b*pi); 
    end
    b = [ones(n,1) , b'];
    b = multiConv(b);
    a = polyval(b,1);
    [h, w] = freqz(b,a,2001);
    total_energy = sum(abs(h).^2);
    E = 0;
    for i = 1 : length(w)
        E = E + abs(h(i).^2);
        if(E/total_energy > percent)
           break
        end
    end
    if(w(i)<=w0*pi)
        disp(['Filter Order = ',num2str(n)]);
        disp([num2str(percent*100),' percent energy normalized frequency = ',num2str(w(i)/pi)]);
        freqz(b,a,2001)
        break
    end
    
    n = n+1;
end