function y = multiConv(b)
    y = conv(b(1,:),b(2,:));
    for i = 3 : size(b,1)
       y =  conv(b(i,:),y);
    end
end