function output = my_reshape(input)
%     input(64,n,20);
    
    if(length(size(input))==3)
        m = size(input,3);
        n = size(input,2);
        output = zeros(m,64*n);
        for i = 1 : m
            output(i,:) = reshape(input(:,:,i),1,64*n);
        end
    end
    if(length(size(input))==4)
        m = size(input,4);
        n = size(input);
        n = n(2:3);
        output = zeros(m,64*prod(n));
        for i = 1 : m
            output(i,:) = reshape(input(:,:,:,i),1,64*prod(n));
        end
    end
end