function Y = mymedfilt(X,n)
    % This function performs the median filter on the input array (image)
    Y = zeros(size(X));
    for i = 1 : size(X,1)
        for j = 1 : size(X,2)
            Z = X(limit_holder(i-n,size(X,1)):limit_holder(i+n,size(X,1)),...
                limit_holder(j-n,size(X,2)):limit_holder(j+n,size(X,2)));
            Y(i,j) = median(median(Z));
        end
    end
end