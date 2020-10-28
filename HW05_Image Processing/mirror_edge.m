function Y = mirror_edge (X,n)
    % This function extends the input array (image) by mirroring the edges
    % in order to enhance the quality of filtering process
    dim = size(X);

    Y(n+1:dim(1)+n, n+1:dim(2)+n, :) = X(:,:,:);
    
    Y(1:n,n+1:dim(2)+n,:) = flip(X(1:n,:,:));
    Y(:,1:n,:) = fliplr(Y(:,n+1:2*n,:));
    
    Y(dim(1)+n+1:dim(1)+2*n,:,:) = flip(Y(end-n+1:end,:,:));
    Y(:,dim(2)+n+1:dim(2)+2*n,:) = fliplr(Y(:,end-n+1:end,:));
    
end