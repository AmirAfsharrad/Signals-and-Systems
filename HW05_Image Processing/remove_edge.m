function Y = remove_edge(X,n)
    % This function removes edges from around the input array (image) which
    % has been added before performing any filter. Obviously, this function
    % is to be used after filtering
    Y = X(n+1:end-n,n+1:end-n,:);
end