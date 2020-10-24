function fc2D=fcReshape(fc1D)

m=size(fc1D,1); % number of unique pairs (handshakes)
n=(1+sqrt(1+4*2*m))/2;

fc2D=ones(n,n); % diagonals will stay at 1 because this is a correlation matrix
inds2fill=find(triu(fc2D,1));
fc2D(inds2fill)=fc1D;

% force matrix to by symmetric
for col = 1:n
    for row=col+1:n
        fc2D(row,col)=fc2D(col,row);
    end
end


