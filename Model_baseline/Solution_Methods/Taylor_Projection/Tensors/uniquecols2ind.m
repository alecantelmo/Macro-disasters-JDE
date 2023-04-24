function [ loc ] = uniquecols2ind( cols,n )
%convert unique columns to index (like u2c). Columns are assumed to satisfy
%a1>=a2>=a3
%
% © Copyright, Oren Levintal, June 13, 2016.

[n_unique,k]=size(cols);

if k>1
    sortM=double(cols(:,end:-1:1)); % the formula below assumes indices a1<=a2<=a3, so i have to change location of indices because i use the convention a3>=a2>=a1

    tempmat=repmat(n+k:-1:n+1,n_unique,1);
    loc=(repmat(prod(n+k-1:-1:n),n_unique,1)-prod(tempmat-repmat(sortM(:,1),1,k),2))/factorial(k)+1;

    for j=2:k-1
        tempmat=repmat(n+k-j+1:-1:n+1,n_unique,1);
        loc=loc+(prod(tempmat-repmat(sortM(:,j-1),1,k-j+1),2)-prod(tempmat-repmat(sortM(:,j),1,k-j+1),2))/factorial(k-j+1);
    end
    loc=loc+sortM(:,k)-sortM(:,k-1);
else
    loc=cols;
end

end

