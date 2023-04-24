function ten=permutecols(ten,ind)
% permute columns of a sparse sptensor
%
% © Copyright, Oren Levintal, June 13, 2016.

ind=ind(:)';
ten.cols=ten.cols(:,ind);
ten.tsize=ten.tsize([1,1+ind]);
end

