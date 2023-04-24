function [A]=changecols(A,newcols,coldim,coli)
%Change the cols of tensor A
%
% © Copyright, Oren Levintal, June 13, 2016.

newcols=intarray(newcols);
A.cols(:,coli)=newcols(A.cols(:,coli));
A.tsize(coli+1)=intarray(coldim);

end