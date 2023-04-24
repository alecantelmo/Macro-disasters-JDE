function [R]=tscalar(A)
% returns a scalar tensor.
%
% © Copyright, Oren Levintal, June 13, 2016.

A=tvec(A);

[i,j,vals]=tfind(A);

Rvals=zeros([A.tsize(1),size(A.vals,1)]);
Rvals(i,:)=vals';

R=sptensor(1);

R.vals=Rvals(:);

end