function [ R,ind ] = chain1_tensor(fv,vx,ind,n_ind,maxload,sum)
%1st order chain rule.
%
% © Copyright, Oren Levintal, June 13, 2016.

[R,ind]=contraction1(fv,vx,ind,n_ind,maxload,sum);

end

