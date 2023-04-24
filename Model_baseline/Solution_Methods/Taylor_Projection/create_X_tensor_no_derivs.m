function [X,ind]=create_X_tensor_no_derivs(N,x_c0,M2,M3,ind,n_ind,maxload,sum)
%
% © Copyright, Oren Levintal, June 13, 2016.

if isempty(ind)
    ind=cell(5,1);
end
n_x=x_c0.tsize(1);

n_s=size(x_c0.vals,1);
temp=sptensor(1);
temp.vals=repmat(temp.vals,n_s,1);
X=vconcat(temp,x_c0);

if N>=2
    [temp,ind{1}]=contraction2(M2,x_c0,x_c0,ind{1},n_ind,maxload,sum);
    X=vconcat(X,unfold(temp));
end

if N>=3
    [temp,ind{3}]=contraction3(M3,x_c0,x_c0,x_c0,ind{3},n_ind,maxload,sum);
    X=vconcat(X,unfold(temp));
end
