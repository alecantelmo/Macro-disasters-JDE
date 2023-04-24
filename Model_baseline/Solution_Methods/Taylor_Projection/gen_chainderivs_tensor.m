function [index,n]=gen_chainderivs_tensor(tilf,v,u,pi_,order,symparams,fname)
%
% © Copyright, Oren Levintal, June 13, 2016.

if order>4
    error('Differentiation order must not exceed 4')
end

n=find_n(pi_,u); % number of substitutions needed to eliminate u

z=[v;u];

% Differentiate with chain rules

[index]=getderivs_tensor(tilf,z,order,symparams,[fname '_tilf']);





