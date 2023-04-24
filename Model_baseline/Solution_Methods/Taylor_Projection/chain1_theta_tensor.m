function [ R,ind,term1,term2 ] = chain1_theta_tensor( fv,fvv,vx,vtheta,vxtheta,ind,n_ind,maxload,sum )
%
% © Copyright, Oren Levintal, June 13, 2016.

if isempty(ind)
    ind=cell(2,1);
end

n_f=fv.tsize(1);

n_x=vx.tsize(2);
n_theta=vtheta.tsize(2);
n_s=max([size(fv.vals,1),size(vx.vals,1),size(vtheta.vals,1)]);

if ~isempty(fvv.vals)
    [term1,ind{1}]=contraction2(fvv,vx,vtheta,ind{1},n_ind,maxload,sum);
    term1=permutecols(term1,[2,1]);
else
    term1=sptensor(n_f,[n_x,n_theta],n_s);
end

[term2,ind{2}]=contraction1(fv,vxtheta,ind{2},n_ind,maxload,sum);

R=tplus(term1,fold(term2,n_x,n_theta),maxload);


end

