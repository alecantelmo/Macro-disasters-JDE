function [ fxxxc,ind ] = chain3c_tensor(fv,fvv,fvvv,vx,vxxc,vxxxc,M2,ind,n_ind,maxload,sum)
%chain3c is compressed.
%
% © Copyright, Oren Levintal, June 13, 2016.

if isempty(ind)
    ind=cell(4,1);
end

%fvvv*kron(vx,vx,vx)*U3

[term1,ind{1}]=AkronBU3(fvvv,vx,ind{1},n_ind,maxload,sum);

[term2,ind{2}]=contraction2(fvv,vx,vxxc,ind{2},n_ind,maxload,sum);

[term2,ind{3}]=contraction1(unfold(term2),M2,ind{3},n_ind,maxload,sum);

sumterm=tplus(term1,term2,maxload); clear term1 term2

[term3,ind{4}]=contraction1(fv,vxxxc,ind{4},n_ind,maxload,sum);

fxxxc=tplus(sumterm,term3,maxload);

if isfield(fv,'ptr2d')
    fxxxc.ptr2d=fv.ptr2d;
end


end

