function [ fxxc,ind ] = chain2c_tensor(fv,fvv,vx,vxxc,ind,n_ind,maxload,sum)
%chain2c is compressed.
%
% © Copyright, Oren Levintal, June 13, 2016.

if isempty(ind)
    ind=cell(2,1);
end

%fvv*kron(vx,vx)*U2
[term1,ind{1}]=AkronBU2(fvv,vx,ind{1},n_ind,maxload,sum);

%fv*vxxc
[term2,ind{2}]=contraction1(fv,vxxc,ind{2},n_ind,maxload,sum);

fxxc=tplus(term1,term2,maxload);

if isfield(fv,'ptr2d')
    fxxc.ptr2d=fv.ptr2d;
end

end

