function [ fxxxxc,ind ] = chain4c_tensor(fv,fvv,fvvv,fvvvv,vx,vxxc,vxxxc,vxxxxc,M2,M3,M4,ind,n_ind,maxload,sum,varargin)
%chain4c is compressed.
%
% © Copyright, Oren Levintal, June 13, 2016.

if ~isempty(varargin)
    convertind=varargin{1};
    maxind=varargin{2};
    doi=1;
else
    doi=0;
end

if isempty(ind)
    ind=cell(8,1);
end

%fvvvv*kron(vx,vx,vx,vx)*U4
if doi==0
    [term1,ind{1}]=AkronBU4(fvvvv,vx,ind{1},n_ind,maxload,sum);
else
    [term1,ind{1}]=AkronBU4i(fvvvv,vx,ind{1},n_ind,maxload,sum,convertind,maxind);
end


%fvvv*kron(kron(vx,vx)*U,vxxc)*kron(W2,W2)*OMEGA2*U4
[term2,ind{2}]=AkronB1U2B2(fvvv,vx,vxxc,ind{2},n_ind,maxload,sum);
[term2,ind{3}]=contraction1(term2,M2,ind{3},n_ind,maxload,sum);


sumterm=tplus(term1,term2); clear term1 term2

%fvv*kron(vx,vxxxc)*kron(I,W3)*OMEGA3*U4
[term3,ind{4}]=contraction2(fvv,vx,vxxxc,ind{4},n_ind,maxload,sum);
[term3,ind{5}]=contraction1(unfold(term3),M3,ind{5},n_ind,maxload,sum);
sumterm=tplus(sumterm,term3); clear term3

%fvv*kron(vxxc,vxxc)*U*W*kron(W2,W2)*OMEGA4*U4
[term4,ind{6}]=AkronBU2(fvv,vxxc,ind{6},n_ind,maxload,sum);

[term4,ind{7}]=contraction1(unfold(term4),M4,ind{7},n_ind,maxload,sum);
sumterm=tplus(sumterm,term4); clear term4

[term5,ind{8}]=contraction1(fv,vxxxxc,ind{8},n_ind,maxload,sum);

fxxxxc=tplus(sumterm,term5);


end

