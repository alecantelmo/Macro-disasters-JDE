function [ R,ind] = chain2c_theta_tensor(fv,fvv,fvvv,vx,vxxc,vtheta,vxtheta,vxxctheta,M2,ind,n_ind,maxload,sum )
%
% © Copyright, Oren Levintal, June 13, 2016.

if isempty(ind)
    ind=cell(5,1);
end

n_f=fv.tsize(1);

n_x=vx.tsize(2);
n_theta=vtheta.tsize(2);
unique2=nchoosek(n_x+1,2);
n_s=max([size(fv.vals,1),size(vx.vals,1),size(vtheta.vals,1)]);


if ~isempty(fvvv.vals)

    [term1,ind{1}]=AkronB1U2B2(fvvv,vx,vtheta,ind{1},n_ind,maxload,sum);
    term1=fold(term1,n_theta,unique2);

    term1=col2ptr(term1,2);
    
else
    term1=sptensor(n_f,[unique2,n_theta],n_s);
    term1=col2ptr(term1,1);
end

[term2,ind{2}]=contraction2(fvv,vx,vxtheta,ind{2},n_ind,maxload,sum);
term2=foldj(term2,1,[n_x,n_theta]);

term2=ptr1d(unfold(col2ptr(term2,2)));

[term2,ind{3}]=contraction1(term2,M2,ind{3},n_ind,maxload,sum);
term2=ptr2d(term2,n_f,n_theta);
term2=ptr2col(term2,2); % tsize=[n_f*n_theta,unique2]
term2=col2ptr(term2,1);


sumterm=tplus(term1,term2,maxload); 

[term3,ind{4}]=contraction2(fvv,vxxc,vtheta,ind{4},n_ind,maxload,sum);

term3=col2ptr(term3,2);
sumterm=tplus(sumterm,term3,maxload); 


[term4,ind{5}]=contraction1(fv,vxxctheta,ind{5},n_ind,maxload,sum);
term4=fold(term4,unique2,n_theta);
term4=col2ptr(term4,1);



R=tplus(sumterm,term4,maxload);

R=ptr2col(R,1);

end

