function [ R ,ind] = chain3c_theta_tensor(fv,fvv,fvvv,fvvvv,...
    vx,vxxc,vxxxc,vtheta,vxtheta,vxxctheta,vxxxctheta,...
    M2,M3,M4,M5,ind,n_ind,maxload,sum )
%
% © Copyright, Oren Levintal, June 13, 2016.

if isempty(ind)
    ind=cell(15,1);
end

n_f=fv.tsize(1);
n_v=fv.tsize(2);
n_x=vx.tsize(2);
n_theta=vtheta.tsize(2);
unique2=nchoosek(n_x+1,2);
unique3=nchoosek(n_x+2,3);

n_s=max([size(fv.vals,1),size(vx.vals,1),size(vtheta.vals,1)]);

if ~isempty(fvvvv.vals)
    [term1,ind{1}]=AkronB1U3B2(fvvvv,vx,vtheta,ind{1},n_ind,maxload,sum);
    term1=fold(term1,n_theta,unique3);
    term1=col2ptr(term1,2);
else
    term1=sptensor(n_f,[unique3,n_theta],n_s);
    term1=col2ptr(term1,1);
end

[term2,ind{2}]=AkronB1U2B2(fvvv,vx,vxtheta,ind{2},n_ind,maxload,sum);

term2=fold(term2,n_x,n_theta,unique2);

term2=ptr1d(unfold(col2ptr(term2,2)));

[term2,ind{3}]=contraction1(term2,M2,ind{3},n_ind,maxload,sum);
term2=ptr2d(term2,n_f,n_theta);% tsize=[n_f*n_theta,unique3]
term2=ptr2col(term2,2);% tsize=[n_f,unique3,n_theta]
term2=col2ptr(term2,1);% tsize=[n_f*unique3,n_theta]

sumterm=tplus(term1,term2,maxload);


fvvv=col2ptr(fvvv,1); %n_f*n_v,n_v,n_v
fvvv=ptr1d(fvvv); %n_f*n_v,n_v,n_v
[term3,ind{4}]=contraction2(fvvv,vx,vxxc,ind{4},n_ind,maxload,'vec'); %n_f*n_v,unique2,n_x
term3=ptr2d(term3,n_f,n_v); %n_f*n_v,unique2,n_x
term3=ptr2col(term3,1); %n_f,n_v,unique2,n_x

term3=fold(term3,n_v,unique2*n_x);
term3=col2ptr(term3,2);
term3=ptr1d(term3);
[term3,ind{12}]=contraction1(term3,vtheta,ind{12},n_ind,maxload,sum);
term3=fold(ptr2col(ptr2d(term3,n_f,unique2*n_x),2),n_theta,unique2,n_x);

term3=col2ptr(term3,1); %n_f*n_theta,unique2,n_x
term3=ptr1d(unfold(term3)); %n_f*n_theta,unique2*n_x
[term3,ind{5}]=contraction1(term3,M3,ind{5},n_ind,maxload,sum);
term3=ptr2d(term3,n_f,n_theta); %n_f*n_theta,unique3
term3=ptr2col(term3,2);%n_f,unique3,n_theta
term3=col2ptr(term3,1);%n_f*unique3,n_theta

sumterm=tplus(sumterm,term3,maxload);

% term4 

rfvv=col2ptr(fvv,1); %n_f*n_v,n_v
rfvv=ptr1d(rfvv); %n_f*n_v,n_v
[term4,ind{6}]=contraction1(rfvv,vx,ind{6},n_ind,maxload,'vec'); %n_f*n_v,n_x
term4=ptr2d(term4,n_f,n_v); %n_f*n_v,n_x
term4=ptr2col(term4,1); %n_f,n_v,n_x

term4=col2ptr(term4,2); %n_f*n_x,n_v
term4=ptr1d(term4);
[term4,ind{13}]=contraction1(term4,col2ptr(fold(vxxctheta,unique2,n_theta),1),ind{13},n_ind,maxload,sum);

% this is special because term4 ptr has 3D
term4.ptr2d=intarray([n_f*n_x,unique2]);
term4=ptr2col(term4,1); %n_f*n_x,unique2,n_theta
term4.ptr2d=intarray([n_f,n_x]);
term4=ptr2d(term4,n_f,n_x);
term4=ptr2col(term4,2); %n_f,unique2,n_x,n_theta
term4=col2ptr(term4,3); %n_f*n_theta,unique2,n_x

term4=ptr1d(unfold(term4)); %n_f*n_theta,unique2*n_x
[term4,ind{7}]=contraction1(term4,M4,ind{7},n_ind,maxload,sum);
term4=ptr2d(term4,n_f,n_theta); %n_f*n_theta,unique3
term4=ptr2col(term4,2);%n_f,unique3,n_theta
term4=col2ptr(term4,1);%n_f*unique3,n_theta

sumterm=tplus(sumterm,term4,maxload);


% term5

[term5,ind{8}]=contraction1(rfvv,vxxc,ind{8},n_ind,maxload,'vec'); %n_f*n_v,unique2
term5=ptr2d(term5,n_f,n_v); %n_f*n_v,unique2
term5=ptr2col(term5,1); %n_f,n_v,unique2
term5=col2ptr(term5,2); %n_f*unique2,n_v
term5=ptr1d(term5);
[term5,ind{14}]=contraction1(term5,col2ptr(fold(vxtheta,n_x,n_theta),1),ind{14},n_ind,maxload,sum);

% this is special because term5 ptr has 3D
term5.ptr2d=intarray([n_f*unique2,n_x]);
term5=ptr2col(term5,1); %n_f*unique2,n_x,n_theta
term5.ptr2d=intarray([n_f,unique2]);
term5=ptr2d(term5,n_f,unique2);
term5=ptr2col(term5,2); %n_f,n_x,unique2,n_theta
term5=col2ptr(term5,3); %n_f*n_theta,n_x,unique2

term5=ptr1d(unfold(term5)); %n_f*n_theta,n_x*unique2
[term5,ind{9}]=contraction1(term5,M5,ind{9},n_ind,maxload,sum);
term5=ptr2d(term5,n_f,n_theta); %n_f*n_theta,unique3
term5=ptr2col(term5,2);%n_f,unique3,n_theta
term5=col2ptr(term5,1);%n_f*unique3,n_theta

sumterm=tplus(sumterm,term5,maxload);

% term6

[term6,ind{10}]=contraction1(rfvv,vxxxc,ind{10},n_ind,maxload,'vec'); %n_f*n_v,unique3
term6=ptr2d(term6,n_f,n_v); %n_f*n_v,unique3
term6=ptr2col(term6,1); %n_f,n_v,unique3
term6=col2ptr(term6,2); %n_f*unique3,n_v
term6=ptr1d(term6);
[term6,ind{15}]=contraction1(term6,vtheta,ind{15},n_ind,maxload,sum);
term6=ptr2d(term6,n_f,unique3); %n_f*unique3,n_theta

sumterm=tplus(sumterm,term6,maxload);

[term7,ind{11}]=contraction1(fv,col2ptr(fold(vxxxctheta,unique3,n_theta),1),ind{11},n_ind,maxload,sum);


R=tplus(sumterm,term7,maxload);

R=ptr2col(R,1);

end

