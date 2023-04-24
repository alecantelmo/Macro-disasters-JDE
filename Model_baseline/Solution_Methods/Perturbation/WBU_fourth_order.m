%
% © Copyright, Oren Levintal, June 13, 2016.

temp1=reshape(W4, unique*n_x^2,n_x^2)*kron_hx_hx;
temp1=reshape(temp1,unique,n_x^4);
temp2=kron_hx_hx*reshape(U4,n_x^2,n_x^2*unique);
temp2=reshape(temp2,n_x^4,unique);
W4BU4=temp1*temp2;


Matrix=hx2_Ezeta2;
perms1=[1,3,2,4;1,3,4,2;3,1,2,4;3,1,4,2;3,4,1,2];
tempresult = permutekronB2(4,perms1,n_x,Matrix,W4);

Matrix=hx_Ezeta3;
perms1=[1,2,4,3;1,4,2,3;4,1,2,3];
tempresult = tempresult+permutekronB2(4,perms1,n_x,Matrix,W4 );

tempresult=tempresult+W4*Ezeta4;

W4BU4=W4BU4+tempresult*U4;
