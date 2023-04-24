%
% © Copyright, Oren Levintal, June 13, 2016.



temp1=(reshape(W5, unique*n_x^3,n_x^2)*kron_hx_hx);
temp1=reshape(temp1,unique,n_x^5);
temp2=kron(hx,kron_hx_hx)*reshape(U5,n_x^3,n_x^2*unique);
temp2=reshape(temp2,n_x^5,unique);
W5BU5=temp1*temp2;

order=5;
Matrix=hx3_Ezeta2;
perms1=[1,3,2,4,5;1,3,4,2,5;1,3,4,5,2;3,1,2,4,5;3,1,4,2,5;3,1,4,5,2;3,4,1,2,5;3,4,1,5,2;3,4,5,1,2];
tempresult = permutekronB2(order,perms1,n_x,Matrix,W5 );

Matrix=hx2_Ezeta3;
perms1=[1,2,4,3,5;1,2,4,5,3;1,4,2,3,5;1,4,2,5,3;1,4,5,2,3;4,1,2,3,5;4,1,2,5,3;4,1,5,2,3;4,5,1,2,3];
tempresult = tempresult+permutekronB2(order,perms1,n_x,Matrix,W5 );

Matrix=hx_Ezeta4;
perms1=[1,2,3,5,4;1,2,5,3,4;1,5,2,3,4;5,1,2,3,4];
tempresult = tempresult+permutekronB2(order,perms1,n_x,Matrix,W5);

tempresult = tempresult+W5*Ezeta5;

W5BU5=W5BU5+tempresult*U5;

clear tempresult