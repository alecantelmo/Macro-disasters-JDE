%
% © Copyright, Oren Levintal, June 13, 2016.

% W3(hx^3)U3
W3BU3=reshape([reshape(reshape([reshape(W3, unique*n_x^2,n_x)*hx]',n_x*unique*n_x,n_x)*hx,n_x,unique*n_x^2)]',unique,n_x^3)...
    *reshape(hx*reshape(U3, n_x,n_x^2*unique),n_x^3,unique);

% W3*E(P(kron(hx,zeta,zeta)))+W3*E(kron(zeta,zeta,zeta))
Matrix=kron(hx,Ezeta2);
perms1=[1,3,2;3,1,2];
tempresult = permutekronB2(3,perms1,n_x,Matrix,W3);

tempresult=tempresult+W3*Ezeta3;

% final result
W3BU3=W3BU3+tempresult*U3;

