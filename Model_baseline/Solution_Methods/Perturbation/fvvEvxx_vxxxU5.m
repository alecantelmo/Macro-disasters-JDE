%
% © Copyright, Oren Levintal, June 13, 2016.

order=5;

% Vxx0(X)Vxxx0

derivs=[2,2,3];
perms1=[];
perms2=perms1;
Matrix=[];

result= permutekron3(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,Ui,fvv,Vxx0,Vxxx0);


% Vxx0(X)[Vxxx1*P(zeta^(X2)(X)hx)]

perms1=[2,1,3,4,5;2,3,1,4,5];
perms2=[1,2;1,2];
Matrix=kron(Ix2,kron(Ezeta2,hx));

result=result+permutekron3(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,Ui,fvv,Vxx0,Vxxx1);

% Vxx0(X)[Vxxx1*(zeta^(X3))]

perms1=[];
perms2=[];
Matrix=Ix2_Ezeta3;

result=result+permutekron3(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,Ui,fvv,Vxx0,Vxxx1);

% (Vxx1*P(hx(X)zeta))(X)[Vxxx1*P(zeta(X)hx^(X2))]

perms1=[1,3,2,4,5;3,1,2,4,5;1,2,3,5,4;1,3,2,5,4;3,1,2,5,4];
perms2=[1,2;1,2;1,2;1,2;1,2];
Matrix=kron(kron(hx,Ezeta2),hx2);

result=result+permutekron3(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,Ui,fvv,Vxx1,Vxxx1);

% (Vxx1*P(hx(X)zeta))(X)[Vxxx1*P(zeta^(X2)(X)hx)]

perms1=[2,1,3,4,5;2,3,1,4,5;1,2,3,5,4;2,1,3,5,4;2,3,1,5,4];
perms2=[1,2;1,2;1,2;1,2;1,2];
Matrix=kron(hx_Ezeta3,hx);

result=result+permutekron3(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,Ui,fvv,Vxx1,Vxxx1);

% (Vxx1*P(hx(X)zeta))(X)[Vxxx1*P(zeta^(X3))]

perms1=[1,2,3,5,4];
perms2=[1,2];
Matrix=kron(hx,Ezeta4);

result=result+permutekron3(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,Ui,fvv,Vxx1,Vxxx1);

% (Vxx1*P(hx(X)zeta))(X)[Vxx1*(zeta(X)hxx)Omega1]

Matrix=kron(kron(hx,Ezeta2),Ix);
tempindex=permute(reshape([1:n_x^4],n_x,n_x,n_x,n_x),[1,2,4,3]);
tempindex=tempindex(:);
Matrix=Matrix+Matrix(tempindex,tempindex);
clear tempindex

fvv_Vxx1_Vxx1=[reshape([reshape([reshape(fvv,n_f*n_v,n_v)*Vxx1]',n_x^2*n_f,n_v)*Vxx1]',n_x^4,n_f)]';
tempresult=fvv_Vxx1_Vxx1*Matrix;

tempresult=[reshape([reshape(tempresult,n_f*n_x,n_x^3)]',n_x^3*n_f,n_x)*hxx]';
tempresult=[reshape(tempresult,n_x^5,n_f)]';

result=result+tempresult*(kron(Ix2,OMEGA_x.OMEGA1)*U5);

clear tempresult

% (Vxx1*zeta^(X2))(X)Vxxx0

perms1=[];
perms2=[];
Matrix=[];

result=result+permutekron3(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,Ui,fvv,Vxx1*Ezeta2,Vxxx0);

% (Vxx1*zeta^(X2))(X)[Vxxx1*P(zeta(X)hx^(X2))]

perms1=[1,3,2,4,5;3,1,2,4,5];
perms2=[1,2;1,2];
Matrix=kron(Ezeta3,hx2);

result=result+permutekron3(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,Ui,fvv,Vxx1,Vxxx1);

% (Vxx1*zeta^(X2))(X)[Vxxx1*P(zeta^(X2)(X)hx)]

perms1=[2,1,3,4,5;2,3,1,4,5];
perms2=[1,2;1,2];
Matrix=kron(Ezeta4,hx);

result=result+permutekron3(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,Ui,fvv,Vxx1,Vxxx1);

% (Vxx1*zeta^(X2))(X)[Vxxx1*zeta^(X3)]

perms1=[];
perms2=[];
Matrix=Ezeta5;

result=result+permutekron3(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,Ui,fvv,Vxx1,Vxxx1);

% (Vxx1*zeta^(X2))(X)[Vxx1*(zeta(X)hxx)Omega1]

Matrix=kron(Ezeta3,hxx);

tempresult=fvv_Vxx1_Vxx1*Matrix;

result=result+tempresult*(kron(Ix2,OMEGA_x.OMEGA1)*U5);

clear tempresult fvv_Vxx1_Vxx1