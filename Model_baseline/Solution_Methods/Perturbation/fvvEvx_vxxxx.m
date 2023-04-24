%
% © Copyright, Oren Levintal, June 13, 2016.

order=5;

% Vx0(X)Vxxxx0
derivs=[2,1,4];
perms1=[];
perms2=perms1;
Matrix=[];

result= permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx0,Vxxxx0);

% Vx0(X)(Vxxxx1*P(hx^(X2)(X)zeta^(X2))

perms1=[1,3,2,4,5;1,3,4,2,5;3,1,2,4,5;3,1,4,2,5;3,4,1,2,5];
perms2=[1,2;1,2;1,2;1,2;1,2];
Matrix=kron(Ix,kron(hx2,Ezeta2));

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx0,Vxxxx1);

% Vx0(X)(Vxxxx1*P(hx(X)zeta^(X3))

perms1=[1,2,4,3,5;1,4,2,3,5;4,1,2,3,5];
perms2=[1,2;1,2;1,2];
Matrix=kron(kron(Ix,hx),Ezeta3);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx0,Vxxxx1);

% Vx0(X)(Vxxxx1*zeta^(X4))

perms1=[];
perms2=[];
Matrix=Ix_Ezeta4;

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx0,Vxxxx1);

% Vx0(X)(Vxxx1*(zeta^(X2)(X)hxx)*Omega2

perms1=[];
perms2=[];
Matrix=[];

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx0,Vxxx1*kron(Ezeta2,hxx)*OMEGA_x.OMEGA2);

% (Vx1*zeta)(X)(Vxxxx1*P(zeta(X)hx^(X3))

perms1=[1,2,4,3,5;1,4,2,3,5;4,1,2,3,5];
perms2=[1,2;1,2;1,2];
Matrix=kron(Ezeta2,hx3);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx1,Vxxxx1);

% (Vx1*zeta)(X)(Vxxxx1*P(zeta^(X2)(X)hx^(X2))

perms1=[1,3,2,4,5;1,3,4,2,5;3,1,2,4,5;3,1,4,2,5;3,4,1,2,5];
perms2=[1,2;1,2;1,2;1,2;1,2];
Matrix=kron(Ezeta3,hx2);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx1,Vxxxx1);

% (Vx1*zeta)(X)(Vxxxx1*P(zeta^(X3)(X)hx)

perms1=[2,1,3,4,5;2,3,1,4,5;2,3,4,1,5];
perms2=[1,2;1,2;1,2];
Matrix=kron(Ezeta4,hx);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx1,Vxxxx1);

% (Vx1*zeta)(X)(Vxxxx1*zeta^(X4))

perms1=[];
perms2=[];
Matrix=Ezeta5;

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx1,Vxxxx1);

% (Vx1*zeta)(X)(Vxxx1*(P(zeta(X)hx)+zeta(X)zeta)*Omega2

Matrix=kron(kron(Ezeta2,hx),Ix);
tempindex=permute(reshape([1:n_x^4],n_x,n_x,n_x,n_x),[1,3,2,4]);
tempindex=tempindex(:);
Matrix=Matrix+Matrix(tempindex,tempindex)+kron(Ezeta3,Ix);
clear tempindex

tempresult=[reshape([reshape([reshape(fvv,n_f*n_v,n_v)*Vx1]',n_x*n_f,n_v)*Vxxx1]',n_x^4,n_f)]';
tempresult=tempresult*Matrix;
tempresult=[reshape([reshape([reshape(tempresult,n_f*n_x,n_x^3)]',n_x^3*n_f,n_x)*hxx]',n_x^5,n_f)]';
tempresult=tempresult*kron(Ix,OMEGA_x.OMEGA2);

result=result+tempresult;

clear tempresult Matrix

% (Vx1*zeta)(X)(Vxx1*(zeta(X)hxxx)*Omega3

Matrix=kron(Ezeta2,hxxx);
tempresult=[reshape([reshape([reshape(fvv,n_f*n_v,n_v)*Vx1]',n_x*n_f,n_v)*Vxx1]',n_x^3,n_f)]';
tempresult=tempresult*Matrix;

tempresult=tempresult*kron(Ix,OMEGA_x.OMEGA3);

result=result+tempresult;


