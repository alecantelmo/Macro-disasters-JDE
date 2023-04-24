%
% © Copyright, Oren Levintal, June 13, 2016.

order=4;

% Vx0(X)Vxxx0
derivs=[2,1,3];
perms1=[];
perms2=perms1;
Matrix=[];

result= permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx0,Vxxx0);

% Vx0(X)[Vxxx1*P(zeta^(X2)(X)hx)]

perms1=[2,1,3,4;2,3,1,4];
perms2=[1,2;1,2];
Matrix=kron(Ix,kron(Ezeta2,hx));

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx0,Vxxx1);

% Vx0(X)[Vxxx1*(zeta^(X3))]

perms1=[];
perms2=[];
Matrix=Ix_Ezeta3;

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx0,Vxxx1);

% (Vx1*zeta)(X)[Vxxx1*P(zeta(X)hx^(X2))]

perms1=[1,3,2,4;3,1,2,4];
perms2=[1,2;1,2];
Matrix=kron(Ezeta2,hx2);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx1,Vxxx1);

% (Vx1*zeta)(X)[Vxxx1*P(zeta^(X2)(X)hx)]

perms1=[2,1,3,4;2,3,1,4];
perms2=[1,2;1,2];
Matrix=kron(Ezeta3,hx);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx1,Vxxx1);

% (Vx1*zeta)(X)[Vxxx1*zeta^(X3)]

perms1=[];
perms2=[];
Matrix=Ezeta4;

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx1,Vxxx1);

% (Vx1*zeta)(X)[Vxx1*(zeta(X)hxx)Omega1]

Matrix=kron(Ezeta2,hxx);

tempresult=[reshape([reshape([reshape(fvv,n_f*n_v,n_v)*Vx1]',n_x*n_f,n_v)*Vxx1]',n_x^3,n_f)]'*Matrix;

tempresult=reshape(tempresult,n_f*n_x^3,n_x)';

tempresult=reshape(tempresult,n_x*n_f,n_x^3);
tempresult=tempresult*OMEGA_x.OMEGA1;
tempresult=reshape(tempresult,n_x,n_f*n_x^3)';
tempresult=reshape(tempresult,n_f,n_x^4);

result=result+tempresult;

clear tempresult 