%
% © Copyright, Oren Levintal, June 13, 2016.

order=3;

% Vx0(X)Vxx0
derivs=[2,1,2];
perms1=[];
perms2=perms1;
Matrix=[];

result= permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx0,Vxx0);

% (Vx1*zeta)(X)[Vxx1*P(zeta(X)hx)]
perms1=[2,1,3];
perms2=[1,2];
Matrix=kron(Ezeta2,hx);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx1,Vxx1);

% Vx0(X)[Vxx1*(zeta^(X2))]
perms1=[];
perms2=[];
Matrix=Ix_Ezeta2;

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx0,Vxx1);

% (Vx1*zeta)(X)[Vxx1*zeta^(X2)]
perms1=[];
perms2=[];
Matrix=Ezeta3;

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vx1,Vxx1);
