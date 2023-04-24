%
% © Copyright, Oren Levintal, June 13, 2016.

order=4;

% Vx0^(X2)(X)Vxx0

derivs=[3,1,1,2];
perms1=[];
perms2=perms1;
Matrix=[];

result= permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx0,Vx0,Vxx0);

% (Vx1*zeta)^(X2)(X)Vxx0

perms1=[];
perms2=[];
Matrix=kron(Ezeta2,Ix2);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx1,Vx1,Vxx0);

% P(Vx0(X)(Vx1*zeta))(X)(Vxx1*P(zeta(X)hx))

perms1=[2,1,3,4;1,2,4,3;2,1,4,3];
perms2=perms1(:,[1,3,4]); perms2(perms2<=2)=1; perms2(perms2>1)=perms2(perms2>1)-1;
Matrix=kron(Ix_Ezeta2,hx);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx0,Vx1,Vxx1);

% (Vx1*zeta)^(X2)(X)(Vxx1*P(zeta(X)hx)

perms1=[2,1,3,4];
perms2=[1,2,3];
Matrix=kron(Ezeta3,hx);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx1,Vx1,Vxx1);

% Vx0^(X2)(X)(Vxx1*zeta^(X2))

perms1=[];
perms2=[];
Matrix=Ix2_Ezeta2;

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx0,Vx0,Vxx1);

% P(Vx0(X)(Vx1*zeta))(X)(Vxx1*zeta^(X2))

perms1=[1,2,4,3];
perms2=[1,3,2];
Matrix=Ix_Ezeta3;

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx0,Vx1,Vxx1);

% (Vx1*zeta)^(X2)(X)(Vxx1*zeta^(X2))

perms1=[];
perms2=[];
Matrix=Ezeta4;

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx1,Vx1,Vxx1);
