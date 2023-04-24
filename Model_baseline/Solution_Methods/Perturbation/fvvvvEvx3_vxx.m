%
% © Copyright, Oren Levintal, June 13, 2016.

order=5;

% Vx0^(X3)(X)Vxx0

derivs=[4,1,1,1,2];
perms1=[];
perms2=perms1;
Matrix=[];

result= permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvv,Vx0,Vx0,Vx0,Vxx0);

% Vx0^(X3)(X)(Vxx1*zeta^(X2))

perms1=[];
perms2=perms1;
Matrix=Ix3_Ezeta2;

result= result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvv,Vx0,Vx0,Vx0,Vxx1);

% P(Vx0^(X2)(X)(Vx1*zeta)(X)(Vxx1*P(zeta(X)hx))

perms1=[1,2,4,3,5;1,2,4,5,3;2,1,3,4,5;2,1,4,3,5;2,1,4,5,3;];
perms2=perms1(:,[1,3,4,5]); perms2(perms2==2)=1; perms2(perms2>1)=perms2(perms2>1)-1;
Matrix=kron(Ix2,kron(Ezeta2,hx));

result= result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvv,Vx0,Vx0,Vx1,Vxx1);

% P(Vx0^(X2)(X)(Vx1*zeta)(X)(Vxx1*zeta^(X2))

perms1=[1,2,4,3,5;1,2,4,5,3];
perms2=perms1(:,[1,3,4,5]); perms2(perms2==2)=1; perms2(perms2>1)=perms2(perms2>1)-1;
Matrix=Ix2_Ezeta3;

result= result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvv,Vx0,Vx0,Vx1,Vxx1);

% P(Vx0(X)(Vx1*zeta^(X2)(X)Vxx0)

perms1=[1,2,3,5,4;1,2,5,3,4];
perms2=perms1(:,[1,3,4,5]); perms2(perms2==2)=1; perms2(perms2>1)=perms2(perms2>1)-1;
Matrix=kron(Ix,kron(Ezeta2,Ix2));

result= result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvv,Vx0,Vx1,Vx1,Vxx0);

% P(Vx0(X)(Vx1*zeta^(X2)(X)(Vxx1*P(zeta(X)hx))

perms1=[1,2,3,5,4;1,2,5,3,4;2,1,3,4,5;2,1,3,5,4;2,1,5,3,4];
perms2=perms1(:,[1,3,4,5]); perms2(perms2==2)=1; perms2(perms2>1)=perms2(perms2>1)-1;
Matrix=kron(Ix,kron(Ezeta3,hx));

result= result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvv,Vx0,Vx1,Vx1,Vxx1);

% P(Vx0(X)(Vx1*zeta^(X2)(X)(Vxx1*zeta^(X2))

perms1=[1,2,3,5,4;1,2,5,3,4];
perms2=perms1(:,[1,3,4,5]); perms2(perms2==2)=1; perms2(perms2>1)=perms2(perms2>1)-1;
Matrix=Ix_Ezeta4;

result= result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvv,Vx0,Vx1,Vx1,Vxx1);

% (Vx1*zeta)^(X3)(X)Vxx0

perms1=[];
perms2=[];
Matrix=kron(Ezeta3,Ix2);

result= result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvv,Vx1,Vx1,Vx1,Vxx0);

% (Vx1*zeta)^(X3)(X)(Vxx1*P(zeta(X)hx)

perms1=[2,1,3,4,5];
perms2=perms1(:,[1,3,4,5]); perms2(perms2==2)=1; perms2(perms2>1)=perms2(perms2>1)-1;
Matrix=kron(Ezeta4,hx);

result= result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvv,Vx1,Vx1,Vx1,Vxx1);

% (Vx1*zeta)^(X3)(X)(Vxx1*zeta2^(X2))

perms1=[];
perms2=[];
Matrix=Ezeta5;

result= result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvv,Vx1,Vx1,Vx1,Vxx1);
