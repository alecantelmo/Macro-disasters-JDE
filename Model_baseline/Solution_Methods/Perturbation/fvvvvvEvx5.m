%
% © Copyright, Oren Levintal, June 13, 2016.

order=5;

% Vx0^(X5)

derivs=[5,1,1,1,1,1];
perms1=[];
perms2=perms1;
Matrix=[];

result= permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvvv,Vx0,Vx0,Vx0,Vx0,Vx0);

% P(Vx0(X)[Vx1*zeta)^(X4)]

perms1=[1,2,3,5,4;1,2,5,3,4;1,5,2,3,4;5,1,2,3,4];
perms2=perms1;
Matrix=Ix_Ezeta4;

result= result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvvv,Vx0,Vx1,Vx1,Vx1,Vx1);

% P(Vx0^(X2)(X)(Vx1*zeta)^(X3))

perms1=[1,2,4,3,5;1,2,4,5,3;1,4,2,3,5;1,4,2,5,3;1,4,5,2,3;4,1,2,3,5;4,1,2,5,3;4,1,5,2,3;4,5,1,2,3];
perms2=perms1;
Matrix=Ix2_Ezeta3;

result= result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvvv,Vx0,Vx0,Vx1,Vx1,Vx1);

% P(Vx0^(X3)(X)(Vx1*zeta)^(X2))

perms1=[1,3,2,4,5;1,3,4,2,5;1,3,4,5,2;3,1,2,4,5;3,1,4,2,5;3,1,4,5,2;3,4,1,2,5;3,4,1,5,2;3,4,5,1,2];
perms2=perms1;
Matrix=Ix3_Ezeta2;

result= result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvvv,Vx0,Vx0,Vx0,Vx1,Vx1);

% (Vx1*zeta)^(X5)

perms1=[];
perms2=perms1;
Matrix=Ezeta5;

result= result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvvv,Vx1,Vx1,Vx1,Vx1,Vx1);

