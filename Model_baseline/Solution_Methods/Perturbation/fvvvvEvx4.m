%
% © Copyright, Oren Levintal, June 13, 2016.

order=4;

% Vx0^(X4)

derivs=[4,1,1,1,1];
perms1=[];
perms2=perms1;
Matrix=[];

result= permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvv,Vx0,Vx0,Vx0,Vx0);

% P(Vx0(X)[Vx1*zeta)^(X3)]

perms1=[1,2,4,3;1,4,2,3;4,1,2,3];
perms2=perms1;
Matrix=Ix_Ezeta3;

result= result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvv,Vx0,Vx1,Vx1,Vx1);

% P(Vx0^(X2)(X)(Vx1*zeta)^(X2))

perms1=[1,3,2,4;1,3,4,2;3,1,2,4;3,1,4,2;3,4,1,2];
perms2=perms1;
Matrix=Ix2_Ezeta2;

result= result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvv,Vx0,Vx0,Vx1,Vx1);

% (Vx1*zeta)^(X4)

perms1=[];
perms2=perms1;
Matrix=Ezeta4;

result= result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvvv,Vx1,Vx1,Vx1,Vx1);

