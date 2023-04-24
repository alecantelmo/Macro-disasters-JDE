%
% © Copyright, Oren Levintal, June 13, 2016.

order=3;

% Vx0^(X3)

derivs=[3,1,1,1];
perms1=[];
perms2=perms1;
Matrix=[];

result= permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx0,Vx0,Vx0);

% P(Vx0(X)[Vx1*zeta)^(X2)]

perms1=[1,3,2;3,1,2];
perms2=perms1;
Matrix=Ix_Ezeta2;

result= result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx0,Vx1,Vx1);

% (Vx1*zeta)^(X3)

perms1=[];
perms2=perms1;
Matrix=Ezeta3;

result= result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx1,Vx1,Vx1);




