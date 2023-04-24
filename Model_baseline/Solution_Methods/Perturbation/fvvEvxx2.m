%
% © Copyright, Oren Levintal, June 13, 2016.

order=4;

% Vxx0^(X2)

derivs=[2,2,2];
perms1=[];
perms2=perms1;
Matrix=[];

result= permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vxx0,Vxx0);

% P(Vxx0(X)(Vxx1*zeta^(X2)))

perms1=[3,4,1,2];
perms2=[2,1];
Matrix=Ix2_Ezeta2;

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vxx0,Vxx1);

% (Vxx1*P(hx(X)zeta))^(X2)

perms1=[2,1,3,4;1,2,4,3;2,1,4,3];
perms2=[1,2;1,2;1,2];
Matrix=hx_Ezeta2_hx;

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vxx1,Vxx1);
% P((Vxx1*P(hx(X)zeta))(X)(Vxx1*zeta^(X2)))

perms1=[1,2,4,3;3,4,1,2;4,3,1,2];
perms2=[1,2;1,2;1,2];
Matrix=hx_Ezeta3;

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vxx1,Vxx1);

% (Vxx1*zeta^(X2))(X)(Vxx1*zeta^(X2))

perms1=[];
perms2=[];
Matrix=Ezeta4;

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvv,Vxx1,Vxx1);

