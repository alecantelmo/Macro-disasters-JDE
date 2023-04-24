%
% © Copyright, Oren Levintal, June 13, 2016.

order=5;

% Vx0(X)Vxx0^(X2)
derivs=[3,1,2,2];

perms1=[];
perms2=perms1;
Matrix=[];

result= permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx0,Vxx0,Vxx0);

% Vx0(X)P(Vxx0(X)(Vxx1*zeta^(X2))

perms1=[3,4,1,2,5];
perms2=[2,1,3];
Matrix=Ix3_Ezeta2;

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx0,Vxx0,Vxx1);

% Vx0(X)P(Vxx1*P(hx(X)zeta))(X)(Vxx1*zeta^(X2))

perms1=[1,2,4,3,5;3,4,1,2,5;4,3,1,2,5];
perms2=[1,2,3;1,2,3;1,2,3];
Matrix=kron(kron(Ix,hx),Ezeta3);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx0,Vxx1,Vxx1);

% Vx0(X)(Vxx1*P(hx(X)zeta))(X)(Vxx1*P(zeta(X)hx))

perms1=[2,1,3,4,5;1,2,4,3,5;2,1,4,3,5];
perms2=[1,2,3;1,2,3;1,2,3];
Matrix=kron(kron(Ix,hx),kron(Ezeta2,hx));

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx0,Vxx1,Vxx1);

% Vx0(X)(Vxx1*zeta^(X2)(X)(Vxx1*zeta^(X2))

perms1=[];
perms2=[];
Matrix=kron(Ix,Ezeta4);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx0,Vxx1,Vxx1);

% (Vx1*zeta)(X)P((Vxx1*P(zeta(X)hx))(X)Vxx0)

perms1=[1,2,4,3,5;3,4,1,2,5;4,3,1,2,5];
perms2=[1,2,3;2,1,3;2,1,3];
Matrix=kron(kron(Ezeta2,hx),Ix2);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx1,Vxx1,Vxx0);

% (Vx1*zeta)(X)P((Vxx1*zeta^(X2)(X)Vxx0)

perms1=[3,4,1,2,5];
perms2=[2,1,3];
Matrix=kron(Ezeta3,Ix2);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx1,Vxx1,Vxx0);

% (Vx1*zeta)(X)P((Vxx1*zeta^(X2))(X)(Vxx1*P(zeta(X)hx)))

perms1=[2,1,3,4,5;3,4,1,2,5;3,4,2,1,5];
perms2=[1,2,3;2,1,3;2,1,3];
Matrix=kron(Ezeta4,hx);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx1,Vxx1,Vxx1);

% (Vx1*zeta)(X)(Vxx1*P(zeta(X)hx))^(X2): this one is different, because the structrue of the
% unpermuted Matrix does not have all zeta's in a row.

% To calculate expected value of the stochastic Matrix, use
% kron(Ix,Ix,zeta,zeta,zeta) and permute it.
Matrix=kron(Ezeta3,hx2);
tempindex=permute(reshape([1:n_x^5],n_x,n_x,n_x,n_x,n_x),[1,3,2,4,5]);
Matrix=Matrix(tempindex,tempindex);
clear tempindex

perms1=[2,1,3,4,5;1,2,4,3,5;2,1,4,3,5;];
perms2=[1,2,3;1,2,3;1,2,3];

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx1,Vxx1,Vxx1);

% (Vx1*zeta)(X)(Vxx1*zeta^(X2))^(X2)

perms1=[];
perms2=[];
Matrix=Ezeta5;

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx1,Vxx1,Vxx1);
