%
% © Copyright, Oren Levintal, June 13, 2016.

order=5;

% Vx0^(X2)(X)Vxxx0

derivs=[3,1,1,3];
perms1=[];
perms2=perms1;
Matrix=[];

result= permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx0,Vx0,Vxxx0);

% Vx0^(X2)(X)(Vxxx1*P(zeta^(X2)(X)hx))

perms1=[2,1,3,4,5;2,3,1,4,5];
perms2=[1,2,3;1,2,3];
Matrix=kron(Ix2_Ezeta2,hx);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx0,Vx0,Vxxx1);

% Vx0^(X2)(X)(Vxxx1*zeta^(X3))

perms1=[];
perms2=[];
Matrix=Ix2_Ezeta3;

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx0,Vx0,Vxxx1);

% P(Vx0(X)(Vx1*zeta))(X)(Vxxx1*P(zeta(X)hx^(X2)))

perms1=[1,3,2,4,5;3,1,2,4,5;1,2,3,5,4;1,3,2,5,4;3,1,2,5,4];
perms2=perms1(:,[1,4,5]); perms2(perms2<=3)=1; perms2(perms2>1)=perms2(perms2>1)-2;
Matrix=kron(Ix_Ezeta2,hx2);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx0,Vx1,Vxxx1);

% P(Vx0(X)(Vx1*zeta))(X)(Vxxx1*P(zeta^(X2)(X)hx)

perms1=[2,1,3,4,5;2,3,1,4,5;1,2,3,5,4;2,1,3,5,4;2,3,1,5,4];
perms2=perms1(:,[1,4,5]); perms2(perms2<=3)=1; perms2(perms2>1)=perms2(perms2>1)-2;
Matrix=kron(Ix_Ezeta3,hx);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx0,Vx1,Vxxx1);

% P(Vx0(X)(Vx1*zeta))(X)(Vxxx1*zeta^(X3))

perms1=[1,2,3,5,4];
perms2=perms1(:,[1,4,5]); perms2(perms2<=3)=1; perms2(perms2>1)=perms2(perms2>1)-2;
Matrix=Ix_Ezeta4;

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx0,Vx1,Vxxx1);

% P(Vx0(X)(Vx1*zeta))(X)(Vxx1*(zeta(X)hxx)*Omega1

perms1=[1,2,4,3];
perms2=perms1(:,[1,3,4]); perms2(perms2<=2)=1; perms2(perms2>1)=perms2(perms2>1)-1;
Matrix=kron(Ix_Ezeta2,Ix);

tempresult=permutekron2([3,1,1,2],4,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx0,Vx1,Vxx1);

tempmat_Omega1=reshape(kron(Ix,hxx),n_x^2,n_x^3)*OMEGA_x.OMEGA1;

clear tempmat

tempmat_Omega1=reshape(tempmat_Omega1,n_x^2,n_x^3);

tempresult=reshape([reshape(tempresult,n_f*n_x^2,n_x^2)]',n_x^2*n_f,n_x^2);
tempresult=[reshape([tempresult*tempmat_Omega1]',n_x^5,n_f)]';

result=result+tempresult;

clear tempresult 

% (Vx1*zeta)^(X2)(X)Vxxx0

perms1=[];
perms2=[];
Matrix=kron(Ezeta2,Ix3);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx1,Vx1,Vxxx0);

% (Vx1*zeta)^(X2)(X)(Vxxx1*P(zeta(X)hx^(X2))

perms1=[1,3,2,4,5;3,1,2,4,5];
perms2=[1,2,3;1,2,3];
Matrix=kron(Ezeta3,hx2);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx1,Vx1,Vxxx1);
% (Vx1*zeta)^(X2)(X)(Vxxx1*P(zeta^(X2)(X)hx)

perms1=[2,1,3,4,5;2,3,1,4,5];
perms2=[1,2,3;1,2,3];
Matrix=kron(Ezeta4,hx);

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx1,Vx1,Vxxx1);

% (Vx1*zeta)^(X2)(X)(Vxxx1*zeta^(X3))

perms1=[];
perms2=[];
Matrix=Ezeta5;

result=result+permutekron2(derivs,order,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx1,Vx1,Vxxx1);

% (Vx1*zeta)^(X2)(X)(Vxx1*(zeta(X)hxx)*Omega1

perms1=[];
perms2=[];
Matrix=kron(Ezeta3,Ix);

tempresult=permutekron2([3,1,1,2],4,perms1,perms2,n_v,n_x,n_f,Matrix,fvvv,Vx1,Vx1,Vxx1);

tempresult=reshape([reshape(tempresult,n_f*n_x^2,n_x^2)]',n_x^2*n_f,n_x^2);
tempresult=[reshape([tempresult*tempmat_Omega1]',n_x^5,n_f)]';

result=result+tempresult;

clear tempresult tempmat_Omega1