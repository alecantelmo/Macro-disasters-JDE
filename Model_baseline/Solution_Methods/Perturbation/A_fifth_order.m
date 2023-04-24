%
% © Copyright, Oren Levintal, June 13, 2016.

fvvvvvEvx5;

A=result;


fvvvvEvx3_vxx;


derivs=[1,4];
perms1=[1,2,4,3;1,4,2,3];
perms2=[1;1];
Matrix=kron(kron(hx,Ezeta2),Ix);
tempresult= permutekron2(derivs,4,perms1,perms2,n_v,n_x,n_f,Matrix,fv,Vxxxx1);

perms1=[];
perms2=[];
Matrix=kron(Ezeta3,Ix);
tempresult=tempresult+permutekron2(derivs,4,perms1,perms2,n_v,n_x,n_f,Matrix,fv,Vxxxx1);

result=result+tempresult*kron(Ix3,hxx);
clear tempresult Matrix

A=A+result*OMEGA_x.OMEGA5;


fvvvEvx2_vxxx;


derivs=[1,3];
struct=[0,0,1];
perms1=[];
perms2=[];
Matrix=kron(Ezeta2,Ix);

tempresult= permutekron2(derivs,3,perms1,perms2,n_v,n_x,n_f,Matrix,fv,Vxxx1);

result=result+tempresult*kron(Ix2,hxxx);
clear tempresult Matrix

A=A+result*OMEGA_x.OMEGA6;


fvvvEvx_vxx2;


A=A+result*OMEGA_x.OMEGA7;


fvvEvx_vxxxx;


A=A+result*OMEGA_x.OMEGA8;


fvvEvxx_vxxx;

A=A+result*OMEGA_x.OMEGA9;

clear result

A=reshape(A,n_f,n_x^5);

A=A+fyp*reshape(chain5(gx,gxx,gxxx,gxxxx,[],hx,hxx,hxxx,hxxxx,[],OMEGA_x.OMEGA5,OMEGA_x.OMEGA6,OMEGA_x.OMEGA7,...
    OMEGA_x.OMEGA8,OMEGA_x.OMEGA9),n_y,n_x^5);
    
clear result