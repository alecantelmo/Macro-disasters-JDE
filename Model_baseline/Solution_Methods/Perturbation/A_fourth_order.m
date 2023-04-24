%
% © Copyright, Oren Levintal, June 13, 2016.

fvvvvEvx4;
A=result;


fvvvEvx2_vxx;

result=result+fv*Vxxx1*kron(Ezeta2,hxx);

A=A+result*OMEGA_x.OMEGA2;

fvvEvx_vxxx;

A=A+result*OMEGA_x.OMEGA3;
saveA=A;

fvvEvxx2;

A=A+result*OMEGA_x.OMEGA4;

A=A+fyp*chain4(gx,gxx,gxxx,[],hx,hxx,hxxx,[],OMEGA_x.OMEGA2,OMEGA_x.OMEGA3,OMEGA_x.OMEGA4);

