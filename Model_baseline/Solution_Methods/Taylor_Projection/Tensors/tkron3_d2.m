function R=tkron3_d2(A,B,C,Ax,Bx,Cx,Axx,Bxx,Cxx)
% Compute first derivative of kron(A,B,C)
% A is m1-by-1
% Ax is m1-by-n_x
% Axx is m1-by-n_x-by-n_x
% B is m2-by-1
% Bx is m2-by-n_x
% Bxx is m2-by-n_x-by-n_x
% C is m3-by-1
% Cx is m3-by-n_x
% Cxx is m3-by-n_x-by-n_x
%
% © Copyright, Oren Levintal, June 13, 2016.

m1=A.tsize(1);
m2=B.tsize(1);
m3=C.tsize(1);

if m1~=prod(A.tsize)
    error('A must be a vector')
end
if m2~=prod(B.tsize)
    error('B must be a vector')
end
if m3~=prod(C.tsize)
    error('C must be a vector')
end

if Ax.tsize(1)~=m1
    error('Incompatible dimensions')
end
if Bx.tsize(1)~=m2
    error('Incompatible dimensions')
end
if Cx.tsize(1)~=m3
    error('Incompatible dimensions')
end

if Ax.tsize(2)~=Bx.tsize(2)
    error('Incompatible dimensions')
elseif Ax.tsize(2)~=Cx.tsize(2)
    error('Incompatible dimensions')
end

Axx=ptr1d(col2ptr(Axx,1));
Bxx=ptr1d(col2ptr(Bxx,1));
Cxx=ptr1d(col2ptr(Cxx,1));


term1=tkron3_d1(tvec(Ax),B,C,Axx,Bx,Cx); %m3*m2*m1*n_x,n_x

term1=ptr2d(term1,m3*m2*m1,n_x);
term1=ptr2col(term1,1); %m3*m2*m1,n_x,n_x


term2=tkron3_d1(A,tvec(Bx),C,Ax,Bxx,Cx); %m3*m2*n_x*m1,n_x

term2=ptr2d(term2,m3*m2,n_x*m1);
term2=ptr2col(term2,1); %m3*m2,n_x*m1,n_x
term2=fold(term2,n_x,m1,n_x);%m3*m2,n_x,m1,n_x
term2=col2ptr(term2,2);
term2=ptr1d(term2);%m3*m2*m1,n_x,n_x

term3=trkron3_d1(A,B,Cx,Ax,Bx,Cxx); %m3*n_x*m2*m1,n_x

term3=ptr2d(term3,m3,n_x*m2*m1);
term3=ptr2col(term3,1); %m3,n_x*m2*m1,n_x
term3=fold(term3,n_x,m2*m1,n_x);%m3,n_x,m2*m1,n_x
term3=col2ptr(term3,2);
term3=ptr1d(term3);%m3*m2*m1,n_x,n_x

R=tplus(term1,term2);

R=tplus(R,term3);

end