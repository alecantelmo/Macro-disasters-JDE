function R=tkron3_x(A,B,C,Ax,Bx,Cx)
% Compute first derivative of kron(A,B,C) w.r.t x
% A is m1-by-1
% Ax is m1-by-n_x
% B is m2-by-1
% Bx is m2-by-n_x
% C is m3-by-1
% Cx is m3-by-n_x
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

term1=tkron3(Ax,B,C); %m3*m2*m1,n_x

term2=tkron3(A,Bx,C); %m3*m2*m1,n_x

term3=trkron3(A,B,Cx); %m3*m2*m1,n_x

R=tplus(term1,term2);

R=tplus(R,term3);

end