function M=get_moments(nep,P,order)
%
% © Copyright, Oren Levintal, June 13, 2016.

order=order(1);
M2=0;M3=0;M4=0;M5=0;
for i=1:length(P)
    if order>=2
        nep2=kron(nep(:,i),nep(:,i));
        M2=M2+nep2*P(i);
    end
    if order>=3
        nep3=kron(nep2,nep(:,i));
        M3=M3+nep3*P(i);
    end
    if order>=4
        nep4=kron(nep3,nep(:,i));
        M4=M4+nep4*P(i);
    end
    if order>=5
        nep5=kron(nep4,nep(:,i));
        M5=M5+nep5*P(i);
    end
end

if order==1
    M.M1=sparse(size(nep,1),1);
end
if order>=2
    M2(abs(M2)<eps)=0;
    M.M2=sparse(M2);
end
if order>=3
    M3(abs(M3)<eps)=0;
    M.M3=sparse(M3);
end
if order>=4
    M4(abs(M4)<eps)=0;
    M.M4=sparse(M4);
end
if order>=5
    M5(abs(M5)<eps)=0;
    M.M5=sparse(M5);
end

end