n_e=size(eta_mat,2);

M2=0;M3=0;M4=0;M5=0;
for i=1:length(P)
    nep2=kron(nep(:,i),nep(:,i));
    nep3=kron(nep2,nep(:,i));
    nep4=kron(nep3,nep(:,i));
    nep5=kron(nep4,nep(:,i));
    M2=M2+nep2*P(i);
    M3=M3+nep3*P(i);
    M4=M4+nep4*P(i);
    M5=M5+nep5*P(i);
end


M2(abs(M2)<eps)=0;
M3(abs(M3)<eps)=0;
M4(abs(M4)<eps)=0;
M5(abs(M5)<eps)=0;

M.M2=sparse(M2);
M.M3=sparse(M3);
M.M4=sparse(M4);
M.M5=sparse(M5);

