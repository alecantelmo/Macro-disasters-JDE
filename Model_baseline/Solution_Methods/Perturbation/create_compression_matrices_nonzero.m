%
% © Copyright, Oren Levintal, June 13, 2016.

if approx>=2
    k=2;
    nonzero=reshape(ones(1,n_x^k),n_x,n_x);
    nonzero(end,1:end-1)=0;
    nonzero(1:end-1,end)=0;
    [U2,W2]=create_UW(n_x,k,nonzero(:));
    nnz=find(nonzero);
    N2=sparse(nnz,1:length(nnz),ones(1,length(nnz)),n_x^k,length(nnz));
    U2=N2*U2;
    W2=W2*N2';
end

if approx>=3
    k=3;
    nonzero=reshape(ones(1,n_x^k),n_x,n_x,n_x);
    nonzero(end,1:end-1,1:end-1)=0;
    nonzero(1:end-1,end,1:end-1)=0;
    nonzero(1:end-1,1:end-1,end)=0;
    [U3,W3]=create_UW(n_x,k,nonzero(:));
    nnz=find(nonzero);
    N3=sparse(nnz,1:length(nnz),ones(1,length(nnz)),n_x^k,length(nnz));
    U3=N3*U3;
    W3=W3*N3';
end

if approx>=4
    k=4;
    nonzero=reshape(ones(1,n_x^k),n_x,n_x,n_x,n_x);
    nonzero(end,1:end-1,1:end-1,1:end-1)=0;
    nonzero(1:end-1,end,1:end-1,1:end-1)=0;
    nonzero(1:end-1,1:end-1,end,1:end-1)=0;
    nonzero(1:end-1,1:end-1,1:end-1,end)=0;
    [U4,W4]=create_UW(n_x,k,nonzero(:));
    nnz=find(nonzero);
    N4=sparse(nnz,1:length(nnz),ones(1,length(nnz)),n_x^k,length(nnz));
    U4=N4*U4;
    W4=W4*N4';
end

if approx>=5
    k=5;
    nonzero=reshape(ones(1,n_x^k),n_x,n_x,n_x,n_x,n_x);
    nonzero(end,1:end-1,1:end-1,1:end-1,1:end-1)=0;
    nonzero(1:end-1,end,1:end-1,1:end-1,1:end-1)=0;
    nonzero(1:end-1,1:end-1,end,1:end-1,1:end-1)=0;
    nonzero(1:end-1,1:end-1,1:end-1,end,1:end-1)=0;
    nonzero(1:end-1,1:end-1,1:end-1,1:end-1,end)=0;
    [U5,W5]=create_UW(n_x,k,nonzero(:));
    nnz=find(nonzero);
    N5=sparse(nnz,1:length(nnz),ones(1,length(nnz)),n_x^k,length(nnz));
    U5=N5*U5;
    W5=W5*N5';
end

if approx==1
    UW=struct([]);
elseif approx==2
    UW.U2=U2; UW.W2=W2;
elseif approx==3
    UW.U2=U2; UW.W2=W2;
    UW.U3=U3; UW.W3=W3;
elseif approx==4
    UW.U2=U2; UW.W2=W2;
    UW.U3=U3; UW.W3=W3;
    UW.U4=U4; UW.W4=W4;
elseif approx==5
    UW.U2=U2; UW.W2=W2;
    UW.U3=U3; UW.W3=W3;
    UW.U4=U4; UW.W4=W4;
    UW.U5=U5; UW.W5=W5;
end

