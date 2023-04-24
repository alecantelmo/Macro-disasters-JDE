function [derivs,sol_time,deriv_time]=solve_dsge(model,params,M,eta,nxss,nyss,approx,algo,varargin)
%solve_dsge(model,params,M,eta,nxss,nyss,approx,algo,gx,hx)
%This function solves the dsge model with perturbation up to a fifth
%order.
%Input variables:
%model: a structure that is generated automatically by the function
%differentiate_dsge.m
%params: a vector of all parameter values ordered in the same order of symparams.
%M: a structure that contains all the cross moments of the shocks. The fields of
%this structure should be M2,M3,M4,M5.
%eta: The matrix eta as defined in Schmitt-Grohe and Uribe (2004).
%nxss,nyss: steady state values of x and y
%algo: algorithm type. 'dlyap' for dlyap or 'vectorize' for vectorization.
%gx,hx: first order solutions. These arguments are optional. If not supplied, the first order solution is calculated
%by the function gx_hx.m written by Schmitt-Grohe and Uribe (2004).
%
%
% © Copyright, Oren Levintal, June 13, 2016.

if strcmpi(algo,'binning')
    kamenik_type=1;
    algo='Kamenik';
elseif strcmpi(algo,'gensylv')
    kamenik_type=2;
    algo='Kamenik';
elseif ~strcmpi(algo,'dlyap') && ~strcmpi(algo,'vectorize')
    error('wrong algorithm')
end

sol_time=0;

n_f=model.n_f;
n_x=model.n_x;
n_y=model.n_y;
n_x2=model.n_x2;
n_x1=model.n_x1;
n_v=model.n_v;

f_ind=model.f_ind;
UW=model.UW; % IMPORTANT: THESE COMPRESSION MATRICES ASSUME A CERTAIN NONZERO PATTERN.
OMEGA_x=model.OMEGA_x;

if isempty(OMEGA_x)
    create_OMEGA_x;
end

n_e=size(eta,2);
eta=[eta;zeros(1,n_e)];

nxss=[nxss(:);0]; % add steady state value of the perturbation variable.
nyss=nyss(:);

nPhi=Phi_fun(nxss,params);
if approx>=1
    nPhix=Phi_d1(nxss',params,model.Phi_ind);
end
if approx>=2
    W2=UW.W2; U2=UW.U2;
    nPhixx=Phi_d2(nxss',params,model.Phi_ind);
end
if approx>=3
    W3=UW.W3; U3=UW.U3;
    nPhixxx=Phi_d3(nxss',params,model.Phi_ind);
end
if approx>=4
    W4=UW.W4; U4=UW.U4;
    nPhixxxx=Phi_d4(nxss',params,model.Phi_ind);
end
if approx>=5
    W5=UW.W5; U5=UW.U5;
    nPhixxxxx=Phi_d5(nxss',params,model.Phi_ind);
end

nv=[nyss;nyss;nxss;nxss];

start=tic;
fv=f_d1(nv',params,model.f_ind);
deriv_time=toc(start);
fv=reshape(fv,n_f,n_v);


fyp=fv(:,1:n_y); fy=fv(:,n_y+1:2*n_y); fxp=fv(:,2*n_y+1:2*n_y+n_x); fx=fv(:,2*n_y+n_x+1:end);

% First Order

% If first order solution not provided, use the function gx_hx of Schmitt-Grohe and Uribe (2004)
if isempty(varargin)
    tic
    [gx,hx,exitflag]=gx_hx(full([fy;zeros(n_x2,n_y)]),full([fx(:,1:end-1);nPhix(:,1:end-1)]),full([fyp;zeros(n_x2,n_y)]),full([fxp(:,1:end-1);[zeros(n_x2,n_x1),-eye(n_x2)]]));
    time=toc;
    sol_time=sol_time+time;
    gx=[gx,zeros(n_y,1)];
    hx=[hx;zeros(1,n_x-1)];
    hx=[hx,zeros(n_x,1)];
    hx(end,end)=1;
    % replace relevant rows of hx with nPhix which is more numerically accurate
    hx(n_x1+1:n_x1+n_x2,:)=nPhix;
else
    gx=zeros(n_y,n_x);
    gx(1:size(varargin{1},1),1:size(varargin{1},2))=varargin{1};
    hx=zeros(n_x,n_x);
    hx(end,end)=1;
    hx(1:size(varargin{2},1),1:size(varargin{2},2))=varargin{2};
end
    
hx=sparse(hx);
gx=sparse(gx);

% Second Order
if approx>=2
    M2=M.M2;
    eta2_M2=reshape([eta*reshape(M2,n_e,n_e)]',n_e,n_x);
    eta2_M2=reshape([eta*eta2_M2]',n_x^2,1);

    start=tic;
    fvv=f_d2(nv',params,model.f_ind);
    deriv_time=deriv_time+toc(start);

    unique=nchoosek(n_x+2-1,2);
    unique=unique-nchoosek(n_x-1+1-1,1);

    Vx0=[gx*hx;gx;hx;speye(n_x)];
    Vx1=[gx;sparse(n_y,n_x);speye(n_x,n_x);sparse(n_x,n_x)];

    Ezeta2=[ sparse(n_x^2,n_x^2-1) , eta2_M2 ];

    A=innerkron(n_f,n_v,fvv,Vx0,Vx0)+innerkron(n_f,n_v,fvv,Vx1,Vx1)*Ezeta2;

    fy_fxp_fypgx=[fv(:,n_y+1:2*n_y) fv(:,2*n_y+1:2*n_y+n_x)+fv(:,1:n_y)*gx];

    G=fy_fxp_fypgx(:,1:n_f);
    H=fy_fxp_fypgx(:,n_f+1:n_f+n_x2);

    D=sparse(n_f,n_f);
    D(:,1:n_y)=fv(:,1:n_y);

    if strcmp(algo,'Kamenik')
        if n_x2==0
            C=A;
        else % do not solve exogenous state variables (see appendix A.6 in the paper). H is the last block of G.
            C=A+H*(nPhixx);
        end
        %Block 1: xx
        spx=sparse([ones(n_x-1,1);0]);
        choosex2=kron(spx,spx);
        choosex2=logical(choosex2==1);
        tempeye=speye(n_x^2);
        Z=tempeye(:,choosex2);
        CZ=C*Z;
        hx_hat=hx(1:end-1,1:end-1);
        if kamenik_type==1
            tic
            Xtemp=kamenik( G,D,hx_hat,-CZ,2 );
            time=toc;
        else
            G=full(G);
            D=full(D);
            hx_hat=full(hx_hat);
            tic
            [~,Xtemp]=gensylv( 2,G,D,hx_hat,full(-CZ) );
            time=toc;
        end
        sol_time=sol_time+time;
        acc=norm(full(CZ+AkronkC(D*Xtemp,hx_hat,2)+G*Xtemp));
        if acc>1e-8
            warning(['Second order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        X=zeros(n_f,n_x^2);
        X(:,choosex2)=full(Xtemp); 
        gxx_hat=Xtemp(1:n_y,:);
        hat_eta2_M2=eta2_M2(choosex2,:);
        %Block 2: ss
        sps=sparse([zeros(n_x-1,1);1]);
        choosex2=kron(sps,sps);
        choosex2=logical(choosex2==1);
        Z=tempeye(:,choosex2);
        CZ=C*Z;
        tic
        Xtemp=-(D+G)\(CZ+(fyp*gxx_hat)*hat_eta2_M2);
        time=toc;
        sol_time=sol_time+time;
        
        acc=norm(full(CZ+(fyp*gxx_hat)*hat_eta2_M2+(D+G)*Xtemp));
        if acc>1e-8
            warning(['Second order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        Xtemp=reshape(Xtemp,n_f,size(Z,2));
        X(:,choosex2)=full(Xtemp); clear Xtemp
    else
        if n_x2==0
            CU2=A*U2;
        else % do not solve exogenous state variables (see appendix A.6 in the paper). H is the last block of G.
            CU2=A*U2+H*(nPhixx*U2);
        end
        W2BU2=(innerkron(unique,n_x,W2,hx,hx)+W2*Ezeta2)*U2;
        %Block 1: xx
        spx=sparse([ones(n_x-1,1);0]);
        choosex2=kron(spx,spx);
        choosex2U=logical(U2'*choosex2==1);
        tempeye=speye(size(U2,2));
        Z=tempeye(:,choosex2U);
        CU2Z=CU2*Z;
        ZTW2BU2Z=Z'*W2BU2*Z;
        if strcmpi(algo,'vectorize')
            Xtemp=reshape(-(kron(ZTW2BU2Z',D)+kron(speye(size(Z,2)),G))\CU2Z(:),n_f,size(Z,2));
        elseif strcmpi(algo,'dlyap')
            Xtemp=dlyap(ZTW2BU2Z',(-G\D)',(-G\CU2Z)');
            Xtemp=Xtemp';
        elseif strcmpi(algo,'slicot')
            Xtemp=HessSchur(full(ZTW2BU2Z'),full((G\D)'),full((-G\CU2Z)'));
            Xtemp=Xtemp';
        end
        acc=norm(full(G*Xtemp+D*Xtemp*ZTW2BU2Z+CU2Z));
        if acc>1e-8
            warning(['Second order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        X=zeros(n_f,size(U2,2));
        X(:,choosex2U)=full(Xtemp); clear Xtemp ZTW2BU2Z
        %Block 2
        sps=sparse([zeros(n_x-1,1);1]);
        choosex2=kron(sps,sps);
        choosex2U=logical(U2'*choosex2==1);
        tempeye=speye(size(U2,2));
        Z=tempeye(:,choosex2U);
        W2BU2Z=W2BU2*Z;
        CU2Z=CU2*Z+D*X*W2BU2Z;
        ZTW2BU2Z=Z'*W2BU2*Z;
        if strcmpi(algo,'vectorize')
            Xtemp=reshape(-(kron(ZTW2BU2Z',D)+kron(speye(size(Z,2)),G))\CU2Z(:),n_f,size(Z,2));
        elseif strcmpi(algo,'dlyap')
            Xtemp=dlyap(ZTW2BU2Z',(-G\D)',(-G\CU2Z)');
            Xtemp=Xtemp';
        elseif strcmpi(algo,'slicot')
            Xtemp=HessSchur(full(ZTW2BU2Z'),full((G\D)'),full((-G\CU2Z)'));
            Xtemp=Xtemp';
        end
        acc=norm(G*Xtemp+D*Xtemp*ZTW2BU2Z+CU2Z);
        if acc>1e-8
            warning(['Second order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        X(:,choosex2U)=full(Xtemp); clear Xtemp W2BU2Z ZTW2BU2Z
        X=X*W2;
    end
    gxx=X(1:n_y,:);
    hxx=[X(n_y+1:end,:);nPhixx;zeros(1,n_x^2)];
end

% Third Order
if approx>=3
    Vxx0=[chain2(gx,gxx,hx,hxx);gxx;hxx;sparse(n_x,n_x^2)];
    Vxx1=[gxx;sparse(n_y+2*n_x,n_x^2)];
    
    M3=M.M3(:);

    eta3_M3=reshape([eta*reshape(M3,n_e,n_e^2)]',n_e,n_e*n_x);
    eta3_M3=reshape([eta*eta3_M3]',n_e,n_x^2);
    eta3_M3=reshape([eta*eta3_M3]',n_x^3,1);

    Ezeta3=[ sparse(n_x^3,n_x^3-1) , eta3_M3 ];
    Ix=speye(n_x);
    Ix_Ezeta2=kron(Ix,Ezeta2);
    
    start=tic;
    fvvv=f_d3(nv',params,model.f_ind);
    deriv_time=deriv_time+toc(start);
    
    unique=nchoosek(n_x+3-1,3);
    unique=unique-nchoosek(n_x-1+2-1,2);
    A_third_order; % create matrix A of a third order solution
    if strcmp(algo,'Kamenik')
        if n_x2==0
            C=A;
        else % do not solve exogenous state variables (see appendix A.6 in the paper). H is the last block of G.
            C=A+H*(nPhixxx);
        end
        %Block 1:xxx
        choosex3=kron(spx,kron(spx,spx));
        choosex3=logical(choosex3==1);
        tempeye=speye(n_x^3);
        Z=tempeye(:,choosex3);
        CZ=C*Z;
        if kamenik_type==1
            tic
            Xtemp=kamenik( G,D,hx_hat,-CZ,3 );
            time=toc;
        else
            tic
            [~,Xtemp]=gensylv( 3,G,D,hx_hat,full(-CZ));
            time=toc;
        end
        sol_time=sol_time+time;
        acc=norm(full(CZ+AkronkC(D*Xtemp,hx_hat,3)+G*Xtemp));
        if acc>1e-8
            warning(['Third order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        X=zeros(n_f,n_x^3);
        X(:,choosex3)=full(Xtemp); 
        gxxx_hat=Xtemp(1:n_y,:);
        hat_eta3_M3=eta3_M3(choosex3,:);
        %Block 2:xss
        choosex3=kron(kron(sps,sps),spx);
        choosex3=logical(choosex3==1);
        Z=tempeye(:,choosex3);
        CZ=C*Z;
        term2=reshape(fyp*gxxx_hat,n_f*(n_x-1),(n_x-1)^2)*hat_eta2_M2;
        term2=reshape(term2',n_f,n_x-1)*hx_hat;
        cons=CZ+term2;
        if kamenik_type==1
            tic
            Xtemp=kamenik( G,D,hx_hat,-cons,1 );
            time=toc;
        else
            tic
            [~,Xtemp]=gensylv( 1,G,D,hx_hat,full(-cons));
            time=toc;
        end
        sol_time=sol_time+time;
        acc=norm(full(cons+AkronkC(D*Xtemp,hx_hat,1)+G*Xtemp));
        if acc>1e-8
            warning(['Third order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        Xtemp=reshape(Xtemp,n_f,size(Z,2));
        X(:,choosex3)=full(Xtemp); clear Xtemp
        %Block 2:sss
        choosex3=kron(kron(sps,sps),sps);
        choosex3=logical(choosex3==1);
        Z=tempeye(:,choosex3);
        CZ=C*Z;
        term2=(fyp*gxxx_hat)*hat_eta3_M3;
        cons=CZ+term2;
        tic
        Xtemp=-(D+G)\cons;
        time=toc;
        sol_time=sol_time+time;

        acc=norm(full(cons+(D+G)*Xtemp));
        if acc>1e-8
            warning(['Third order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        Xtemp=reshape(Xtemp,n_f,size(Z,2));
        X(:,choosex3)=full(Xtemp); clear Xtemp
        % add symmetric entries, but first permute the indices, because the
        % compression matrices assume i1>=i2>=i3, while block xss has
        % i1<i2=i3.
        X=reshape(permute(reshape(X,[],n_x,n_x,n_x),[1,4,3,2]),[],n_x^3)*U3*W3; 
    else
        if n_x2==0
            CU3=A*U3;
        else
            CU3=A*U3+H*(nPhixxx*U3);
        end
        WBU_third_order; % create the matrix W3BU3
        %Block 1:xxx
        choosex3=kron(spx,kron(spx,spx));
        choosex3U=logical(U3'*choosex3~=0);
        tempeye=speye(size(U3,2));
        Z=tempeye(:,choosex3U);
        CU3Z=CU3*Z;
        ZTW3BU3Z=Z'*W3BU3*Z;
        if strcmp(algo,'vectorize')
            Xtemp=reshape(-(kron(ZTW3BU3Z',D)+kron(speye(size(Z,2)),G))\CU3Z(:),n_f,size(Z,2));
        elseif strcmp(algo,'dlyap')
            Xtemp=dlyap(ZTW3BU3Z',(-G\D)',(-G\CU3Z)');
            Xtemp=Xtemp';
        elseif strcmp(algo,'slicot')
            Xtemp=HessSchur(full(ZTW3BU3Z'),full((G\D)'),full((-G\CU3Z)'));
            Xtemp=Xtemp';
        end
        acc=norm(G*Xtemp+D*Xtemp*ZTW3BU3Z+CU3Z);
        if acc>1e-8
            warning(['Third order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        X=zeros(n_f,size(U3,2));
        X(:,choosex3U)=full(Xtemp); clear Xtemp ZTW3BU3Z
        %Block 2:xss
        choosex3=kron(kron(sps,sps),spx);
        choosex3=choosex3+kron(kron(sps,spx),sps);
        choosex3=choosex3+kron(kron(spx,sps),sps);
        choosex3U=logical(U3'*choosex3~=0);
        Z=tempeye(:,choosex3U);
        W3BU3Z=W3BU3*Z;
        CU3Z=CU3*Z+D*X*W3BU3Z;
        ZTW3BU3Z=Z'*W3BU3*Z;
        if strcmp(algo,'vectorize')
            Xtemp=reshape(-(kron(ZTW3BU3Z',D)+kron(speye(size(Z,2)),G))\CU3Z(:),n_f,size(Z,2));
        elseif strcmp(algo,'dlyap')
            Xtemp=dlyap(ZTW3BU3Z',(-G\D)',(-G\CU3Z)');
            Xtemp=Xtemp';
        elseif strcmp(algo,'slicot')
            Xtemp=HessSchur(full(ZTW3BU3Z'),full((G\D)'),full((-G\CU3Z)'));
            Xtemp=Xtemp';
        end
        acc=norm(G*Xtemp+D*Xtemp*ZTW3BU3Z+CU3Z);
        if acc>1e-8
            warning(['Third order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        X(:,choosex3U)=full(Xtemp); clear Xtemp ZTW3BU3Z
        %Block 3:sss
        choosex3=kron(kron(sps,sps),sps);
        choosex3U=logical(U3'*choosex3~=0);
        Z=tempeye(:,choosex3U);
        W3BU3Z=W3BU3*Z;
        CU3Z=CU3*Z+D*X*W3BU3Z;
        ZTW3BU3Z=Z'*W3BU3*Z;
        if strcmp(algo,'vectorize')
            Xtemp=reshape(-(kron(ZTW3BU3Z',D)+kron(speye(size(Z,2)),G))\CU3Z(:),n_f,size(Z,2));
        elseif strcmp(algo,'dlyap')
            Xtemp=dlyap(ZTW3BU3Z',(-G\D)',(-G\CU3Z)');
            Xtemp=Xtemp';
        elseif strcmp(algo,'slicot')
            Xtemp=HessSchur(full(ZTW3BU3Z'),full((G\D)'),full((-G\CU3Z)'));
            Xtemp=Xtemp';
        end
        acc=norm(G*Xtemp+D*Xtemp*ZTW3BU3Z+CU3Z);
        if acc>1e-8
            warning(['Third order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        X(:,choosex3U)=full(Xtemp); clear Xtemp ZTW3BU3Z W3BU3
        X=X*W3;
    end

    gxxx=X(1:n_y,:);
    hxxx=[X(n_y+1:end,:);nPhixxx;zeros(1,n_x^3)];
end

% Fourth Order
if approx>=4
    Vxxx0=[chain3(gx,gxx,gxxx,hx,hxx,hxxx,OMEGA_x.OMEGA1);gxxx;hxxx;sparse(n_x,n_x^3)];
    Vxxx1=[gxxx;sparse(n_y+2*n_x,n_x^3)];
    
    M4=M.M4(:);

    eta4_M4=reshape([eta*reshape(M4,n_e,n_e^3)]',n_e,n_e^2*n_x);
    eta4_M4=reshape([eta*eta4_M4]',n_e,n_e*n_x^2);
    eta4_M4=reshape([eta*eta4_M4]',n_e,n_x^3);
    eta4_M4=reshape([eta*eta4_M4]',n_x^4,1);

    Ezeta4=[ sparse(n_x^4,n_x^4-1) , eta4_M4 ];
    Ix_Ezeta3=kron(Ix,Ezeta3);
    Ix2=speye(n_x^2);
    Ix2_Ezeta2=kron(Ix2,Ezeta2);
    
    hx2=kron(hx,hx);
    hx_Ezeta2_hx=kron(kron(hx,Ezeta2),hx);
    hx_Ezeta3=kron(hx,Ezeta3);
    hx2_Ezeta2=kron(hx2,Ezeta2);
    
    start=tic;
    fvvvv=f_d4(nv',params,model.f_ind);
    deriv_time=deriv_time+toc(start);
    
    unique=nchoosek(n_x+4-1,4);
    unique=unique-nchoosek(n_x-1+3-1,3);

    clear result R
    A_fourth_order; % create matrix A of a fourth order solution

    if strcmp(algo,'Kamenik')
        if n_x2==0
            C=A;
        else % do not solve exogenous state variables (see appendix A.6 in the paper). H is the last block of G.
            C=A+H*(nPhixxxx);
        end
        %Block 1:xxxx
        choosex4=kron(kron(spx,spx),kron(spx,spx));
        choosex4=logical(choosex4==1);
        tempeye=speye(n_x^4);
        Z=tempeye(:,choosex4);
        CZ=C*Z;
        if kamenik_type==1
            tic
            Xtemp=kamenik( G,D,hx_hat,-CZ,4 );
            time=toc;
        else
            tic
            [~,Xtemp]=gensylv( 4,G,D,hx_hat,full(-CZ));
            time=toc;
        end
        sol_time=sol_time+time;
        acc=norm(full(CZ+AkronkC(D*Xtemp,hx_hat,4)+G*Xtemp));
        if acc>1e-8
            warning(['Fourth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        X=zeros(n_f,n_x^4);
        X(:,choosex4)=full(Xtemp); 
        gxxxx_hat=Xtemp(1:n_y,:);
        hat_eta4_M4=eta4_M4(choosex4,:);
        %Block 2:xxss
        choosex4=kron(kron(sps,sps),kron(spx,spx));
        choosex4=logical(choosex4==1);
        Z=tempeye(:,choosex4);
        CZ=C*Z;
        term2=reshape(fyp*gxxxx_hat,n_f*(n_x-1)^2,(n_x-1)^2)*hat_eta2_M2;
        term2=reshape(term2',n_f*(n_x-1),(n_x-1))*hx_hat;
        term2=reshape(term2',(n_x-1)*n_f,(n_x-1))*hx_hat;
        term2=reshape(term2',(n_x-1)^2,n_f);
        term2=term2';
        cons=CZ+term2;
        if kamenik_type==1
            tic
            Xtemp=kamenik( G,D,hx_hat,-cons,2 );
            time=toc;
        else
            tic
            [~,Xtemp]=gensylv( 2,G,D,hx_hat,full(-cons));
            time=toc;
        end
        sol_time=sol_time+time;

        acc=norm(full(cons+AkronkC(D*Xtemp,hx_hat,2)+G*Xtemp));
        if acc>1e-8
            warning(['Fourth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        Xtemp=reshape(Xtemp,n_f,size(Z,2));
        X(:,choosex4)=full(Xtemp);
        gxxss_hat=Xtemp(1:n_y,:);
        %Block 2:xsss
        choosex4=kron(kron(sps,sps),kron(sps,spx));
        choosex4=logical(choosex4==1);
        Z=tempeye(:,choosex4);
        CZ=C*Z;
        term2=reshape(fyp*gxxxx_hat,n_f*(n_x-1),(n_x-1)^3)*hat_eta3_M3;
        term2=reshape(term2',n_f,(n_x-1))*hx_hat;
        cons=CZ+term2;
        if kamenik_type==1
            tic
            Xtemp=kamenik( G,D,hx_hat,-cons,1 );
            time=toc;
        else
            tic
            [~,Xtemp]=gensylv( 1,G,D,hx_hat,full(-cons));
            time=toc;
        end
        sol_time=sol_time+time;
        acc=norm(full(cons+AkronkC(D*Xtemp,hx_hat,1)+G*Xtemp));
        if acc>1e-8
            warning(['Fourth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        Xtemp=reshape(Xtemp,n_f,size(Z,2));
        X(:,choosex4)=full(Xtemp); 
        %Block 2:ssss
        choosex4=kron(kron(sps,sps),kron(sps,sps));
        choosex4=logical(choosex4==1);
        Z=tempeye(:,choosex4);
        CZ=C*Z;
        term2=(6*fyp*gxxss_hat)*hat_eta2_M2;
        term3=(fyp*gxxxx_hat)*hat_eta4_M4;
        cons=CZ+term2+term3;
        tic
        Xtemp=-(D+G)\cons;
        time=toc;
        sol_time=sol_time+time;
        acc=norm(full(cons+(D+G)*Xtemp));
        if acc>1e-8
            warning(['Fourth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        Xtemp=reshape(Xtemp,n_f,size(Z,2));
        X(:,choosex4)=full(Xtemp); clear Xtemp
        % add symmetric entries, but first permute the indices.
        X=reshape(permute(reshape(X,[],n_x,n_x,n_x,n_x),[1,5,4,3,2]),[],n_x^4)*U4*W4;     
    else
        if n_x2==0
            CU4=A*U4;
        else
            CU4=A*U4+H*(nPhixxxx*U4);
        end
        kron_hx_hx=kron(hx,hx);
        WBU_fourth_order; % create the matrix W4BU4
        %Block 1:xxxx
        choosex4=kron(kron(spx,spx),kron(spx,spx));
        choosex4U=logical(U4'*choosex4~=0);
        tempeye=speye(size(U4,2));
        Z=tempeye(:,choosex4U);
        CU4Z=CU4*Z;
        ZTW4BU4Z=Z'*W4BU4*Z;
        if strcmp(algo,'vectorize')
            Xtemp=reshape(-(kron(ZTW4BU4Z',D)+kron(speye(size(Z,2)),G))\CU4Z(:),n_f,size(Z,2));
        elseif strcmp(algo,'dlyap')
            Xtemp=dlyap(ZTW4BU4Z',(-G\D)',(-G\CU4Z)');
            Xtemp=Xtemp';
        elseif strcmp(algo,'slicot')
            Xtemp=HessSchur(full(ZTW4BU4Z'),full((G\D)'),full((-G\CU4Z)'));
            Xtemp=Xtemp';
        end
        acc=norm(G*Xtemp+D*Xtemp*ZTW4BU4Z+CU4Z);
        if acc>1e-8
            warning(['Fourth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        X=zeros(n_f,size(U4,2));
        X(:,choosex4U)=full(Xtemp); clear Xtemp ZTW4BU4Z
        %Block 2:xxss
        choosex4=kron(kron(sps,sps),kron(spx,spx));
        choosex4=choosex4+kron(kron(sps,spx),kron(sps,spx));
        choosex4=choosex4+kron(kron(sps,spx),kron(spx,sps));
        choosex4=choosex4+kron(kron(spx,sps),kron(sps,spx));
        choosex4=choosex4+kron(kron(spx,sps),kron(spx,sps));
        choosex4=choosex4+kron(kron(spx,spx),kron(sps,sps));
        choosex4U=logical(U4'*choosex4~=0);
        Z=tempeye(:,choosex4U);
        W4BU4Z=W4BU4*Z;
        CU4Z=CU4*Z+D*X*W4BU4Z;
        ZTW4BU4Z=Z'*W4BU4*Z;
        if strcmp(algo,'vectorize')
            Xtemp=reshape(-(kron(ZTW4BU4Z',D)+kron(speye(size(Z,2)),G))\CU4Z(:),n_f,size(Z,2));
        elseif strcmp(algo,'dlyap')
            Xtemp=dlyap(ZTW4BU4Z',(-G\D)',(-G\CU4Z)');
            Xtemp=Xtemp';
        elseif strcmp(algo,'slicot')
            Xtemp=HessSchur(full(ZTW4BU4Z'),full((G\D)'),full((-G\CU4Z)'));
            Xtemp=Xtemp';
        end
        acc=norm(G*Xtemp+D*Xtemp*ZTW4BU4Z+CU4Z);
        if acc>1e-8
            warning(['Fourth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        X(:,choosex4U)=full(Xtemp); clear Xtemp ZTW4BU4Z
        %Block 3:xsss
        choosex4=kron(kron(sps,sps),kron(sps,spx));
        choosex4=choosex4+kron(kron(sps,sps),kron(spx,sps));
        choosex4=choosex4+kron(kron(sps,spx),kron(sps,sps));
        choosex4=choosex4+kron(kron(spx,sps),kron(sps,sps));
        choosex4U=logical(U4'*choosex4~=0);
        Z=tempeye(:,choosex4U);
        W4BU4Z=W4BU4*Z;
        CU4Z=CU4*Z+D*X*W4BU4Z;
        ZTW4BU4Z=Z'*W4BU4*Z;
        if strcmp(algo,'vectorize')
            Xtemp=reshape(-(kron(ZTW4BU4Z',D)+kron(speye(size(Z,2)),G))\CU4Z(:),n_f,size(Z,2));
        elseif strcmp(algo,'dlyap')
            Xtemp=dlyap(ZTW4BU4Z',(-G\D)',(-G\CU4Z)');
            Xtemp=Xtemp';
        elseif strcmp(algo,'slicot')
            Xtemp=HessSchur(full(ZTW4BU4Z'),full((G\D)'),full((-G\CU4Z)'));
            Xtemp=Xtemp';
        end
        acc=norm(G*Xtemp+D*Xtemp*ZTW4BU4Z+CU4Z);
        if acc>1e-8
            warning(['Fourth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        X(:,choosex4U)=full(Xtemp); clear Xtemp ZTW4BU4Z
        %Block 4:ssss
        choosex4=kron(kron(sps,sps),kron(sps,sps));
        choosex4U=logical(U4'*choosex4~=0);
        Z=tempeye(:,choosex4U);
        W4BU4Z=W4BU4*Z;
        CU4Z=CU4*Z+D*X*W4BU4Z;
        ZTW4BU4Z=Z'*W4BU4*Z;
        if strcmp(algo,'vectorize')
            Xtemp=reshape(-(kron(ZTW4BU4Z',D)+kron(speye(size(Z,2)),G))\CU4Z(:),n_f,size(Z,2));
        elseif strcmp(algo,'dlyap')
            Xtemp=dlyap(ZTW4BU4Z',(-G\D)',(-G\CU4Z)');
            Xtemp=Xtemp';
        elseif strcmp(algo,'slicot')
            Xtemp=HessSchur(full(ZTW4BU4Z'),full((G\D)'),full((-G\CU4Z)'));
            Xtemp=Xtemp';
        end
        acc=norm(G*Xtemp+D*Xtemp*ZTW4BU4Z+CU4Z);
        if acc>1e-8
            warning(['Fourth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        X(:,choosex4U)=full(Xtemp); clear Xtemp ZTW4BU4Z
        X=X*W4;
    end
    
    gxxxx=X(1:n_y,:);
    hxxxx=[X(n_y+1:end,:);nPhixxxx;zeros(1,n_x^4)];
end

% Fifth Order
if approx>=5
    Vxxxx0=[chain4(gx,gxx,gxxx,gxxxx,hx,hxx,hxxx,hxxxx,OMEGA_x.OMEGA2,OMEGA_x.OMEGA3,OMEGA_x.OMEGA4);gxxxx;hxxxx;sparse(n_x,n_x^4)];
    Vxxxx1=[gxxxx;sparse(n_y+2*n_x,n_x^4)];

    M5=M.M5(:);
    eta5_M5=reshape([eta*reshape(M5,n_e,n_e^4)]',n_e,n_e^3*n_x);
    eta5_M5=reshape([eta*eta5_M5]',n_e,n_e^2*n_x^2);
    eta5_M5=reshape([eta*eta5_M5]',n_e,n_e^1*n_x^3);
    eta5_M5=reshape([eta*eta5_M5]',n_e,n_x^4);
    eta5_M5=reshape([eta*eta5_M5]',n_x^5,1);

    start=tic;
    fvvvvv=f_d5(nv',params,model.f_ind);
    deriv_time=deriv_time+toc(start);

    clearvars -except n_f n_v n_x n_y n_x2 eta5_M5 Ix* *Ezeta* OMEGA_x U5 W5 gx* hx* nPhix* Vx* fv* fyp algo D G H approx sps spx kron_hx_hx hat_eta2_M2 hat_eta3_M3 hat_eta4_M4 sol_time deriv_time kamenik_type

    Ezeta5=[ sparse(n_x^5,n_x^5-1) , eta5_M5 ];
    Ix_Ezeta4=kron(Ix,Ezeta4);
    Ix2_Ezeta3=kron(Ix2,Ezeta3);
    Ix3=speye(n_x^3);
    Ix3_Ezeta2=kron(Ix3,Ezeta2);

    hx3=kron(hx2,hx);
    hx3_Ezeta2=kron(hx3,Ezeta2);
    hx2_Ezeta3=kron(hx2,Ezeta3);
    hx_Ezeta4=kron(hx,Ezeta4);

    unique=nchoosek(n_x+5-1,5);
    unique=unique-nchoosek(n_x-1+4-1,4);

    A_fifth_order; % create matrix A of a fifth order solution
    if strcmp(algo,'Kamenik')
        if n_x2==0
            C=A;
        else % do not solve exogenous state variables (see appendix A.6 in the paper). H is the last block of G.
            C=A+H*(nPhixxxxx);
        end
        %Block 1:xxxxx
        choosex5=kron(kron(spx,kron(spx,spx)),kron(spx,spx));
        choosex5=logical(choosex5==1);
        tempeye=speye(n_x^5);
        Z=tempeye(:,choosex5);
        CZ=C*Z;
        if kamenik_type==1
            tic
            Xtemp=kamenik( G,D,hx_hat,-CZ,5 );
            time=toc;
        else
            tic
            [~,Xtemp]=gensylv( 5,G,D,hx_hat,full(-CZ));
            time=toc;
        end
        sol_time=sol_time+time;
        acc=norm(full(CZ+AkronkC(D*Xtemp,hx_hat,5)+G*Xtemp));
        if acc>1e-8
            warning(['Fifth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        X=zeros(n_f,n_x^5);
        X(:,choosex5)=full(Xtemp); 
        gxxxxx_hat=Xtemp(1:n_y,:);
        hat_eta5_M5=eta5_M5(choosex5,:);
        %Block 2:xxxss
        choosex5=kron(kron(sps,sps),kron(spx,kron(spx,spx)));
        choosex5=logical(choosex5==1);
        Z=tempeye(:,choosex5);
        CZ=C*Z;
        term2=reshape(fyp*gxxxxx_hat,n_f*(n_x-1)^3,(n_x-1)^2)*hat_eta2_M2;
        term2=reshape(term2',n_f*(n_x-1)^2,(n_x-1))*hx_hat;
        term2=reshape(term2',(n_x-1)*n_f*(n_x-1),(n_x-1))*hx_hat;
        term2=reshape(term2',(n_x-1)^2*n_f,(n_x-1))*hx_hat;
        term2=reshape(term2',(n_x-1)^3,n_f);
        term2=term2';
        cons=CZ+term2;
        if kamenik_type==1
            tic
            Xtemp=kamenik( G,D,hx_hat,-cons,3 );
            time=toc;
        else
            tic
            [~,Xtemp]=gensylv( 3,G,D,hx_hat,full(-cons));
            time=toc;
        end
        sol_time=sol_time+time;
        acc=norm(full(cons+AkronkC(D*Xtemp,hx_hat,3)+G*Xtemp));
        if acc>1e-8
            warning(['Fifth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        Xtemp=reshape(Xtemp,n_f,size(Z,2));
        X(:,choosex5)=full(Xtemp);
        gxxxss_hat=Xtemp(1:n_y,:);
        %Block 3:xxsss
        choosex5=kron(kron(sps,sps),kron(sps,kron(spx,spx)));
        choosex5=logical(choosex5==1);
        Z=tempeye(:,choosex5);
        CZ=C*Z;
        term2=reshape(fyp*gxxxxx_hat,n_f*(n_x-1)^2,(n_x-1)^3)*hat_eta3_M3;
        term2=reshape(term2',n_f*(n_x-1),(n_x-1))*hx_hat;
        term2=reshape(term2',(n_x-1)*n_f,(n_x-1))*hx_hat;
        term2=reshape(term2',(n_x-1)^2,n_f);
        term2=term2';
        cons=CZ+term2;
        if kamenik_type==1
            tic
            Xtemp=kamenik( G,D,hx_hat,-cons,2 );
            time=toc;
        else
            tic
            [~,Xtemp]=gensylv( 2,G,D,hx_hat,full(-cons));
            time=toc;
        end
        sol_time=sol_time+time;
        acc=norm(full(cons+AkronkC(D*Xtemp,hx_hat,2)+G*Xtemp));
        if acc>1e-8
            warning(['Fifth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        Xtemp=reshape(Xtemp,n_f,size(Z,2));
        X(:,choosex5)=full(Xtemp); 
        gxxsss_hat=Xtemp(1:n_y,:);
        %Block 2:xssss
        choosex5=kron(kron(sps,sps),kron(kron(sps,sps),spx));
        choosex5=logical(choosex5==1);
        Z=tempeye(:,choosex5);
        CZ=C*Z;
        term2=reshape(6*fyp*gxxxss_hat,n_f*(n_x-1),(n_x-1)^2)*hat_eta2_M2;
        term2=reshape(term2',n_f,(n_x-1))*hx_hat;
        term3=reshape(fyp*gxxxxx_hat,n_f*(n_x-1),(n_x-1)^4)*hat_eta4_M4;
        term3=reshape(term3',n_f,(n_x-1))*hx_hat;
        cons=CZ+term2+term3;
        if kamenik_type==1
            tic
            Xtemp=kamenik( G,D,hx_hat,-cons,1 );
            time=toc;
        else
            tic
            [~,Xtemp]=gensylv( 1,G,D,hx_hat,full(-cons));
            time=toc;
        end
        sol_time=sol_time+time;
        acc=norm(full(cons+AkronkC(D*Xtemp,hx_hat,1)+G*Xtemp));
        if acc>1e-8
            warning(['Fifth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        Xtemp=reshape(Xtemp,n_f,size(Z,2));
        X(:,choosex5)=full(Xtemp);         
        %Block 5:sssss
        choosex5=kron(kron(sps,sps),kron(kron(sps,sps),sps));
        choosex5=logical(choosex5==1);
        Z=tempeye(:,choosex5);
        CZ=C*Z;
        term2=(10*fyp*gxxsss_hat)*hat_eta2_M2;
        term3=(10*fyp*gxxxss_hat)*hat_eta3_M3;
        term4=(fyp*gxxxxx_hat)*hat_eta5_M5;
        cons=CZ+term2+term3+term4;
        tic
        Xtemp=-(D+G)\cons;
        time=toc;
        sol_time=sol_time+time;
        acc=norm(full(cons+(D+G)*Xtemp));
        if acc>1e-8
            warning(['Fifth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        Xtemp=reshape(Xtemp,n_f,size(Z,2));
        X(:,choosex5)=full(Xtemp); clear Xtemp
        % add symmetric entries, but first permute the indices.
        X=reshape(permute(reshape(X,[],n_x,n_x,n_x,n_x,n_x),[1,6,5,4,3,2]),[],n_x^5)*U5*W5;     
    else
        if n_x2==0
            CU5=A*U5;
        else
            CU5=A*U5+H*(nPhixxxxx*U5);
        end
        WBU_fifth_order; % create the matrix W5BU5
        clearvars -except G D W5BU5 CU5 W5 gx* hx hxx* algo n_y nPhixxxxx n_x n_f unique approx spx sps U5
        %Block 1:xxxxx
        choosex5=kron(kron(kron(spx,spx),kron(spx,spx)),spx);
        choosex5U=logical(U5'*choosex5~=0);
        tempeye=speye(size(U5,2));
        Z=tempeye(:,choosex5U);
        CU5Z=CU5*Z;
        ZTW5BU5Z=Z'*W5BU5*Z;
        if strcmp(algo,'vectorize')
            Xtemp=reshape(-(kron(ZTW5BU5Z',D)+kron(speye(size(Z,2)),G))\CU5Z(:),n_f,size(Z,2));
        elseif strcmp(algo,'dlyap')
            Xtemp=dlyap(ZTW5BU5Z',(-G\D)',(-G\CU5Z)');
            Xtemp=Xtemp';
        elseif strcmp(algo,'slicot')
            Xtemp=HessSchur(full(ZTW5BU5Z'),full((G\D)'),full((-G\CU5Z)'));
            Xtemp=Xtemp';
        end
        acc=norm(G*Xtemp+D*Xtemp*ZTW5BU5Z+CU5Z);
        if acc>1e-8
            warning(['Fifth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        X=zeros(n_f,size(U5,2));
        X(:,choosex5U)=full(Xtemp); clear Xtemp ZTW5BU5Z
        %Block 2:xxxss
        choosex5=kron(kron(kron(sps,sps),kron(spx,spx)),spx);
        choosex5=choosex5+kron(kron(kron(sps,spx),kron(sps,spx)),spx);
        choosex5=choosex5+kron(kron(kron(sps,spx),kron(spx,sps)),spx);
        choosex5=choosex5+kron(kron(kron(sps,spx),kron(spx,spx)),sps);
        choosex5=choosex5+kron(kron(kron(spx,sps),kron(sps,spx)),spx);
        choosex5=choosex5+kron(kron(kron(spx,sps),kron(spx,sps)),spx);
        choosex5=choosex5+kron(kron(kron(spx,sps),kron(spx,spx)),sps);
        choosex5=choosex5+kron(kron(kron(spx,spx),kron(sps,sps)),spx);
        choosex5=choosex5+kron(kron(kron(spx,spx),kron(sps,spx)),sps);
        choosex5=choosex5+kron(kron(kron(spx,spx),kron(spx,sps)),sps);
        choosex5U=logical(U5'*choosex5~=0);
        Z=tempeye(:,choosex5U);
        W5BU5Z=W5BU5*Z;
        CU5Z=CU5*Z+D*X*W5BU5Z;
        ZTW5BU5Z=Z'*W5BU5*Z;
        if strcmp(algo,'vectorize')
            Xtemp=reshape(-(kron(ZTW5BU5Z',D)+kron(speye(size(Z,2)),G))\CU5Z(:),n_f,size(Z,2));
        elseif strcmp(algo,'dlyap')
            Xtemp=dlyap(ZTW5BU5Z',(-G\D)',(-G\CU5Z)');
            Xtemp=Xtemp';
        elseif strcmp(algo,'slicot')
            Xtemp=HessSchur(full(ZTW5BU5Z'),full((G\D)'),full((-G\CU5Z)'));
            Xtemp=Xtemp';
        end
        acc=norm(G*Xtemp+D*Xtemp*ZTW5BU5Z+CU5Z);
        if acc>1e-8
            warning(['Fifth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        X(:,choosex5U)=full(Xtemp); clear Xtemp ZTW5BU5Z
        %Block 3:xxsss
        choosex5=kron(kron(kron(spx,spx),kron(sps,sps)),sps);
        choosex5=choosex5+kron(kron(kron(spx,sps),kron(spx,sps)),sps);
        choosex5=choosex5+kron(kron(kron(spx,sps),kron(sps,spx)),sps);
        choosex5=choosex5+kron(kron(kron(spx,sps),kron(sps,sps)),spx);
        choosex5=choosex5+kron(kron(kron(sps,spx),kron(spx,sps)),sps);
        choosex5=choosex5+kron(kron(kron(sps,spx),kron(sps,spx)),sps);
        choosex5=choosex5+kron(kron(kron(sps,spx),kron(sps,sps)),spx);
        choosex5=choosex5+kron(kron(kron(sps,sps),kron(spx,spx)),sps);
        choosex5=choosex5+kron(kron(kron(sps,sps),kron(spx,sps)),spx);
        choosex5=choosex5+kron(kron(kron(sps,sps),kron(sps,spx)),spx);
        choosex5U=logical(U5'*choosex5~=0);
        Z=tempeye(:,choosex5U);
        W5BU5Z=W5BU5*Z;
        CU5Z=CU5*Z+D*X*W5BU5Z;
        ZTW5BU5Z=Z'*W5BU5*Z;
        if strcmp(algo,'vectorize')
            Xtemp=reshape(-(kron(ZTW5BU5Z',D)+kron(speye(size(Z,2)),G))\CU5Z(:),n_f,size(Z,2));
        elseif strcmp(algo,'dlyap')
            Xtemp=dlyap(ZTW5BU5Z',(-G\D)',(-G\CU5Z)');
            Xtemp=Xtemp';
        elseif strcmp(algo,'slicot')
            Xtemp=HessSchur(full(ZTW5BU5Z'),full((G\D)'),full((-G\CU5Z)'));
            Xtemp=Xtemp';
        end
        acc=norm(G*Xtemp+D*Xtemp*ZTW5BU5Z+CU5Z);
        if acc>1e-8
            warning(['Fifth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        X(:,choosex5U)=full(Xtemp); clear Xtemp ZTW5BU5Z
        %Block 4:xssss
        choosex5=kron(kron(kron(spx,sps),kron(sps,sps)),sps);
        choosex5=choosex5+kron(kron(kron(sps,spx),kron(sps,sps)),sps);
        choosex5=choosex5+kron(kron(kron(sps,sps),kron(spx,sps)),sps);
        choosex5=choosex5+kron(kron(kron(sps,sps),kron(sps,spx)),sps);
        choosex5=choosex5+kron(kron(kron(sps,sps),kron(sps,sps)),spx);
        choosex5U=logical(U5'*choosex5~=0);
        Z=tempeye(:,choosex5U);
        W5BU5Z=W5BU5*Z;
        CU5Z=CU5*Z+D*X*W5BU5Z;
        ZTW5BU5Z=Z'*W5BU5*Z;
        if strcmp(algo,'vectorize')
            Xtemp=reshape(-(kron(ZTW5BU5Z',D)+kron(speye(size(Z,2)),G))\CU5Z(:),n_f,size(Z,2));
        elseif strcmp(algo,'dlyap')
            Xtemp=dlyap(ZTW5BU5Z',(-G\D)',(-G\CU5Z)');
            Xtemp=Xtemp';
        elseif strcmp(algo,'slicot')
            Xtemp=HessSchur(full(ZTW5BU5Z'),full((G\D)'),full((-G\CU5Z)'));
            Xtemp=Xtemp';
        end
        acc=norm(G*Xtemp+D*Xtemp*ZTW5BU5Z+CU5Z);
        if acc>1e-8
            warning(['Fifth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        X(:,choosex5U)=full(Xtemp); clear Xtemp ZTW5BU5Z
        %Block 5:sssss
        choosex5=kron(kron(kron(sps,sps),kron(sps,sps)),sps);
        choosex5U=logical(U5'*choosex5~=0);
        Z=tempeye(:,choosex5U);
        W5BU5Z=W5BU5*Z;
        CU5Z=CU5*Z+D*X*W5BU5Z;
        ZTW5BU5Z=Z'*W5BU5*Z;
        if strcmp(algo,'vectorize')
            Xtemp=reshape(-(kron(ZTW5BU5Z',D)+kron(speye(size(Z,2)),G))\CU5Z(:),n_f,size(Z,2));
        elseif strcmp(algo,'dlyap')
            Xtemp=dlyap(ZTW5BU5Z',(-G\D)',(-G\CU5Z)');
            Xtemp=Xtemp';
        elseif strcmp(algo,'slicot')
            Xtemp=HessSchur(full(ZTW5BU5Z'),full((G\D)'),full((-G\CU5Z)'));
            Xtemp=Xtemp';
        end
        acc=norm(G*Xtemp+D*Xtemp*ZTW5BU5Z+CU5Z);
        if acc>1e-8
            warning(['Fifth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
        end
        X(:,choosex5U)=full(Xtemp); clear Xtemp ZTW5BU5Z W5BU5
        X=X*W5;
    end
    
    gxxxxx=X(1:n_y,:);
    hxxxxx=[X(n_y+1:end,:);nPhixxxxx;zeros(1,n_x^5)];
end

clear derivs
derivs.gx=full(gx);
derivs.hx=full(hx(1:end-1,:));
if approx>=2
    derivs.gxx=full(gxx);
    derivs.hxx=full(hxx(1:end-1,:));
end
if approx>=3
    derivs.gxxx=full(gxxx);
    derivs.hxxx=full(hxxx(1:end-1,:));
end
if approx>=4
    derivs.gxxxx=full(gxxxx);
    derivs.hxxxx=full(hxxxx(1:end-1,:));
end
if approx>=5
    derivs.gxxxxx=full(gxxxxx);
    derivs.hxxxxx=full(hxxxxx(1:end-1,:));
end

end



