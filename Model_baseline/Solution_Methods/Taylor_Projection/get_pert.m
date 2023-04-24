function [derivs,stoch_pert,nonstoch_pert,model]=get_pert(model,params,M,eta,nxss,nyss,algo)
%get_pert(model,params,M,eta,nxss,nyss,algo)
%This function solves the model with perturbation up to fourth
%order.
%Input variables:
%model: a structure that is generated automatically by the function
%prepare_model.m
%params: a vector of all parameter values ordered in the same order of symparams.
%M: a structure that contains all the cross moments of the shocks. The fields of
%this structure should be M2,M3,M4.
%eta: The matrix eta as defined in Schmitt-Grohe and Uribe (2004).
%nxss,nyss: steady state values of x and y
%algo: algorithm type. 'dlyap' for dlyap or 'vectorize' for vectorization.
%
% © Copyright, Oren Levintal, June 13, 2016.


approx=model.order(2);

if approx>4
    error('perturbation order cannot exceed 4')
end

tic
n_f=model.n_f;

n_y=model.n_y;
n_x2=model.n_x2;
n_x1=model.n_x1;
n_z=model.n_z; % no. of auxiliary variables

% note that pert system has one more state variable than tp system.
n_x_tp=model.n_x;
n_v_tp=model.n_v;
n_x=model.n_x+1;
n_v=model.n_v+2;

if isfield(model,'UW')
    UW=model.UW; % IMPORTANT: THESE COMPRESSION MATRICES ASSUME A CERTAIN NONZERO PATTERN.
else
    UW=create_compression_matrices_nonzero2(approx,n_x); % compression matrices for the perturbation solution.
    model.UW=UW;
end

if approx==2
    W2=UW.W2; U2=UW.U2;
elseif approx==3
    W2=UW.W2; U2=UW.U2;
    W3=UW.W3; U3=UW.U3;
elseif approx==4
    W2=UW.W2; U2=UW.U2;
    W3=UW.W3; U3=UW.U3;
    W4=UW.W4; U4=UW.U4;
end

if isfield(model,'OMEGA_x')
    OMEGA_x=model.OMEGA_x;
else
    create_OMEGA_x;
    model.OMEGA_x=OMEGA_x;
end


n_e=size(eta,2);
eta=[eta;zeros(1,n_e)];

nx=nxss;
% nxss=[nxss(:);0]; % add steady state value of the perturbation variable.
nxss=nxss(:);
nyss=nyss(:);
% nuss=eval_u([nyss;nyss;nxss;nxss],params);
if ~isfield(model,'n_ind')
    n_ind=1;
else
    n_ind=model.n_ind;
end
if ~isfield(model,'maxload')
    maxload=intarray(1);
else
    maxload=intarray(model.maxload);
end
if ~isfield(model,'pertind')
    model.pertind=cell(model.totindi,1);
end
pertindi=1;

% nPhi=Phi_fun(nx,params);
get_Phi_derivs;
   

get_f_derivs;

% numu_derivs=evalu_derivs(params,nxss,nyss,approx);
% call_evalf_derivs;
% 
% sym2script_fv;




%%%%%%%%%%%%%
fyp=fv(:,1:n_y); fy=fv(:,n_y+1:2*n_y); fxp=fv(:,2*n_y+1:2*n_y+n_x); fx=fv(:,2*n_y+n_x+1:end);

% First Order

% use the function gx_hx of Schmitt-Grohe and Uribe (2004) to solve the
% first order solution

% tic
[gx,hx,exitflag]=gx_hx(full([fy;zeros(n_x2,n_y)]),full([fx(:,1:end-1);nPhix(:,1:end-1)]),full([fyp;zeros(n_x2,n_y)]),full([fxp(:,1:end-1);[zeros(n_x2,n_x1),-eye(n_x2)]]));

gx=[gx,zeros(n_y,1)];
hx=[hx;zeros(1,n_x-1)];
hx=[hx,zeros(n_x,1)];
hx(end,end)=1;
% replace relevant rows of hx with nPhix which is more numerically accurate
hx(n_x1+1:n_x1+n_x2,:)=nPhix;

% disp(['1st order completed in ' num2str(toc,'%15.2f') ' seconds'])
    
hx=sparse(hx);
gx=sparse(gx);

% Second Order
if approx>=2
%     tic
    M2=M.M2;
    eta2_M2=reshape([eta*reshape(M2,n_e,n_e)]',n_e,n_x);
    eta2_M2=reshape([eta*eta2_M2]',n_x^2,1);

%     sym2script_fvv;

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
    if strcmp(algo,'vectorize')
        Xtemp=reshape(-(kron(ZTW2BU2Z',D)+kron(speye(size(Z,2)),G))\CU2Z(:),n_f,size(Z,2));
    elseif strcmp(algo,'dlyap')
        Xtemp=dlyap(ZTW2BU2Z',(-G\D)',(-G\CU2Z)');
        Xtemp=Xtemp';
    end
    acc=norm(full(G*Xtemp+D*Xtemp*ZTW2BU2Z+CU2Z));
    if acc>1e-10
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
    if strcmp(algo,'vectorize')
        Xtemp=reshape(-(kron(ZTW2BU2Z',D)+kron(speye(size(Z,2)),G))\CU2Z(:),n_f,size(Z,2));
    elseif strcmp(algo,'dlyap')
        Xtemp=dlyap(ZTW2BU2Z',(-G\D)',(-G\CU2Z)');
        Xtemp=Xtemp';
    end
    acc=norm(G*Xtemp+D*Xtemp*ZTW2BU2Z+CU2Z);
    if acc>1e-10
        warning(['Second order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
    end
    X(:,choosex2U)=full(Xtemp); clear Xtemp W2BU2Z ZTW2BU2Z
    X=X*W2;

    gxx=X(1:n_y,:);
    hxx=[X(n_y+1:end,:);nPhixx;zeros(1,n_x^2)];
%     disp(['2nd order completed in ' num2str(toc,'%15.2f') ' seconds'])
end

% Third Order
if approx>=3
%     tic
  
    Vxx0=[chain2(gx,gxx,hx,hxx);gxx;hxx;sparse(n_x,n_x^2)];
    Vxx1=[gxx;sparse(n_y+2*n_x,n_x^2)];
    
    M3=M.M3(:);

    eta3_M3=reshape([eta*reshape(M3,n_e,n_e^2)]',n_e,n_e*n_x);
    eta3_M3=reshape([eta*eta3_M3]',n_e,n_x^2);
    eta3_M3=reshape([eta*eta3_M3]',n_x^3,1);

    Ezeta3=[ sparse(n_x^3,n_x^3-1) , eta3_M3 ];
    Ix=speye(n_x);
    Ix_Ezeta2=kron(Ix,Ezeta2);
    
    
%     sym2script_fvvv;

    unique=nchoosek(n_x+3-1,3);
    unique=unique-nchoosek(n_x-1+2-1,2);
%     [Ui,~]=find(U3);
    A_third_order; % create matrix A of a third order solution

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
    end
    acc=norm(G*Xtemp+D*Xtemp*ZTW3BU3Z+CU3Z);
    if acc>1e-10
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
    end
    acc=norm(G*Xtemp+D*Xtemp*ZTW3BU3Z+CU3Z);
    if acc>1e-10
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
    if acc>1e-10
        warning(['Third order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
    end
    X(:,choosex3U)=full(Xtemp); clear Xtemp ZTW3BU3Z W3BU3
    X=X*W3;

    gxxx=X(1:n_y,:);
    hxxx=[X(n_y+1:end,:);nPhixxx;zeros(1,n_x^3)];
%     disp(['3rd order completed in ' num2str(toc,'%15.2f') ' seconds'])
end

% Fourth Order
if approx>=4
%     tic

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
    
%     sym2script_fvvvv;

    unique=nchoosek(n_x+4-1,4);
    unique=unique-nchoosek(n_x-1+3-1,3);

    clear result R

    A_fourth_order; % create matrix A of a fourth order solution

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
    end
    acc=norm(G*Xtemp+D*Xtemp*ZTW4BU4Z+CU4Z);
    if acc>1e-10
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
    end
    acc=norm(G*Xtemp+D*Xtemp*ZTW4BU4Z+CU4Z);
    if acc>1e-10
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
    end
    acc=norm(G*Xtemp+D*Xtemp*ZTW4BU4Z+CU4Z);
    if acc>1e-10
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
    end
    acc=norm(G*Xtemp+D*Xtemp*ZTW4BU4Z+CU4Z);
    if acc>1e-10
        warning(['Fourth order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
    end
    X(:,choosex4U)=full(Xtemp); clear Xtemp ZTW4BU4Z
    X=X*W4;
    
    gxxxx=X(1:n_y,:);
    hxxxx=[X(n_y+1:end,:);nPhixxxx;zeros(1,n_x^4)];
%     disp(['4th order completed in ' num2str(toc,'%15.2f') ' seconds'])
end



clear derivs
derivs.gx=full(gx);
derivs.hx=full(hx(1:end-1,:));
derivs.gxx=sparse(n_y,n_x^2);
derivs.hxx=sparse(n_x-1,n_x^2);
derivs.gxxx=sparse(n_y,n_x^3);
derivs.hxxx=sparse(n_x-1,n_x^3);

if approx>=2
    derivs.gxx=reshape(full(gxx),n_y,n_x,n_x);
    derivs.hxx=reshape(full(hxx(1:end-1,:)),n_x-1,n_x,n_x);
end
if approx>=3
    derivs.gxxx=reshape(full(gxxx),n_y,n_x,n_x,n_x);
    derivs.hxxx=reshape(full(hxxx(1:end-1,:)),n_x-1,n_x,n_x,n_x);
end
if approx>=4
    derivs.gxxxx=reshape(full(gxxxx),n_y,n_x,n_x,n_x,n_x);
    derivs.hxxxx=reshape(full(hxxxx(1:end-1,:)),n_x-1,n_x,n_x,n_x,n_x);
end

[ stoch_pert, nonstoch_pert ] = get_initial( model,model.order(1),derivs,nyss,nxss );





