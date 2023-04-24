function [T,J,model]=tp_no_shift(coeffs,x0,model,params,eta,c0,nep,P,varargin)

% This function is like tp, but the center of the power series can be
% different from x0.
%
% © Copyright, Oren Levintal, June 13, 2016.

order=model.order(1);

% make sure all inputs are full
coeffs=full(coeffs);
x0=full(x0);
params=full(params);
eta=full(eta);
c0=full(c0);
nep=full(nep);
P=full(P);

c0_old=c0;
c0_new=x0;

% if c0~=x0 % shift the center of the power series to x0
%     warning('shifting center, Jacobian may be incorrect')
%     [ coeffs ] = shift_center( coeffs,c0_old,c0_new,order,model );
%     coeffs=full(coeffs);
%     c0=c0_new;
% end

if isempty(varargin)
   [ g_theta,h_theta,gx_theta,hx_theta,gxxc_theta,hxxc_theta,gxxxc_theta,hxxxc_theta,model ] = precompute_no_shift(x0,c0,model,order ); 
else
    g_theta=varargin{1};
    h_theta=varargin{2};
    gx_theta=varargin{3};
    hx_theta=varargin{4};
    gxxc_theta=varargin{5};
    hxxc_theta=varargin{6};
    gxxxc_theta=varargin{7};
    hxxxc_theta=varargin{8};
end

indi=5; %count indi from 5, due to precompute.

if ~isfield(model,'jacobian')
    model.jacobian='exact';
end
if ~isfield(model,'n_ind')
    n_ind=1;
else
    n_ind=model.n_ind;
end
if ~isfield(model,'maxload')
    maxload=intarray(length(P));
else
    maxload=intarray(model.maxload);
end

params=params(:);
P=P(:);

n_f=model.n_f; % no. of model conditions
n_x=model.n_x; % no. of state variables
n_y=model.n_y; % no. of control variables
n_x2=model.n_x2; % no. of exogenous state variables
n_x1=model.n_x1; % no. of endogenous state variabels
n_v=model.n_v; % total no. of variables v=(yp,y,xp,x)
n_theta=model.n_theta; % no. of Polynomial coefficients
n_b=model.n_b; % no. of terms in the basis function
n_nodes=size(nep,2); % no. nodes in the discrete shocks
n_z=model.n_z; % no. of auxiliary variables

coeffs=reshape(coeffs,n_f,n_b);
g_coeffs=coeffs(1:n_y,:);

%create stochy: index of future control variables that affect stochastic equations
tempv=zeros(n_v,1);
tempv(model.stochfvars)=model.stochfvars;
tempv(n_y+1:end)=0;
stochy=nonzeros(tempv);
n_stochy=length(stochy);
stochg_coeffs=coeffs(stochy,:); 

% U=model.U;

unique2=[]; unique3=[];

prefvars=model.prefvars; % index of variables that affect nonstochastic equations.
stochfvars=model.stochfvars;% index of variables that affect stochastic equations.
n_prefvars=length(prefvars);
n_stochfvars=length(stochfvars);

% lower block of v_theta, vx_theta and so on

v_theta_lower_terms=vconcat(vconcat(g_theta,h_theta),sptensor(n_x,n_theta));

if order>=1
    vx_theta_lower_terms=vconcat(vconcat(gx_theta,hx_theta),sptensor(n_x,[n_x,n_theta]));
end
if order>=2
    unique2=model.unique2;
    vxxc_theta_lower_terms=vconcat(vconcat(gxxc_theta,hxxc_theta),sptensor(n_x,[unique2,n_theta]));
end
if order>=3
    unique3=model.unique3;
    vxxxc_theta_lower_terms=vconcat(vconcat(gxxxc_theta,hxxxc_theta),sptensor(n_x,[unique3,n_theta]));
end

% current state
nx=x0;

% evaluate policy functions
[X,Xx,Xxx,Xxxx,Xxxxx]=create_X(order,n_x,x0,c0,model.W{2},model.unique2,model.W{3},model.unique3);
g=coeffs(1:n_y,:)*X;
h1=coeffs(n_y+1:end,:)*X;
h=zeros(n_x,1);

% derivatives of policy functions at x0 (assuming c0=x0)
if order>=1
    gx=coeffs(1:n_y,:)*Xx;
    h1x=coeffs(n_y+1:end,:)*Xx;
end

if order>=2
    gxxc=coeffs(1:n_y,:)*Xxx*model.U{2};
    h1xxc=coeffs(n_y+1:end,:)*Xxx*model.U{2};
end
if order>=3
    gxxxc=coeffs(1:n_y,:)*Xxxx*model.U{3};
    h1xxxc=coeffs(n_y+1:end,:)*Xxxx*model.U{3};
end

% control vars
ny=g;

% STEP 1: COMPUTE NONSTOCHASTIC EQUATIONS AND DERIVATIVES, AND THEIR JACOBIAN

% build h(x)
h(1:n_x1,1)=h1; % predetermined endogenous state vars
h(n_x1+1:end,1)=Phi_fun(nx,params); % expected value of exogenous state vars. shocks are added later

% expected value of next period state variables.
nxp=h;

% derivatives of h at x0
if order>=1
    hx=vconcat(sptensor(h1x),Phi_d_c1(nx(:)',params,model.Phi_indc));
end
if order>=2
    hxxc=vconcat(sptensor(h1xxc),Phi_d_c2(nx(:)',params,model.Phi_indc));
end
if order>=3
    hxxxc=vconcat(sptensor(h1xxxc),Phi_d_c3(nx(:)',params,model.Phi_indc));
end

% evaluate residuals R0
nv=[zeros(n_y,1);ny;nxp;nx]; % model variables, yp are zero for now, and nxp is an expected value. stochastic components will be added later.

n_u=model.n_u;
nu=zeros(n_u,1);
npreu=preu_fun(nv(model.preuvars),params); % all predetermined u
nu(model.preurows)=npreu;
nz=[nv;nu]; % auxiliary variables
pref=pretilf_fun(nz(model.pretilfzvars),params); % since i already have z, i use pretilf to evaluate f.


% compute derivatives of prePI and pretilf w.r.t z (prePI=nonstochastic rows of PI, pretilf=nonstochastic rows of f)

prePIz_full=prePI_d1(nz',params,model.prePI_ind_u);
pretilfz=pretilf_tilf_d1([nv(prefvars);nu(model.prefuvars)]',params,model.pretilf_ind_u);
if order>=1
    prePIzz_full=prePI_d2(nz',params,model.prePI_ind_u);
    pretilfzz=pretilf_tilf_d2([nv(prefvars);nu(model.prefuvars)]',params,model.pretilf_ind_u);
end
if order>=2
    prePIzzz_full=prePI_d3(nz',params,model.prePI_ind_u);
    pretilfzzz=pretilf_tilf_d3([nv(prefvars);nu(model.prefuvars)]',params,model.pretilf_ind_u);
end
if order>=3
    prePIzzzz_full=prePI_d4(nz',params,model.prePI_ind_u);
    pretilfzzzz=pretilf_tilf_d4([nv(prefvars);nu(model.prefuvars)]',params,model.pretilf_ind_u);
end

% extract derivatives of prePI w.r.t v from prePIz, prePIzz, ...
if order>=0
    [prePIz,model.ind{indi}]=extract(prePIz_full,model.prefzvars,model.prefzvars,0,model.ind{indi});
    indi=indi+1;
    [prePIv,model.ind{indi}]=extract(prePIz_full,model.prefzvars,model.prefvars,0,model.ind{indi});
    indi=indi+1;
end
if order>=1
    [prePIzz,model.ind{indi}]=extract(prePIzz_full,model.prefzvars,model.prefzvars,0,model.ind{indi});
    indi=indi+1;
    [prePIvvc,model.ind{indi}]=extract(prePIzz_full,model.prefzvars,model.prefvars,1,model.ind{indi});
    indi=indi+1;
end
if order>=2
    [prePIzzz,model.ind{indi}]=extract(prePIzzz_full,model.prefzvars,model.prefzvars,0,model.ind{indi});
    indi=indi+1;
    [prePIvvvc,model.ind{indi}]=extract(prePIzzz_full,model.prefzvars,model.prefvars,1,model.ind{indi});
    indi=indi+1;
    if ~isfield(model,'prefvars_chain3c_M2')
        [ tempM ] = chainsM( n_prefvars,3 );
        model.prefvars_chain3c_M2=tempM{2};
        clear tempM
    end
end
if order>=3
    [prePIzzzz,model.ind{indi}]=extract(prePIzzzz_full,model.prefzvars,model.prefzvars,0,model.ind{indi});
    indi=indi+1;
    [prePIvvvvc,model.ind{indi}]=extract(prePIzzzz_full,model.prefzvars,model.prefvars,1,model.ind{indi});
    indi=indi+1;
    if ~isfield(model,'prefvars_chain4c_M2')
        [ tempM ] = chainsM( n_prefvars,4 );
        model.prefvars_chain4c_M2=tempM{2};
        model.prefvars_chain4c_M3=tempM{3};
        model.prefvars_chain4c_M4=tempM{4};
        clear tempM
    end
end


% Use high order chain rules to compute derivatives of prePI w.r.t v
totindi=5+8;
if order==0
    indi=totindi;
    for i=2:model.pre_n
        [prePIv,model.ind{indi}]=chain1_tensor(prePIz,prePIv,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    end
    [prefv,model.ind{indi}]=chain1_tensor(pretilfz,prePIv,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
elseif order==1
    indi=totindi+(model.pre_n-1)*(1)+1;
    for i=2:model.pre_n
        [prePIv,model.ind{indi}]=colsort(prePIv,model.ind{indi});
        indi=indi+1;
        [prePIvvc,model.ind{indi}]=chain2c_tensor(prePIz,prePIzz,prePIv,prePIvvc,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [prePIv,model.ind{indi}]=chain1_tensor(prePIz,prePIv,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    end
    [prePIv,model.ind{indi}]=colsort(prePIv,model.ind{indi});
    indi=indi+1;

    [prefv,model.ind{indi}]=chain1_tensor(pretilfz,prePIv,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [prefvvc,model.ind{indi}]=chain2c_tensor(pretilfz,pretilfzz,prePIv,prePIvvc,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;

    [ prefvv,model.ind{indi} ] = uncompressderivs( prefvvc,2,n_prefvars,model.fv(model.prefrows,prefvars),model.ind{indi} );
    indi=indi+1;

elseif order==2
    indi=totindi+(model.pre_n-1)*(1+3)+1+4;
    for i=2:model.pre_n
        [prePIv,model.ind{indi}]=colsort(prePIv,model.ind{indi});
        indi=indi+1;
        [prePIvvvc,model.ind{indi}]=chain3c_tensor(prePIz,prePIzz,prePIzzz,...
            prePIv,prePIvvc,prePIvvvc,model.prefvars_chain3c_M2,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [prePIvvc,model.ind{indi}]=chain2c_tensor(prePIz,prePIzz,prePIv,prePIvvc,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [prePIv,model.ind{indi}]=chain1_tensor(prePIz,prePIv,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    end
    [prePIv,model.ind{indi}]=colsort(prePIv,model.ind{indi});
    indi=indi+1;

    [prefv,model.ind{indi}]=chain1_tensor(pretilfz,prePIv,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [prefvvc,model.ind{indi}]=chain2c_tensor(pretilfz,pretilfzz,prePIv,prePIvvc,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [prefvvvc,model.ind{indi}]=chain3c_tensor(pretilfz,pretilfzz,pretilfzzz,...
        prePIv,prePIvvc,prePIvvvc,...
        model.prefvars_chain3c_M2,...
        model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;

    [ prefvv,model.ind{indi} ] = uncompressderivs( prefvvc,2,n_prefvars,model.fv(model.prefrows,prefvars),model.ind{indi} );
    indi=indi+1;

    [ prefvvv,model.ind{indi} ] = uncompressderivs( prefvvvc,3,n_prefvars,model.fv(model.prefrows,prefvars),model.ind{indi} );
    indi=indi+1;

elseif order==3
    indi=totindi+(model.pre_n-1)*(1+3+4)+1+4+6;
    for i=2:model.pre_n
        [prePIv,model.ind{indi}]=colsort(prePIv,model.ind{indi});
        indi=indi+1;
        [prePIvvc,model.ind{indi}]=colsort(prePIvvc,model.ind{indi});
        indi=indi+1;
        [prePIvvvvc,model.ind{indi}]=chain4c_tensor(prePIz,prePIzz,prePIzzz,prePIzzzz,...
            prePIv,prePIvvc,prePIvvvc,prePIvvvvc,...
            model.prefvars_chain4c_M2,model.prefvars_chain4c_M3,model.prefvars_chain4c_M4,...
            model.ind{indi},n_ind,maxload,'vec',model.prezz,model.maxzz);
        indi=indi+1;
        [prePIvvvc,model.ind{indi}]=chain3c_tensor(prePIz,prePIzz,prePIzzz,...
            prePIv,prePIvvc,prePIvvvc,model.prefvars_chain3c_M2,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [prePIvvc,model.ind{indi}]=chain2c_tensor(prePIz,prePIzz,prePIv,prePIvvc,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [prePIv,model.ind{indi}]=chain1_tensor(prePIz,prePIv,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    end
    [prePIv,model.ind{indi}]=colsort(prePIv,model.ind{indi});
    indi=indi+1;
    [prePIvvc,model.ind{indi}]=colsort(prePIvvc,model.ind{indi});
    indi=indi+1;

    [prefv,model.ind{indi}]=chain1_tensor(pretilfz,prePIv,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;

    [prefvvc,model.ind{indi}]=chain2c_tensor(pretilfz,pretilfzz,prePIv,prePIvvc,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [prefvvvc,model.ind{indi}]=chain3c_tensor(pretilfz,pretilfzz,pretilfzzz,...
        prePIv,prePIvvc,prePIvvvc,...
        model.prefvars_chain3c_M2,...
        model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [prefvvvvc,model.ind{indi}]=chain4c_tensor(pretilfz,pretilfzz,pretilfzzz,pretilfzzzz,...
        prePIv,prePIvvc,prePIvvvc,prePIvvvvc,...
        model.prefvars_chain4c_M2,model.prefvars_chain4c_M3,model.prefvars_chain4c_M4,...
        model.ind{indi},n_ind,maxload,'vec',model.pretilfz,model.maxtilfz);
    indi=indi+1;

    [ prefvv,model.ind{indi} ] = uncompressderivs( prefvvc,2,n_prefvars,model.fv(model.prefrows,prefvars),model.ind{indi} );
    indi=indi+1;

    [ prefvvv,model.ind{indi} ] = uncompressderivs( prefvvvc,3,n_prefvars,model.fv(model.prefrows,prefvars),model.ind{indi} );
    indi=indi+1;

    [ prefvvvv,model.ind{indi} ] = uncompressderivs( prefvvvvc,4,n_prefvars,model.fv(model.prefrows,prefvars),model.ind{indi} );
    indi=indi+1;
end
totindi=totindi+(model.pre_n-1)*(1+3+4+6)+1+4+6+9;
indi=totindi;


preR0=pref;

% compute R1 (first order derivatives of R w.r.t x)
if order>=1
    vx=vconcat(sptensor(n_y,n_x),sptensor(gx));
    vx=vconcat(vx,hx);
    vx=vconcat(vx,spteye(n_x));
    
    prevx=takerows(vx,prefvars);
    if isempty(pref)
        preR1=zeros(0,n_x);
    else
        [preR1,model.ind{indi}]=chain1_tensor(prefv,prevx,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        preR1=sptensor2spmat(preR1);
    end
end

% compute R2 (second order derivatives of R w.r.t x)
if order>=2
    vxxc=vconcat(sptensor(n_y,unique2),sptensor(gxxc));
    vxxc=vconcat(vxxc,hxxc);
    vxxc=vconcat(vxxc,sptensor(n_x,unique2));
    prevxxc=takerows(vxxc,prefvars);
    prevx=colsort(prevx);
    if isempty(pref)
        preR2=zeros(0,unique2);
    else
    [preR2,model.ind{indi}]=chain2c_tensor(prefv,prefvv,prevx,prevxxc,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    preR2=sptensor2spmat(preR2);
    end
end

% compute R3 (third order derivatives of R w.r.t x)
if order>=3
    vxxxc=vconcat(sptensor(n_y,unique3),sptensor(gxxxc));
    vxxxc=vconcat(vxxxc,hxxxc);
    vxxxc=vconcat(vxxxc,sptensor(n_x,unique3));
    prevxxxc=takerows(vxxxc,prefvars);

    if ~isfield(model,'x_chain3c_M2')
        [~,W2x]=create_UW(n_x,2);
        [U3x,~]=create_UW(n_x,3);
        OMEGAx=create_OMEGA(n_x,3);
        model.x_chain3c_M2=spmat2sptensor(kron(speye(n_x),W2x)*OMEGAx.OMEGA1*U3x);
        clear W2x U3x
    end
    
    if isempty(pref)
        preR3=zeros(0,unique3);
    else
        [preR3,model.ind{indi}]=chain3c_tensor(prefv,prefvv,prefvvv,prevx,prevxxc,prevxxxc,model.x_chain3c_M2,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        preR3=sptensor2spmat(preR3);
    end

end
totindi=totindi+3;
indi=totindi;
% Jacobian of the nonstochastic equations and derivatives

v_theta=vconcat(sptensor(n_y,n_theta),v_theta_lower_terms);
prev_theta=takerows(v_theta,prefvars);

[ preJ0,model.ind{indi} ] = chain0_theta_tensor( prefv,prev_theta,model.ind{indi},n_ind,maxload,'vec' );
indi=indi+1;
preJ0=sptensor2spmat(preJ0);

if order>=1
    vx_theta=vconcat(sptensor(n_y,[n_x,n_theta]),vx_theta_lower_terms);
    prevx_theta=takerows(vx_theta,prefvars);
    prevx_theta=unfold(prevx_theta);
    if isempty(pref)
        preJ1=zeros(0,n_x*n_theta);
    else
        [ preJ1,model.ind{indi} ] = chain1_theta_tensor( prefv,prefvv,...
            prevx,prev_theta,prevx_theta,...
            model.ind{indi},n_ind,maxload,'vec' );
        indi=indi+1;
        preJ1=sptensor2spmat(preJ1);
    end
end
if order>=2
    vxxc_theta=vconcat(sptensor(n_y,[unique2,n_theta]),vxxc_theta_lower_terms);
    prevxxc_theta=takerows(vxxc_theta,prefvars);
    prevxxc_theta=unfold(prevxxc_theta);

    if ~isfield(model,'chain2c_theta_M2')
        GAMMA=create_GAMMA(n_x,n_theta,2);
        U2=create_UW(n_x,2);
        model.chain2c_theta_M2=spmat2sptensor(GAMMA.GAMMA1*U2);
        clear GAMMA U2
    end
    
    if isempty(pref)
        preJ2=zeros(0,unique2*n_theta);
    else
    [ preJ2,model.ind{indi} ] = chain2c_theta_tensor( prefv,prefvv,prefvvv,...
        prevx,prevxxc,...
        prev_theta,prevx_theta,prevxxc_theta,...
        model.chain2c_theta_M2,...
        model.ind{indi},n_ind,maxload,'vec' );
    indi=indi+1;
    preJ2=sptensor2spmat(preJ2);
    end
end
if order>=3
    vxxxc_theta=vconcat(sptensor(n_y,[unique3,n_theta]),vxxxc_theta_lower_terms);
    prevxxxc_theta=takerows(vxxxc_theta,prefvars);
    prevxxxc_theta=unfold(prevxxxc_theta);

    if ~isfield(model,'chain3c_theta_M2')
        GAMMA=create_GAMMA(n_x,n_theta,3);
        [~,W2]=create_UW(n_x,2);
        U3=create_UW(n_x,3);
        model.chain3c_theta_M2=spmat2sptensor(kron(W2,speye(n_x))*GAMMA.GAMMA2*U3);
        model.chain3c_theta_M3=spmat2sptensor(kron(speye(n_x),W2)*GAMMA.GAMMA3*U3);
        model.chain3c_theta_M4=spmat2sptensor(kron(speye(n_x),W2)*GAMMA.GAMMA4*U3);
        tempind=1:unique2*n_x;
        tempmat=speye(unique2*n_x);
        model.chain3c_theta_M5=spmat2sptensor(tempmat(:,permute(reshape(tempind,n_x,unique2),[2,1]))*kron(speye(n_x),W2)*GAMMA.GAMMA5*U3);
        clear GAMMA W2 U3 GAMMA

    end
    
    if isempty(pref)
        preJ3=zeros(0,unique3*n_theta);
    else
        [ preJ3,model.ind{indi} ] = chain3c_theta_tensor( prefv,prefvv,prefvvv,prefvvvv,...
            prevx,prevxxc,prevxxxc,...
            prev_theta,prevx_theta,prevxxc_theta,prevxxxc_theta,...
            model.chain3c_theta_M2,model.chain3c_theta_M3,model.chain3c_theta_M4,model.chain3c_theta_M5,...
            model.ind{indi},n_ind,maxload,'vec' );
        indi=indi+1;
        preJ3=sptensor2spmat(preJ3);
    end
end
totindi=totindi+4;
indi=totindi;

% STEP 2: COMPUTE STOCHASTIC EQUATIONS AND DERIVATIVES, AND THEIR JACOBIAN,
% AND TAKE EXPECTED VALUE

% prepare vectorized expressions for the stochastic part. 

nx_vec=repmat(nx,1,n_nodes); % n_nodes=the number of nodes of the discrete shocks.
    
% next period state
nxp_vec=repmat(h,1,n_nodes)+eta*nep;
nxp_c0_vec=nxp_vec-repmat(c0,1,n_nodes);

nxp_c0=sptensor(nxp_c0_vec(:,1)); % convert to 1D tensor
nxp_c0.vals=nxp_c0_vec'; % assign the correct values in vals.

if ~isfield(model,'M2')
    model.M2=[];
    model.M3=[];
    model.M3x=[];
    model.TW2=[];
    model.TW3=[];
    if order>=2
        model.M2=fold(spmat2sptensor(model.W{2}*model.W{2}'*model.U{2}'),n_x,n_x);
        model.TW2=fold(spmat2sptensor(model.W{2}),n_x,n_x);
        model.TW2U2=spmat2sptensor(model.W{2}*model.U{2});
    end
    if order>=3
        model.M3=fold(spmat2sptensor(model.W{3}*model.W{3}'*model.U{3}'),n_x,n_x,n_x);
        model.M3x=fold(spmat2sptensor(3*model.W{3}*kron(model.W{2}'*model.U{2}',speye(n_x))),n_x,n_x,n_x);
        model.TW3=fold(spmat2sptensor(model.W{3}),n_x,n_x,n_x);
        model.TW3U3=spmat2sptensor(model.W{3}*model.U{3});
    end
end

if ~isfield(model,'M3')
    model.M3=[];
    model.M3x=[];
    model.TW3=[];
end

[Xp_vecT,Xpx_vecT,Xpxx_vecT,Xpxxx_vecT,Xpxxxx_vecT,model.ind{indi}]=create_X_tensor(order,nxp_c0,...
    model.M2,model.M3,model.M3x,model.TW2,model.TW3,unique2,unique3,model.ind{indi},n_ind,maxload,'vec');
indi=indi+1;

if strcmp(model.jacobian,'approximate')
      nx_c0=sptensor(n_x,1);
      nx_c0.vals=repmat(nx_c0.vals,n_nodes,1);
      [X_vecT,Xx_vecT,Xxx_vecT,Xxxx_vecT,Xxxxx_vecT,model.ind{indi}]=create_X_tensor(order,nx_c0,...
            model.M2,model.M3,model.M3x,model.TW2,model.TW3,unique2,unique3,model.ind{indi},n_ind,maxload,'vec');
      indi=indi+1;
end
totindi=totindi+2;
indi=totindi;


gp_vec=g_coeffs*Xp_vecT.vals';

stochg_coeffsT=sptensor(stochg_coeffs);
[stochgpx_vecT,model.ind{indi}]=contraction1(stochg_coeffsT,Xpx_vecT,model.ind{indi},n_ind,maxload,'vec');
indi=indi+1;

% control vars in t+1
nyp_vec=gp_vec;

% evaluate residuals R0
nv_vec=[nyp_vec;repmat(ny,1,n_nodes);nxp_vec;nx_vec];
nstochv_vec=nv_vec(stochfvars,:); % stochastic vars.

nstochu_vec=stochu_fun(nv_vec(model.stochuvars,:),params); % all stochastic u vars
nu_vec=repmat(nu,1,n_nodes);
nu_vec(model.stochurows,:)=nstochu_vec;
nz_vec=[nv_vec;nu_vec];
stochf_vec=stochtilf_fun(nz_vec(model.stochtilfzvars,:),params);

% compute derivatives of stochPI and stochtilf w.r.t z

stochPIz_vec=stochPI_d1(nz_vec',params,model.stochPI_ind_u);
prePIz_full.vals=repmat(prePIz_full.vals,n_nodes,1);

i1=tfind(prePIz_full);
i2=tfind(stochPIz_vec);
PIz_vec=sptensor([i1;i2],[prePIz_full.cols;stochPIz_vec.cols],[prePIz_full.vals,stochPIz_vec.vals],n_z,[n_z]);
stochtilfz_vec=stochtilf_tilf_d1([nv_vec(stochfvars,:);nu_vec(model.stochfuvars,:)]',params,model.stochtilf_ind_u);


if order>=1
    stochPIzz_vec=stochPI_d2(nz_vec',params,model.stochPI_ind_u);
    prePIzz_full.vals=repmat(prePIzz_full.vals,n_nodes,1);

    i1=tfind(prePIzz_full);
    i2=tfind(stochPIzz_vec);
    PIzz_vec=sptensor([i1;i2],[prePIzz_full.cols;stochPIzz_vec.cols],[prePIzz_full.vals,stochPIzz_vec.vals],n_z,[n_z,n_z]);
    stochtilfzz_vec=stochtilf_tilf_d2([nv_vec(stochfvars,:);nu_vec(model.stochfuvars,:)]',params,model.stochtilf_ind_u);
end
if order>=2
    stochPIzzz_vec=stochPI_d3(nz_vec',params,model.stochPI_ind_u);
    prePIzzz_full.vals=repmat(prePIzzz_full.vals,n_nodes,1);

    i1=tfind(prePIzzz_full);
    i2=tfind(stochPIzzz_vec);
    PIzzz_vec=sptensor([i1;i2],[prePIzzz_full.cols;stochPIzzz_vec.cols],[prePIzzz_full.vals,stochPIzzz_vec.vals],n_z,[n_z,n_z,n_z]);
    stochtilfzzz_vec=stochtilf_tilf_d3([nv_vec(stochfvars,:);nu_vec(model.stochfuvars,:)]',params,model.stochtilf_ind_u);
end
if order>=3
    stochPIzzzz_vec=stochPI_d4(nz_vec',params,model.stochPI_ind_u);
    prePIzzzz_full.vals=repmat(prePIzzzz_full.vals,n_nodes,1);

    i1=tfind(prePIzzzz_full);
    i2=tfind(stochPIzzzz_vec);
    PIzzzz_vec=sptensor([i1;i2],[prePIzzzz_full.cols;stochPIzzzz_vec.cols],[prePIzzzz_full.vals,stochPIzzzz_vec.vals],n_z,[n_z,n_z,n_z,n_z]);
    stochtilfzzzz_vec=stochtilf_tilf_d4([nv_vec(stochfvars,:);nu_vec(model.stochfuvars,:)]',params,model.stochtilf_ind_u);
end

totindi=totindi+1;
indi=totindi;

EstochR0=stochf_vec*P; %expected value of stochR0
h_thetaT=h_theta;
[chain0_theta_stochgpT,model.ind{indi}]=contraction1(stochgpx_vecT,h_thetaT,model.ind{indi},n_ind,maxload,'vec');
indi=indi+1;

totindi=totindi+1;
indi=totindi;

% compute derivatives of g w.r.t x at next period state (nxp)
if order>=1
    if order>1
        [stochgpxx_vecT,model.ind{indi}]=contraction1(stochg_coeffsT,Xpxx_vecT,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    else
        stochgpxx_vecT=sptensor(n_stochy,n_x^2,n_nodes);
    end
    hxT=hx;
    [chain1_stochgpT,model.ind{indi}]=contraction1(stochgpx_vecT,hxT,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;

    hx_thetaT=unfold(hx_theta);

    stochgpxx_vecT=fold(stochgpxx_vecT,n_x,n_x);
    [chain1_theta_stochgpT,model.ind{indi}]=chain1_theta_tensor(stochgpx_vecT,stochgpxx_vecT,hxT,h_thetaT,hx_thetaT,...
        model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
end
if order>=2
    if order>2
        [stochgpxxx_vecT,model.ind{indi}]=contraction1(stochg_coeffsT,Xpxxx_vecT,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    else
        stochgpxxx_vecT=sptensor(n_stochy,n_x^3,n_nodes);
    end
    hxxcT=hxxc;
    hxT=colsort(hxT);
    [chain2_stochgpT,model.ind{indi}]=chain2c_tensor(stochgpx_vecT,stochgpxx_vecT,hxT,hxxcT,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;

    hxxc_thetaT=unfold(hxxc_theta);
    stochgpxxx_vecT=fold(stochgpxxx_vecT,n_x,n_x,n_x);

    [chain2c_theta_stochgpT,model.ind{indi}]=chain2c_theta_tensor(stochgpx_vecT,stochgpxx_vecT,stochgpxxx_vecT,hxT,hxxcT,h_thetaT,hx_thetaT,hxxc_thetaT,...
        model.chain2c_theta_M2,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
end
if order>=3
    if order>3
        [stochgpxxxx_vecT,model.ind{indi}]=contraction1(stochg_coeffsT,Xpxxxx_vecT,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    else
        stochgpxxxx_vecT=sptensor(n_stochy,n_x^4,n_nodes);
    end
    hxxxcT=hxxxc;

    [chain3_stochgpT,model.ind{indi}]=chain3c_tensor(stochgpx_vecT,stochgpxx_vecT,stochgpxxx_vecT,...
        hxT,hxxcT,hxxxcT,...
        model.x_chain3c_M2,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;

    hxxxc_thetaT=unfold(hxxxc_theta);
    stochgpxxxx_vecT=fold(stochgpxxxx_vecT,n_x,n_x,n_x,n_x);

    [chain3c_theta_stochgpT,model.ind{indi}]=chain3c_theta_tensor(stochgpx_vecT,stochgpxx_vecT,stochgpxxx_vecT,stochgpxxxx_vecT,...
    hxT,hxxcT,hxxxcT,h_thetaT,hx_thetaT,hxxc_thetaT,hxxxc_thetaT,...
        model.chain3c_theta_M2,model.chain3c_theta_M3,model.chain3c_theta_M4,model.chain3c_theta_M5,...
        model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;   
end
totindi=totindi+9;
indi=totindi;


% extract PIv,Pivv,... from PIz,PIzz,...
if order>=0
    [stochPIz_vec,model.ind{indi}]=extract(PIz_vec,model.stochfzvars,model.stochfzvars,0,model.ind{indi});
    indi=indi+1;
    [stochPIv_vec,model.ind{indi}]=extract(PIz_vec,model.stochfzvars,model.stochfvars,0,model.ind{indi});
    indi=indi+1;  
end
if order>=1
    [stochPIzz_vec,model.ind{indi}]=extract(PIzz_vec,model.stochfzvars,model.stochfzvars,0,model.ind{indi});
    indi=indi+1;
    [stochPIvvc_vec,model.ind{indi}]=extract(PIzz_vec,model.stochfzvars,model.stochfvars,1,model.ind{indi});
    indi=indi+1;
end
if order>=2
    [stochPIzzz_vec,model.ind{indi}]=extract(PIzzz_vec,model.stochfzvars,model.stochfzvars,0,model.ind{indi});
    indi=indi+1;
    [stochPIvvvc_vec,model.ind{indi}]=extract(PIzzz_vec,model.stochfzvars,model.stochfvars,1,model.ind{indi});
    indi=indi+1;
    if ~isfield(model,'stochfvars_chain3c_M2')
        [ tempM ] = chainsM( n_stochfvars,3 );
        model.stochfvars_chain3c_M2=tempM{2};
        clear tempM
    end
end
if order>=3
    [stochPIzzzz_vec,model.ind{indi}]=extract(PIzzzz_vec,model.stochfzvars,model.stochfzvars,0,model.ind{indi});
    indi=indi+1;
    [stochPIvvvvc_vec,model.ind{indi}]=extract(PIzzzz_vec,model.stochfzvars,model.stochfvars,1,model.ind{indi});
    indi=indi+1;
    if ~isfield(model,'stochfvars_chain4c_M2')
        [ tempM ] = chainsM( n_stochfvars,4 );
        model.stochfvars_chain4c_M2=tempM{2};
        model.stochfvars_chain4c_M3=tempM{3};
        model.stochfvars_chain4c_M4=tempM{4};
        clear tempM
    end
end
totindi=totindi+8;
indi=totindi;


% compute derivatives of stochf w.r.t v by high order chain rules, and
% multiply by the weights P (probabilities).

if order==0
    indi=totindi;
    for i=2:model.stoch_n
        [stochPIv_vec,model.ind{indi}]=chain1_tensor(stochPIz_vec,stochPIv_vec,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    end
    [stochfv_vec,model.ind{indi}]=chain1_tensor(stochtilfz_vec,stochPIv_vec,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    stochfv_vec.vals=multcol(stochfv_vec.vals,P);
elseif order==1
    indi=totindi+(model.stoch_n-1)*(1)+1;
    for i=2:model.stoch_n
        [stochPIv_vec,model.ind{indi}]=colsort(stochPIv_vec,model.ind{indi});
        indi=indi+1;
        [stochPIvvc_vec,model.ind{indi}]=chain2c_tensor(stochPIz_vec,stochPIzz_vec,stochPIv_vec,stochPIvvc_vec,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [stochPIv_vec,model.ind{indi}]=chain1_tensor(stochPIz_vec,stochPIv_vec,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    end
    [stochPIv_vec,model.ind{indi}]=colsort(stochPIv_vec,model.ind{indi});
    indi=indi+1;

    [stochfv_vec,model.ind{indi}]=chain1_tensor(stochtilfz_vec,stochPIv_vec,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [stochfvvc_vec,model.ind{indi}]=chain2c_tensor(stochtilfz_vec,stochtilfzz_vec,stochPIv_vec,stochPIvvc_vec,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;

    stochfv_vec.vals=multcol(stochfv_vec.vals,P);
    stochfvvc_vec.vals=multcol(stochfvvc_vec.vals,P);

    [ stochfvv_vec,model.ind{indi} ] = uncompressderivs( stochfvvc_vec,2,n_stochfvars,model.fv(model.stochfrows,stochfvars),model.ind{indi} );
    indi=indi+1;
    clear stochfvvc_vec
    
elseif order==2
    indi=totindi+(model.stoch_n-1)*(1+3)+1+4;
    for i=2:model.stoch_n
        [stochPIv_vec,model.ind{indi}]=colsort(stochPIv_vec,model.ind{indi});
        indi=indi+1;
        [stochPIvvvc_vec,model.ind{indi}]=chain3c_tensor(stochPIz_vec,stochPIzz_vec,stochPIzzz_vec,...
            stochPIv_vec,stochPIvvc_vec,stochPIvvvc_vec,model.stochfvars_chain3c_M2,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [stochPIvvc_vec,model.ind{indi}]=chain2c_tensor(stochPIz_vec,stochPIzz_vec,stochPIv_vec,stochPIvvc_vec,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [stochPIv_vec,model.ind{indi}]=chain1_tensor(stochPIz_vec,stochPIv_vec,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    end
    [stochPIv_vec,model.ind{indi}]=colsort(stochPIv_vec,model.ind{indi});
    indi=indi+1;

    [stochfv_vec,model.ind{indi}]=chain1_tensor(stochtilfz_vec,stochPIv_vec,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [stochfvvc_vec,model.ind{indi}]=chain2c_tensor(stochtilfz_vec,stochtilfzz_vec,stochPIv_vec,stochPIvvc_vec,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [stochfvvvc_vec,model.ind{indi}]=chain3c_tensor(stochtilfz_vec,stochtilfzz_vec,stochtilfzzz_vec,...
        stochPIv_vec,stochPIvvc_vec,stochPIvvvc_vec,...
        model.stochfvars_chain3c_M2,...
        model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;

    stochfv_vec.vals=multcol(stochfv_vec.vals,P);
    stochfvvc_vec.vals=multcol(stochfvvc_vec.vals,P);
    stochfvvvc_vec.vals=multcol(stochfvvvc_vec.vals,P);

    [ stochfvv_vec,model.ind{indi} ] = uncompressderivs( stochfvvc_vec,2,n_stochfvars,model.fv(model.stochfrows,stochfvars),model.ind{indi} );
    indi=indi+1;
    clear stochfvvc_vec

    [ stochfvvv_vec,model.ind{indi} ] = uncompressderivs( stochfvvvc_vec,3,n_stochfvars,model.fv(model.stochfrows,stochfvars),model.ind{indi} );
    indi=indi+1;
    clear stochfvvvc_vec

elseif order==3
    indi=totindi+(model.stoch_n-1)*(1+3+4)+1+4+6;
    for i=2:model.stoch_n
        [stochPIv_vec,model.ind{indi}]=colsort(stochPIv_vec,model.ind{indi});
        indi=indi+1;
        [stochPIvvc_vec,model.ind{indi}]=colsort(stochPIvvc_vec,model.ind{indi});
        indi=indi+1;
        [stochPIvvvvc_vec,model.ind{indi}]=chain4c_tensor(stochPIz_vec,stochPIzz_vec,stochPIzzz_vec,stochPIzzzz_vec,...
            stochPIv_vec,stochPIvvc_vec,stochPIvvvc_vec,stochPIvvvvc_vec,...
            model.stochfvars_chain4c_M2,model.stochfvars_chain4c_M3,model.stochfvars_chain4c_M4,...
            model.ind{indi},n_ind,maxload,'vec',model.stochzz,model.maxzz);
        indi=indi+1;
        [stochPIvvvc_vec,model.ind{indi}]=chain3c_tensor(stochPIz_vec,stochPIzz_vec,stochPIzzz_vec,...
            stochPIv_vec,stochPIvvc_vec,stochPIvvvc_vec,model.stochfvars_chain3c_M2,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [stochPIvvc_vec,model.ind{indi}]=chain2c_tensor(stochPIz_vec,stochPIzz_vec,stochPIv_vec,stochPIvvc_vec,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
        [stochPIv_vec,model.ind{indi}]=chain1_tensor(stochPIz_vec,stochPIv_vec,model.ind{indi},n_ind,maxload,'vec');
        indi=indi+1;
    end
    [stochPIv_vec,model.ind{indi}]=colsort(stochPIv_vec,model.ind{indi});
    indi=indi+1;
    [stochPIvvc_vec,model.ind{indi}]=colsort(stochPIvvc_vec,model.ind{indi});
    indi=indi+1;

    [stochfv_vec,model.ind{indi}]=chain1_tensor(stochtilfz_vec,stochPIv_vec,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;

    [stochfvvc_vec,model.ind{indi}]=chain2c_tensor(stochtilfz_vec,stochtilfzz_vec,stochPIv_vec,stochPIvvc_vec,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [stochfvvvc_vec,model.ind{indi}]=chain3c_tensor(stochtilfz_vec,stochtilfzz_vec,stochtilfzzz_vec,...
        stochPIv_vec,stochPIvvc_vec,stochPIvvvc_vec,...
        model.stochfvars_chain3c_M2,...
        model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    [stochfvvvvc_vec,model.ind{indi}]=chain4c_tensor(stochtilfz_vec,stochtilfzz_vec,stochtilfzzz_vec,stochtilfzzzz_vec,...
        stochPIv_vec,stochPIvvc_vec,stochPIvvvc_vec,stochPIvvvvc_vec,...
        model.stochfvars_chain4c_M2,model.stochfvars_chain4c_M3,model.stochfvars_chain4c_M4,...
        model.ind{indi},n_ind,maxload,'vec',model.stochtilfz,model.maxtilfz);
    indi=indi+1;

    stochfv_vec.vals=multcol(stochfv_vec.vals,P);
    stochfvvc_vec.vals=multcol(stochfvvc_vec.vals,P);
    stochfvvvc_vec.vals=multcol(stochfvvvc_vec.vals,P);
    stochfvvvvc_vec.vals=multcol(stochfvvvvc_vec.vals,P);

    [ stochfvv_vec,model.ind{indi} ] = uncompressderivs( stochfvvc_vec,2,n_stochfvars,model.fv(model.stochfrows,stochfvars),model.ind{indi} );
    indi=indi+1;
    clear stochfvvc_vec

    [ stochfvvv_vec,model.ind{indi} ] = uncompressderivs( stochfvvvc_vec,3,n_stochfvars,model.fv(model.stochfrows,stochfvars),model.ind{indi} );
    indi=indi+1;
    clear stochfvvvc_vec

    [ stochfvvvv_vec,model.ind{indi} ] = uncompressderivs( stochfvvvvc_vec,4,n_stochfvars,model.fv(model.stochfrows,stochfvars),model.ind{indi} );
    indi=indi+1;
    clear stochfvvvvc_vec

end
totindi=totindi+8+(model.stoch_n-1)*(1+3+4+6)+1+4+6+9;
indi=totindi;

temp=reshape(1:n_theta,n_f,n_b);
stochtheta=temp(stochy,:);
stochtheta=stochtheta(:);
n_stochtheta=length(stochtheta);

% exact/approximate Jacobian of gp,gpx,... 
IstochyT=spteye(n_stochy);
if strcmp(model.jacobian,'exact')
      [stochgp_stochtheta_vecT,model.ind{indi}]=tkron(ttranspose(Xp_vecT),IstochyT,model.ind{indi},n_ind,maxload,'vec');
      indi=indi+1;
else
      [stochgp_stochtheta_vecT,model.ind{indi}]=tkron(ttranspose(X_vecT),IstochyT,model.ind{indi},n_ind,maxload,'vec');
      indi=indi+1;
end

nnzgp_theta_vecT=unfold(changecols(changerows(stochgp_stochtheta_vecT,stochy,n_y),stochy,n_f,1)); %n_y,n_theta with only nnz values
nnzchain0_theta_gpT=changerows(chain0_theta_stochgpT,stochy,n_y);

v_theta_lower_terms_vecT=v_theta_lower_terms;
v_theta_lower_terms_vecT.vals=repmat(v_theta_lower_terms_vecT.vals,n_nodes,1);

v_theta_vecT=vconcat(tplus(nnzchain0_theta_gpT,nnzgp_theta_vecT,maxload),v_theta_lower_terms_vecT);
clear nnzchain0_theta_gpT nnzgp_theta_vecT
stochv_theta_vecT=takerows(v_theta_vecT,stochfvars);

clear v_theta_vecT

totindi=totindi+2;
indi=totindi;
% derivatives of v w.r.t x
if order>=1
      vxT=vx; 
      stochvxT=takerows(vxT,stochfvars(n_stochy+1:end));
      stochvxT.vals=repmat(stochvxT.vals,n_nodes,1);
      stochvxT=vconcat(chain1_stochgpT,stochvxT);

      if strcmp(model.jacobian,'exact')
            [stochgpx_stochtheta_vecT,model.ind{indi}]=tkron(ttranspose(Xpx_vecT),IstochyT,model.ind{indi},n_ind,maxload,'vec');
            indi=indi+1;
      else
            [stochgpx_stochtheta_vecT,model.ind{indi}]=tkron(ttranspose(Xx_vecT),IstochyT,model.ind{indi},n_ind,maxload,'vec');
            indi=indi+1;
      end
      stochgpx_stochtheta_vecT=ptr1d(col2ptr(ptr2col(ptr2d(unfold(stochgpx_stochtheta_vecT),n_stochy,n_x),1),2));%n_stochy*n_stochtheta,n_x
      %
      [chain1_gp_theta_vecT,model.ind{indi}]=contraction1(stochgpx_stochtheta_vecT,hxT,model.ind{indi},n_ind,maxload,'vec');    
      indi=indi+1;

      chain1_gp_theta_vecT=ptr2col(ptr2d(chain1_gp_theta_vecT,n_stochy,n_stochtheta),2);%n_stochy,n_x,n_stochtheta
      nnzchain1_gp_theta_vecT=col2ptr(changecols(changerows(chain1_gp_theta_vecT,stochy,n_y),stochtheta,n_theta,2),1); %n_y*n_x,n_theta with only nnz values

      nnzchain1_theta_gpT=col2ptr(changerows(chain1_theta_stochgpT,stochy,n_y),1);% n_y,n_x,n_theta

      vx_theta_lower_terms_vecT=vx_theta_lower_terms;
      vx_theta_lower_terms_vecT.vals=repmat(vx_theta_lower_terms_vecT.vals,n_nodes,1);

      vx_theta_upper_terms_vecT=unfold(ptr2col(tplus(nnzchain1_theta_gpT,nnzchain1_gp_theta_vecT,maxload),1));
      clear nnzchain1_theta_gpT nnzchain1_gp_theta_vecT
      vx_theta_vecT=vconcat(vx_theta_upper_terms_vecT,unfold(vx_theta_lower_terms_vecT));
      
      clear vx_theta_upper_terms_vecT vx_theta_lower_terms_vecT
      
      stochvx_theta_vecT=takerows(vx_theta_vecT,stochfvars);
      
      clear vx_theta_vecT
end
if order>=2
      vxxcT=vxxc;
      stochvxxcT=takerows(vxxcT,stochfvars(n_stochy+1:end));
      stochvxxcT.vals=repmat(stochvxxcT.vals,n_nodes,1);
      stochvxxcT=vconcat(chain2_stochgpT,stochvxxcT);

      if strcmp(model.jacobian,'exact')
            [stochgpxx_stochtheta_vecT,model.ind{indi}]=tkron(ttranspose(Xpxx_vecT),IstochyT,model.ind{indi},n_ind,maxload,'vec');
            indi=indi+1;
      else
            [stochgpxx_stochtheta_vecT,model.ind{indi}]=tkron(ttranspose(Xxx_vecT),IstochyT,model.ind{indi},n_ind,maxload,'vec');
            indi=indi+1;
      end

      stochgpxx_stochtheta_vecT=ptr1d(col2ptr(ptr2col(ptr2d(unfold(stochgpxx_stochtheta_vecT),n_stochy,n_x^2),1),2));%n_stochy*n_stochtheta,n_x^2
      stochgpxx_stochtheta_vecT=fold(stochgpxx_stochtheta_vecT,n_x,n_x);

      [chain2c_stochgp_stochtheta_vecT,model.ind{indi}]=chain2c_tensor(stochgpx_stochtheta_vecT,stochgpxx_stochtheta_vecT,...
        hxT,hxxcT,model.ind{indi},n_ind,maxload,'vec');    
      indi=indi+1;

      chain2c_stochgp_stochtheta_vecT=ptr2col(ptr2d(chain2c_stochgp_stochtheta_vecT,n_stochy,n_stochtheta),2);%n_stochy,unique2,n_stochtheta
      nnzchain2c_gp_theta_vecT=col2ptr(changecols(changerows(chain2c_stochgp_stochtheta_vecT,stochy,n_y),stochtheta,n_theta,2),1); %n_y*unique2,n_theta with only nnz values
      clear chain2c_stochgp_stochtheta_vecT

      nnzchain2c_theta_gpT=col2ptr(changerows(chain2c_theta_stochgpT,stochy,n_y),1);% n_y,unique2,n_theta

      vxxc_theta_lower_terms_vecT=unfold(vxxc_theta_lower_terms);
      vxxc_theta_lower_terms_vecT.vals=repmat(vxxc_theta_lower_terms_vecT.vals,n_nodes,1);

      vxxc_theta_upper_terms_vecT=unfold(ptr2col(tplus(nnzchain2c_theta_gpT,nnzchain2c_gp_theta_vecT,maxload),1));
      clear nnzchain2c_theta_gpT nnzchain2c_gp_theta_vecT
      vxxc_theta_vecT=vconcat(vxxc_theta_upper_terms_vecT,vxxc_theta_lower_terms_vecT);
      
      clear vxxc_theta_upper_terms_vecT vxxc_theta_lower_terms_vecT
      
      stochvxxc_theta_vecT=takerows(vxxc_theta_vecT,stochfvars);  
      
      clear vxxc_theta_vecT
end
if order>=3
      vxxxcT=vxxxc;
      stochvxxxcT=takerows(vxxxcT,stochfvars(n_stochy+1:end));
      stochvxxxcT.vals=repmat(stochvxxxcT.vals,n_nodes,1);
      stochvxxxcT=vconcat(chain3_stochgpT,stochvxxxcT);

      if strcmp(model.jacobian,'exact')
            [stochgpxxx_stochtheta_vecT,model.ind{indi}]=tkron(ttranspose(Xpxxx_vecT),IstochyT,model.ind{indi},n_ind,maxload,'vec');
            indi=indi+1;
      else
            [stochgpxxx_stochtheta_vecT,model.ind{indi}]=tkron(ttranspose(Xxxx_vecT),IstochyT,model.ind{indi},n_ind,maxload,'vec');
            indi=indi+1;
      end

      stochgpxxx_stochtheta_vecT=ptr1d(col2ptr(ptr2col(ptr2d(unfold(stochgpxxx_stochtheta_vecT),n_stochy,n_x^3),1),2));%n_stochy*n_stochtheta,n_x^2
      stochgpxxx_stochtheta_vecT=fold(stochgpxxx_stochtheta_vecT,n_x,n_x,n_x);

      [chain3c_stochgp_stochtheta_vecT,model.ind{indi}]=chain3c_tensor(stochgpx_stochtheta_vecT,stochgpxx_stochtheta_vecT,stochgpxxx_stochtheta_vecT,...
            hxT,hxxcT,hxxxcT,...
            model.x_chain3c_M2,...
            model.ind{indi},n_ind,maxload,'vec');    
      indi=indi+1;

      chain3c_stochgp_stochtheta_vecT=ptr2col(ptr2d(chain3c_stochgp_stochtheta_vecT,n_stochy,n_stochtheta),2);%n_stochy,unique3,n_stochtheta
      nnzchain3c_gp_theta_vecT=col2ptr(changecols(changerows(chain3c_stochgp_stochtheta_vecT,stochy,n_y),stochtheta,n_theta,2),1); %n_y*unique3,n_theta with only nnz values
      
      clear chain3c_stochgp_stochtheta_vecT
      
      nnzchain3c_theta_gpT=col2ptr(changerows(chain3c_theta_stochgpT,stochy,n_y),1);% n_y,unique3,n_theta

      clear chain3c_theta_stochgpT

      vxxxc_theta_lower_terms_vecT=unfold(vxxxc_theta_lower_terms);
      vxxxc_theta_lower_terms_vecT.vals=repmat(vxxxc_theta_lower_terms_vecT.vals,n_nodes,1);

      vxxxc_theta_upper_terms_vecT=unfold(ptr2col(tplus(nnzchain3c_theta_gpT,nnzchain3c_gp_theta_vecT,maxload),1));
      
      clear nnzchain3c_theta_gpT nnzchain3c_gp_theta_vecT
      vxxxc_theta_vecT=vconcat(vxxxc_theta_upper_terms_vecT,vxxxc_theta_lower_terms_vecT);
      
      clear vxxxc_theta_upper_terms_vecT vxxxc_theta_lower_terms_vecT
      
      stochvxxxc_theta_vecT=takerows(vxxxc_theta_vecT,stochfvars);    
      
      clear vxxxc_theta_vecT
end
totindi=totindi+9;
indi=totindi;

[ EstochJ0,model.ind{indi} ] = chain0_theta_tensor( stochfv_vec,stochv_theta_vecT,model.ind{indi},n_ind,maxload,'sum' );
indi=indi+1;

EstochJ0=sptensor2spmat(EstochJ0);

% EstochR1,EstochR2,EstochR3,EstochJ1,EstochJ2,EstochJ3
if order>=1
    [ EstochR1,model.ind{indi} ] = chain1_tensor( stochfv_vec,stochvxT,model.ind{indi},n_ind,maxload,'sum' );
    indi=indi+1;
    EstochR1=sptensor2spmat(EstochR1);

    [ EstochJ1,model.ind{indi} ] = chain1_theta_tensor( stochfv_vec,stochfvv_vec,...
        stochvxT,stochv_theta_vecT,stochvx_theta_vecT,...
        model.ind{indi},n_ind,maxload,'sum' );
    indi=indi+1;
    EstochJ1=sptensor2spmat(EstochJ1);
end
if order>=2
    stochvxT=colsort(stochvxT);
    [ EstochR2,model.ind{indi} ] = chain2c_tensor( stochfv_vec,stochfvv_vec,...
        stochvxT,stochvxxcT,model.ind{indi},n_ind,maxload,'sum' );
    indi=indi+1;
    EstochR2=sptensor2spmat(EstochR2);

    [ EstochJ2,model.ind{indi} ] = chain2c_theta_tensor( stochfv_vec,stochfvv_vec,stochfvvv_vec,...
        stochvxT,stochvxxcT,...
        stochv_theta_vecT,stochvx_theta_vecT,stochvxxc_theta_vecT,...
        model.chain2c_theta_M2,...
        model.ind{indi},n_ind,maxload,'sum' );
    indi=indi+1;
    EstochJ2=sptensor2spmat(EstochJ2);
    
end
if order>=3
    [ EstochR3,model.ind{indi} ] = chain3c_tensor( stochfv_vec,stochfvv_vec,stochfvvv_vec,...
        stochvxT,stochvxxcT,stochvxxxcT,...
        model.x_chain3c_M2,...
        model.ind{indi},n_ind,maxload,'sum' );
    indi=indi+1;
    EstochR3=sptensor2spmat(EstochR3);
    
    [ EstochJ3,model.ind{indi} ] = chain3c_theta_tensor( stochfv_vec,stochfvv_vec,stochfvvv_vec,stochfvvvv_vec,...
        stochvxT,stochvxxcT,stochvxxxcT,...
        stochv_theta_vecT,stochvx_theta_vecT,stochvxxc_theta_vecT,stochvxxxc_theta_vecT,...
        model.chain3c_theta_M2,model.chain3c_theta_M3,model.chain3c_theta_M4,model.chain3c_theta_M5,...
        model.ind{indi},n_ind,maxload,'sum' );
    indi=indi+1;
    EstochJ3=sptensor2spmat(EstochJ3);
end
totindi=totindi+7;

% build the system and the Jacobian

R0=spalloc(n_f,1,n_f);
J0=spalloc(n_f,n_theta,n_f*n_theta);

R0(model.prefrows)=preR0;
J0(model.prefrows,:)=preJ0;
R0(model.stochfrows)=EstochR0;
J0(model.stochfrows,:)=EstochJ0;

T=R0;
J=J0;

if order>=1
    R1=spalloc(n_f,n_x,nnz(preR1)+nnz(EstochR1));
    J1=spalloc(n_f,n_x*n_theta,nnz(preJ1)+nnz(EstochJ1));
    R1(model.prefrows,:)=preR1;
    J1(model.prefrows,:)=preJ1;
    R1(model.stochfrows,:)=EstochR1;
    J1(model.stochfrows,:)=EstochJ1;
    
    T=[T;R1(:)];
    J=[J; reshape(J1,n_f*n_x,n_theta)];
end
if order>=2
    R2=spalloc(n_f,size(preR2,2),nnz(preR2)+nnz(EstochR2));
    J2=spalloc(n_f,size(preJ2,2),nnz(preJ2)+nnz(EstochJ2));
    R2(model.prefrows,:)=preR2;
    J2(model.prefrows,:)=preJ2;
    R2(model.stochfrows,:)=EstochR2;
    J2(model.stochfrows,:)=EstochJ2;
    
    T=[T;R2(:)];
    J=[J;reshape(J2,n_f*unique2,n_theta)];
end
if order>=3
    R3=spalloc(n_f,size(preR3,2),nnz(preR3)+nnz(EstochR3));
    J3=spalloc(n_f,size(preJ3,2),nnz(preJ3)+nnz(EstochJ3));
    R3(model.prefrows,:)=preR3;
    J3(model.prefrows,:)=preJ3;
    R3(model.stochfrows,:)=EstochR3;
    J3(model.stochfrows,:)=EstochJ3;
    
    T=[T;R3(:)];
    J=[J;reshape(J3,n_f*unique3,n_theta)];
end


end