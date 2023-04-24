function [R0,g,h1,nu_vec]=eval_eqs(coeffs,x0,model,params,eta,c0,nep,P,N,poly)
% model is used to calculate the model equations and it is common for all solutions.
% if basis function is a power series, model should include compression
% matrices.

params=params(:);
P=P(:);

n_f=model.n_f;
n_x=model.n_x;
n_y=model.n_y;
% n_x2=model.n_x2;
n_x1=model.n_x1;
% n_v=model.n_v;

n_theta=numel(coeffs);
nparams=n_theta/n_f;

n_s=size(nep,2);

% pref_ind=model.pref_ind;
% stochf_ind=model.stochf_ind;

prefvars=model.prefvars;
stochfvars=model.stochfvars;

% Phi_ind=model.Phi_ind;

% current state
nx=x0;

if strcmp(poly.type,'power')
    % evaluate policy functions and their derivatives at current state
    % [g,h1]=fun_gh_model(nx,c0,DELTA,coeffs);
    xpowers=[1;x0-c0];
    lastkron=x0-c0;
    for i=2:N
        lastkron=kron(lastkron,x0-c0);
        xpowers=[xpowers;model.W{i}*lastkron];
    end
    gh_coeffs=reshape(coeffs,n_f,nparams);
    g_coeffs=gh_coeffs(1:n_y,:);
    gh1=gh_coeffs*xpowers;
elseif strcmp(poly.type,'smol')
    gh_coeffs=reshape(coeffs,n_f,nparams);
    addpath(poly.Tcheb_fun_folder);
%     Tcheb=feval(poly.Tcheb_fun,(x0-c0)./poly.DELTA,[]);
    Tcheb=Tcheb_fun((x0-c0)./poly.DELTA,[]);
    g_coeffs=gh_coeffs(1:n_y,:);
    gh1=gh_coeffs*Tcheb;
elseif strcmp(poly.type,'cheb')
%     x=[loghatkstar;
%     xvar;
%     logd;
%     logdback;
%     u];
    MU=params(9);
    u=x0(5);
    xvar=x0(2);
    logdback=x0(4);
    loghatkstar=x0(1);
    logd=x0(3);
    loghatz=MU+u+xvar*logdback;
    logtilk=loghatkstar-loghatz+xvar*logdback;
    x0_3D=[logtilk;logd;loghatz]; % 3D state-space
    Tcheb=feval(poly.Tcheb_fun,(x0_3D-c0)./poly.DELTA,[]);
    gh_coeffs=reshape(coeffs,poly.n_f,poly.nparams);
    g_coeffs=gh_coeffs(1:n_y,:);
    gh1=gh_coeffs*Tcheb;
    
    % need to add loghatkstarp
    hatkstarp=hatkstarp_fun([gh1(1:n_y);x0_3D],params);
    gh1=[gh1;log(hatkstarp)];
    
end

g=gh1(1:n_y,1);
h1=gh1(n_y+1:end,1);

% control vars
ny=g;

% expected values: 

% The residual function R is calculated in two steps. The first step
% calculates the nonstochastic (predetermined) rows of f (pref). The second step calculates
% the expected value of the stochastic rows of f (stochf).

% Step 1. predetermined equations (pref)

% preallocate
h=zeros(n_x,1);

% build h(x)
h(1:n_x1,1)=h1; % predetermined endogenous state vars
h(n_x1+1:end,1)=Phi_fun(nx,params); % expected value of exogenous state vars. shocks are added later

% next period state
nxp=h;

% evaluate residuals R0
nv=[zeros(n_y,1);ny;nxp;nx]; % all variables with zeros for the stochastic vars.

n_u=model.n_u;
nu=zeros(n_u,1);
npreu=preu_fun(nv(model.preuvars),params); % all predetermined u
nu(model.preurows)=npreu;
nz=[nv;nu];
pref=pretilf_fun(nz(model.pretilfzvars),params); % since i already have z, i use pretilf to evaluate f.

preR0=pref;

% Step 2. expected value of stochastic equations (stochf) 

% vectorized expressions

nx_vec=repmat(nx,1,n_s);
    
% next period state
nxp_vec=repmat(h,1,n_s)+eta*nep;
% nxp_c0_vec=nxp_vec-repmat(c0,1,n_s);

% evaluate g(x) at next period values
% [gp_vec]=fun_gp_model(nxp_vec,c0,DELTA,g_coeffs);



if strcmp(poly.type,'power')
    Xp_c0_vec=zeros(nparams,n_s);
    for j=1:n_s
        xpowers=[1;nxp_vec(:,j)-c0];
        lastkron=nxp_vec(:,j)-c0;
        for i=2:N
            lastkron=kron(lastkron,nxp_vec(:,j)-c0);
            xpowers=[xpowers;model.W{i}*lastkron];
        end
        Xp_c0_vec(:,j)=xpowers;
    end
    gp_vec=g_coeffs*Xp_c0_vec;
elseif strcmp(poly.type,'smol')
    Tchebp=zeros(length(Tcheb),n_s);
    for j=1:n_s
        Tchebp(:,j)=feval('Tcheb_fun',(nxp_vec(:,j)-c0)./poly.DELTA,[]);
    end
    rmpath(poly.Tcheb_fun_folder);
    gp_vec=g_coeffs*Tchebp;
elseif strcmp(poly.type,'cheb')
%     x=[loghatkstar;
%     xvar;
%     logd;
%     logdback;
%     u];
    Tchebp=zeros(length(Tcheb),n_s);
    MU=params(9);
    for j=1:n_s
        up=nxp_vec(5,j);
        xvarp=nxp_vec(2,j);
        logdbackp=nxp_vec(4,j);
        loghatkstarp=nxp_vec(1,j);
        logdp=nxp_vec(3,j);
        loghatzp=MU+up+xvarp*logdbackp;
        logtilkp=loghatkstarp-loghatzp+xvarp*logdbackp;
        xp_3D=[logtilkp;logdp;loghatzp]; % 3D state-space
        Tchebp(:,j)=feval(poly.Tcheb_fun,(xp_3D-c0)./poly.DELTA,[]);
    end
    gp_vec=g_coeffs*Tchebp;
end

% control vars in t+1
nyp_vec=gp_vec; 

% evaluate residuals R0
nv_vec=[nyp_vec;repmat(ny,1,n_s);nxp_vec;nx_vec];
nstochv_vec=nv_vec(stochfvars,:); % stochastic vars.

nstochu_vec=stochu_fun(nv_vec(model.stochuvars,:),params); % all stochastic u vars
nu_vec=repmat(nu,1,n_s);
nu_vec(model.stochurows,:)=nstochu_vec;
nz_vec=[nv_vec;nu_vec];
stochf_vec=stochtilf_fun(nz_vec(model.stochtilfzvars,:),params);

EstochR0=stochf_vec*P;
R0=zeros(n_f,1);
R0(model.prefrows)=preR0;
R0(model.stochfrows)=EstochR0;