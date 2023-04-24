function [R0,g,h1,nPhi,nu_vec]=eval_model(coeffs,x0_vec,model,params,eta,c0,nep,P)
% The function evaluates the model residuals at points x0_vec, and returns
% the residuals (R0), policy functions of control variables (g) and
% endogeneous state variables (h1), and expected exogeneous state variables
% (nPhi)
%
% © Copyright, Oren Levintal, June 13, 2016.

order=model.order(1);

params=params(:);
P=P(:);

n_f=model.n_f;
n_x=model.n_x;
n_y=model.n_y;
n_x2=model.n_x2;
n_x1=model.n_x1;
n_v=model.n_v;
n_theta=model.n_theta;
n_b=model.n_b;
n_s=size(nep,2);

prefvars=model.prefvars;
stochfvars=model.stochfvars;


% current state
nx=x0_vec;
n_grid=size(x0_vec,2);

nx_c0_vec=nx-repmat(c0,1,n_grid);
nx_c0=sptensor(nx_c0_vec(:,1)); % convert to 1D tensor
nx_c0.vals=nx_c0_vec'; % assign the correct values in vals.

[X_vecT]=create_X_tensor_no_derivs(order,nx_c0,...
    model.M2,model.M3,[],model.n_ind,n_grid,'vec');

gh_coeffs=reshape(coeffs,n_f,n_b);
g_coeffs=gh_coeffs(1:n_y,:);
gh1=gh_coeffs*X_vecT.vals';
g=gh1(1:n_y,:);
h1=gh1(n_y+1:end,:);

% control vars
ny=g;

% expected values: 

% The residual function R is calculated in two steps. The first step
% calculates the nonstochastic (predetermined) rows of f (pref). The second step calculates
% the expected value of the stochastic rows of f (stochf).

% Step 1. predetermined equations (pref)

% preallocate
h=zeros(n_x,n_grid);

% build h(x)
h(1:n_x1,:)=h1; % predetermined endogenous state vars
nPhi=Phi_fun(nx,params);
h(n_x1+1:end,:)=nPhi; % expected value of exogenous state vars. shocks are added later

% next period state
nxp=h;

% evaluate residuals R0
nv=[zeros(n_y,n_grid);ny;nxp;nx]; % all variables with zeros for the stochastic vars.

n_u=model.n_u;
nu=zeros(n_u,n_grid);
npreu=preu_fun(nv(model.preuvars,:),params); % all predetermined u
nu(model.preurows,:)=npreu;
nz=[nv;nu];
pref=pretilf_fun(nz(model.pretilfzvars,:),params); % since i already have z, i use pretilf to evaluate f.

preR0=pref;

% Step 2. expected value of stochastic equations (stochf) 

% vectorized expressions

nx_vec=repmat(nx,1,n_s);
    
% next period state
nxp_vec=repmat(h,1,n_s)+reshape(repmat(eta*nep,n_grid,1),n_x,n_grid*n_s);
nxp_c0_vec=nxp_vec-repmat(c0,1,n_grid*n_s);

nxp_c0=sptensor(nxp_c0_vec(:,1)); % convert to 1D tensor
nxp_c0.vals=nxp_c0_vec'; % assign the correct values in vals.

[Xp_vecT]=create_X_tensor_no_derivs(order,nxp_c0,...
    model.M2,model.M3,[],model.n_ind,n_grid*n_s,'vec');


gp_vec=g_coeffs*Xp_vecT.vals';

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


EstochR0=reshape(reshape(stochf_vec,[],n_s)*P,[],n_grid);
R0=zeros(n_f,n_grid);
R0(model.prefrows,:)=preR0;
R0(model.stochfrows,:)=EstochR0;

nu_vec=reshape(nu_vec,[],n_grid,n_s);