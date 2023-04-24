function [R,J,model]=fun_collocation(coeffs,model,params,x_grid,c0_vec,DELTA_vec,nep,P,Tcheb,...
    g_theta,h_theta,Phi,eta) 

% The function computes the residual functiona and Jacobian over a grid of
% points.
%
% © Copyright, Jesus Fernandez-Villaverde and Oren Levintal, June 13, 2016.

if ~isfield(model,'ind')
    model.ind=cell(200,1);
    indi=1;
else
    indi=2;
end

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

N=0;

params=params(:);
P=P(:);


n_f=model.n_f;
n_x=model.n_x;
n_y=model.n_y;
n_x2=model.n_x2;
n_x1=model.n_x1;
n_v=model.n_v;
n_theta=model.n_theta;
nparams=model.nparams;
n_s=size(nep,2);
n_z=model.n_z;
coeffs=reshape(coeffs,n_f,nparams);
g_coeffs=coeffs(1:n_y,:);

tempv=zeros(n_v,1);
tempv(model.stochfvars)=model.stochfvars;
tempv(n_y+1:end)=0;
stochy=nonzeros(tempv);
n_stochy=length(stochy);
stochg_coeffs=coeffs(stochy,:); 

prefvars=model.prefvars;
stochfvars=model.stochfvars;
n_prefvars=length(prefvars);
n_stochfvars=length(stochfvars);

n_prefrows=length(model.prefrows);
n_stochfrows=length(model.stochfrows);


% lower terms of v_theta
grid_points=size(x_grid,2);
v_theta_lower_terms=vconcat(vconcat(g_theta,h_theta),sptensor(n_x,n_theta,grid_points));


coeffs=reshape(coeffs,n_f,nparams);

gh1=coeffs*Tcheb;

ny=gh1(1:n_y,:);
h=[gh1(n_y+1:end,:);Phi];

nv=zeros(n_v,grid_points);

nv(2*n_y+n_x+1:end,:)=x_grid;
nv(n_y+1:2*n_y,:)=ny;
nv(2*n_y+1:2*n_y+n_x,:)=h;

% current state
nx=x_grid;



nxp=h;

% evaluate residuals R0

n_u=model.n_u;
nu=zeros(n_u,grid_points);
npreu=preu_fun(nv(model.preuvars,:),params); % all predetermined u
nu(model.preurows,:)=npreu;
nz=[nv;nu];
pref=pretilf_fun(nz(model.pretilfzvars,:),params); % since i already have z, i use pretilf to evaluate f.

prePIz_full=prePI_d1(nz',params,model.prePI_ind_u);
pretilfz=pretilf_tilf_d1([nv(prefvars,:);nu(model.prefuvars,:)]',params,model.pretilf_ind_u);

[prePIz,model.ind{indi}]=extract(prePIz_full,model.prefzvars,model.prefzvars,0,model.ind{indi});
indi=indi+1;
[prePIv,model.ind{indi}]=extract(prePIz_full,model.prefzvars,model.prefvars,0,model.ind{indi});
indi=indi+1;

for i=2:model.pre_n
    [prePIv,model.ind{indi}]=chain1_tensor(prePIz,prePIv,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
end
[prefv,model.ind{indi}]=chain1_tensor(pretilfz,prePIv,model.ind{indi},n_ind,maxload,'vec');
indi=indi+1;

preR0=pref;

% Jacobian

v_theta=vconcat(sptensor(n_y,n_theta,grid_points),v_theta_lower_terms);
prev_theta=takerows(v_theta,prefvars);

[ preJ0,model.ind{indi} ] = chain0_theta_tensor( prefv,prev_theta,model.ind{indi},n_ind,maxload,'vec' );
indi=indi+1;

preJ0_temp=preJ0;
preJ0_temp.vals=preJ0_temp.vals(1,:);
[i,j]=tfind(preJ0_temp);
i=double(i);
j=double(j);
vals=preJ0.vals';
vals=vals(:);
preJ0=sparse(repmat(i,grid_points,1)+kron((0:grid_points-1)',n_prefrows*ones(length(i),1)),...
    repmat(j,grid_points,1),...
    vals,...
    n_prefrows*grid_points,n_theta);
clear vals preJ0_temp

% vectorized expressions

nx_vec=repmat(nx,1,n_s);
    
% next period state
nxp_vec=repmat(h,1,n_s)+reshape(repmat(eta*nep,grid_points,1),n_x,[]);
nxp_c0_vec=nxp_vec-c0_vec;

nxp_c0=sptensor(nxp_c0_vec(:,1)); % convert to 1D tensor
nxp_c0.vals=nxp_c0_vec'; % assign the correct values in vals.

scaled_nxp_c0_vec=nxp_c0_vec./DELTA_vec; 

Tchebp=Tcheb_fun(scaled_nxp_c0_vec,[]);
Xp_vecT=sptensor(ones(model.nparams,1));
Xp_vecT.vals=Tchebp';

% Note that X is a function of DELTA^-1*(x-c0). When we differentiate X
% w.r.t to x we need to use the chain rule.
Xpx_vecT=Tcheb_d1(scaled_nxp_c0_vec',[],model.Tcheb_ind); % direct Jacobian

invDELTA=sptensor(1:n_x,(1:n_x)',1./DELTA_vec(:,1)',n_x,n_x); % internal Jacobian
[Xpx_vecT,model.ind{indi}]=contraction1(Xpx_vecT,invDELTA,model.ind{indi},n_ind,maxload,'vec'); % chain rule
indi=indi+1;

gp_vec=g_coeffs*Xp_vecT.vals';

stochg_coeffsT=sptensor(stochg_coeffs);



[stochgpx_vecT,model.ind{indi}]=contraction1(stochg_coeffsT,Xpx_vecT,model.ind{indi},n_ind,maxload,'vec');
indi=indi+1;

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

stochPIz_vec=stochPI_d1(nz_vec',params,model.stochPI_ind_u);
prePIz_full.vals=repmat(prePIz_full.vals,n_s,1);

i1=tfind(prePIz_full);
i2=tfind(stochPIz_vec);
PIz_vec=sptensor([i1;i2],[prePIz_full.cols;stochPIz_vec.cols],[prePIz_full.vals,stochPIz_vec.vals],n_z,[n_z]);
stochtilfz_vec=stochtilf_tilf_d1([nv_vec(stochfvars,:);nu_vec(model.stochfuvars,:)]',params,model.stochtilf_ind_u);

EstochR0=reshape(stochf_vec,n_stochfrows*grid_points,n_s)*P;
EstochR0=reshape(EstochR0,n_stochfrows,grid_points);

h_thetaT_vec=h_theta;
h_thetaT_vec.vals=repmat(h_thetaT_vec.vals,n_s,1);
chain0_theta_stochgpT=contraction1(stochgpx_vecT,h_thetaT_vec,model.ind{indi},n_ind,maxload,'vec');
indi=indi+1;

% Jacobian

[stochPIz_vec,model.ind{indi}]=extract(PIz_vec,model.stochfzvars,model.stochfzvars,0,model.ind{indi});
indi=indi+1;
[stochPIv_vec,model.ind{indi}]=extract(PIz_vec,model.stochfzvars,model.stochfvars,0,model.ind{indi});
indi=indi+1;  

for i=2:model.stoch_n
    [stochPIv_vec,model.ind{indi}]=chain1_tensor(stochPIz_vec,stochPIv_vec,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
end
[stochfv_vec,model.ind{indi}]=chain1_tensor(stochtilfz_vec,stochPIv_vec,model.ind{indi},n_ind,maxload,'vec');
indi=indi+1;

    
temp=reshape(1:n_theta,n_f,nparams);
stochtheta=temp(stochy,:);
stochtheta=stochtheta(:);
n_stochtheta=length(stochtheta);

IstochyT=spteye(n_stochy);
[stochgp_stochtheta_vecT,model.ind{indi}]=tkron(ttranspose(Xp_vecT),IstochyT,model.ind{indi},n_ind,maxload,'vec');
indi=indi+1;


nnzgp_theta_vecT=unfold(changecols(changerows(stochgp_stochtheta_vecT,stochy,n_y),stochy,n_f,1)); %n_y,n_theta with only nnz values
nnzchain0_theta_gpT=changerows(chain0_theta_stochgpT,stochy,n_y);

v_theta_lower_terms_vecT=v_theta_lower_terms;
v_theta_lower_terms_vecT.vals=repmat(v_theta_lower_terms_vecT.vals,n_s,1);

v_theta_vecT=vconcat(tplus(nnzchain0_theta_gpT,nnzgp_theta_vecT,maxload),v_theta_lower_terms_vecT);
clear nnzchain0_theta_gpT nnzgp_theta_vecT
stochv_theta_vecT=takerows(v_theta_vecT,stochfvars);


clear v_theta_vecT

[ stochJ0,model.ind{indi} ] = chain0_theta_tensor( stochfv_vec,stochv_theta_vecT,model.ind{indi},n_ind,maxload,'vec' );
indi=indi+1;

stochJ0.vals=(reshape(stochJ0.vals',(stochJ0.ptr(end)-1)*grid_points,n_s)*P)';

EstochJ0=stochJ0;
EstochJ0.vals=EstochJ0.vals(1,:);
[i,j]=tfind(EstochJ0);
i=double(i);
j=double(j);

EstochJ0=sparse(repmat(i,grid_points,1)+kron((0:grid_points-1)',n_stochfrows*ones(length(i),1)),...
    repmat(j,grid_points,1),...
    stochJ0.vals,...
    n_stochfrows*grid_points,n_theta);

R=zeros(n_f,grid_points);
R(model.prefrows,:)=preR0;
R(model.stochfrows,:)=EstochR0;

J=sparse(n_f,grid_points*n_theta);
J(model.prefrows,:)=reshape(preJ0,n_prefrows,grid_points*n_theta);
J(model.stochfrows,:)=reshape(EstochJ0,n_stochfrows,grid_points*n_theta);

R=R(:);
J=reshape(J,n_f*grid_points,n_theta);


