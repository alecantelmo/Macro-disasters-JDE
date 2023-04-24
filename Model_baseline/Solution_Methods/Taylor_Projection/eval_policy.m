function [g,h1,nPhi]=eval_policy(coeffs,x0_vec,model,params,c0)
%
% © Copyright, Oren Levintal, June 13, 2016.
     
order=model.order(1);

params=params(:);


n_f=model.n_f;
n_y=model.n_y;

n_b=model.n_b;

% current state
nx=x0_vec;
n_grid=size(x0_vec,2);

nx_c0_vec=nx-repmat(c0,1,n_grid);
nx_c0=sptensor(nx_c0_vec(:,1)); % convert to 1D tensor
nx_c0.vals=nx_c0_vec'; % assign the correct values in vals.

[X_vecT]=create_X_tensor_no_derivs(order,nx_c0,...
    model.M2,model.M3,[],model.n_ind,n_grid,'vec');

gh_coeffs=reshape(coeffs,n_f,n_b);

gh1=gh_coeffs*X_vecT.vals';
g=gh1(1:n_y,:);
h1=gh1(n_y+1:end,:);

nPhi=Phi_fun(nx,params); % expected value of exogenous state vars. shocks are added later

