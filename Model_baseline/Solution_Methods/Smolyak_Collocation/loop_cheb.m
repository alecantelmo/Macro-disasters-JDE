function [R,J,model]=loop_cheb(coeffs,model,params,x_grid,c0,DELTA,nep,P,eta,maxload)

% Computes residual function and Jacobian by looping over the grid points.
% Each iteration computes simulatenously multiple grid points by vectorized
% operations. Maximum number of grid points in each iteration is controlled
% by maxload. In addition, model.maxload also controls the size of
% vectorized expressions inside the MEX functions. Both, maxload and
% model.maxload affect memory (and to some extent also speed).
%
% © Copyright, Jesus Fernandez-Villaverde and Oren Levintal, June 13, 2016.

%%%%%%%%%%%%
grid_points=size(x_grid,2);
n_s=length(P);
n_e=size(eta,2);
n_f=model.n_f;
n_theta=model.n_theta;
n_x=model.n_x;
n_x1=model.n_x1;
n_x2=model.n_x2;
n_y=model.n_y;

R=zeros(model.n_f*grid_points,1);
J=zeros(model.n_f*grid_points,model.n_theta); % i'm using full matrix, becuase it's faster to compute the Newton step. The Jacobian is not sparse enough to justify sparse operations. nnz(J)/numel(J) is .27 for the largest model.
grid_vec=1:maxload:grid_points;
grid_vec=[grid_vec,grid_points(end)+1];
loc=0;

for i=1:length(grid_vec)-1
%     disp([num2str(i) '/' num2str(length(grid_vec)-1)])
    x_grid_i=x_grid(:,grid_vec(i):grid_vec(i+1)-1);

    Phi=Phi_fun(x_grid_i,params);
    scaled_grid_i=(x_grid_i-repmat(c0,1,size(x_grid_i,2)))./repmat(DELTA,1,size(x_grid_i,2)); 
    grid_points_i=grid_vec(i+1)-grid_vec(i);
    c0_vec=repmat(c0,1,grid_points_i*n_s);
    DELTA_vec=repmat(DELTA,1,grid_points_i*n_s);

    Tcheb=Tcheb_fun(scaled_grid_i,[]);
    X=sptensor(ones(model.nparams,1));
    X.vals=Tcheb';
    if ~isfield(model,'ind')
    model.ind=cell(200,1);
    end
    indi=1;

    if ~isfield(model,'n_ind')
        n_ind=2;
    else
        n_ind=model.n_ind;
    end
    IfT=spteye(n_f);
    [gh_theta,model.ind{indi}]=tkron(ttranspose(X),IfT,model.ind{indi},n_ind,maxload,'vec');
    indi=indi+1;
    gh_theta=unfold(gh_theta);
    g_theta=takerows(gh_theta,1:n_y);
    h_theta=vconcat(takerows(gh_theta,n_y+1:n_f),sptensor(n_x2,n_theta,grid_points_i));

    [tempR,tempJ,model]=fun_collocation(coeffs,model,params,x_grid_i,c0_vec,DELTA_vec,nep,P,Tcheb,...
        g_theta,h_theta,Phi,eta);
    
    R(loc+1:loc+length(tempR))=tempR;
    J(loc+1:loc+length(tempR),:)=tempJ;
    loc=loc+length(tempR);
end


