function [ncoeffs,model,T,J,iter]=tpsolve(coeffs,x0,model,params,eta,c0,nep,P,tolX,tolF,maxiter,varargin)

% The function computes the polynomial coefficients by the Newton method or
% by fsolve/lsqnonlin.
%
% © Copyright, Oren Levintal, June 13, 2016.

display=1;
if isfield(model,'display')
    if strcmp(model.display,'off')
        display=0;
    end
end

order=model.order(1);

% make sure all inputs are full
coeffs=full(coeffs);
x0=full(x0);
params=full(params);
eta=full(eta);
c0=full(c0);
P=full(P);

n_x=model.n_x;

% shift the center of the Polynomials to x0

c0_old=c0;
c0_new=x0;
if ~isequal(c0_old,x0)
    [ coeffs ] = shift_center( coeffs,c0_old,c0_new,order,model );
    coeffs=full(coeffs);
    c0=c0_new;
end

% precompute expressions that are independent of coeffs
[ g_theta,h_theta,gx_theta,hx_theta,gxxc_theta,hxxc_theta,gxxxc_theta,hxxxc_theta,model ] = precompute(x0,c0,model,order );

iter=1;

% start without weights.

Wcols=[]; Wrows=[]; coeffs=full(coeffs);

[T,J,Wcols,Wrows,model]=tpscale(coeffs,x0,model,params,eta,c0,nep,P,Wcols,Wrows,...
      g_theta,h_theta,gx_theta,hx_theta,gxxc_theta,hxxc_theta,gxxxc_theta,hxxxc_theta);  
normT=norm(T);

bestT=normT;
best_weighted_coeffs=Wcols.*coeffs;

ncoeffs=coeffs-full(J\T)./Wcols;
d=ncoeffs-coeffs;
normd=norm(d);
optimality=norm(T'*J,inf);

if display==1
    disp(['iter    norm(T)   norm(step)    optimality measure'])

    disp([' ' num2str(iter) '    ' num2str(normT,'%3.2e') '   ' num2str(normd,'%3.2e') '          ' num2str(optimality,'%3.2e')] )
end

% continue to convergence
while normd > max(tolX,eps) && normT>max(tolF,eps) && iter<maxiter 
    iter=iter+1;
    coeffs=ncoeffs;
    weighted_coeffs=full(Wcols.*coeffs);
    [T,J]=tpscale(weighted_coeffs,x0,model,params,eta,c0,nep,P,Wcols,Wrows,...
        g_theta,h_theta,gx_theta,hx_theta,gxxc_theta,hxxc_theta,gxxxc_theta,hxxxc_theta);  

    normT=norm(T);
    if normT<bestT
        bestT=normT;
        best_weighted_coeffs=weighted_coeffs;
    end
    
    ncoeffs=coeffs-full(J\T)./Wcols;

    d=ncoeffs-coeffs;
    optimality=norm(T'*J,inf);
    normd=norm(d);
    if display==1
        disp([' ' num2str(iter) '    ' num2str(normT,'%3.2e') '   ' num2str(normd,'%3.2e') '          ' num2str(optimality,'%3.2e')] )
    end
end

% test if the last iteration improves best_weighted_coeffs, by computing accuracy in the
% neighborhood of x0.

xgrid=repmat(x0,1,50)+(floor(rand(n_x,50)*5)-2)*.05;

resids_best=eval_model(best_weighted_coeffs./Wcols,xgrid,model,params,eta,c0,nep,P);
resids=eval_model(ncoeffs,xgrid,model,params,eta,c0,nep,P);

if norm(resids)<norm(resids_best)
    best_weighted_coeffs=Wcols.*ncoeffs;
    coeffs=ncoeffs;
else
    coeffs=best_weighted_coeffs./Wcols;
end

% if Newton method did not converge with the maximum number of iterations,
% switch to fsolve (default) or lsqnonlin.
if logical(normd < max(tolX,eps) || normT<tolF)==0
    
    if isempty(varargin)
        if display==1
            OPTIONS = optimoptions('fsolve','Jacobian','on','TolX',tolX,'TolF',tolF,'display','iter-detailed','MaxIter',maxiter);
        else
            OPTIONS = optimoptions('fsolve','Jacobian','on','TolX',tolX,'TolF',tolF,'display','off','MaxIter',maxiter);
        end
    else
        OPTIONS=varargin{1};
        OPTIONS=optimoptions(OPTIONS,'Jacobian','on');
    end
    
    weighted_coeffs=best_weighted_coeffs;
    
    if isa(OPTIONS,'optim.options.Fsolve')
        fprintf(2,'Could not find an exact solution, switching to fsolve ...\n');
        [ncoeffs,T,exitflag,output,J]=fsolve(@(coeffs) tpscale(coeffs,x0,model,params,eta,c0,nep,P,Wcols,Wrows,...
          g_theta,h_theta,gx_theta,hx_theta,gxxc_theta,hxxc_theta,gxxxc_theta,hxxxc_theta),weighted_coeffs,OPTIONS);
    elseif isa(OPTIONS,'optim.options.Lsqnonlin')
        fprintf(2,'Could not find an exact solution, switching to lsqnonlin ...\n');
        [ncoeffs,~,T,exitflag,output,~,J]=lsqnonlin(@(coeffs) tpscale(coeffs,x0,model,params,eta,c0,nep,P,Wcols,Wrows,...
          g_theta,h_theta,gx_theta,hx_theta,gxxc_theta,hxxc_theta,gxxxc_theta,hxxxc_theta),weighted_coeffs,[],[],OPTIONS);
    else
        error('solver must be fsolve or lsqnonlin')
    end

    ncoeffs=ncoeffs./Wcols;
    coeffs=ncoeffs;
    
end

T=T.*Wrows;

% shift the center back to the original c0
if ~isequal(c0_old,x0)
    [ coeffs ] = shift_center( coeffs,c0_new,c0_old,order,model );
end
ncoeffs=full(coeffs);

