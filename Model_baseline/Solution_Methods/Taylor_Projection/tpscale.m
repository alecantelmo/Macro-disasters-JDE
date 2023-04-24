function [ scaledT,scaledJ,Wcols,Wrows,model ] = tpscale( weighted_coeffs,x0,model,params,eta,c0,nep,P,Wcols,Wrows,varargin )
%This function computes the scaled system and its Jacobian. If weights are
%empty, the function computes and returns the weights. Otherwise, the
%supplied weights are used.
%
% © Copyright, Oren Levintal, June 13, 2016.

order=model.order(1);

if isempty(Wcols) && isempty(Wrows) 
    % if weights are not supplied assume no weights.
    coeffs=weighted_coeffs;
elseif ~isempty(Wcols) && ~isempty(Wrows)
    % if all weights are supplied, unscale weighted_coeffs
    coeffs=weighted_coeffs./Wcols;
else
    error('missing weights');
end

if isempty(varargin) % precomputed variables are not available
    [T,J,model]=tp(coeffs,x0,model,params,eta,c0,nep,P);
else
    g_theta=varargin{1};
    h_theta=varargin{2};
    gx_theta=varargin{3};
    hx_theta=varargin{4};
    gxxc_theta=varargin{5};
    hxxc_theta=varargin{6};
    gxxxc_theta=varargin{7};
    hxxxc_theta=varargin{8};
    [T,J,model]=tp(coeffs,x0,model,params,eta,c0,nep,P,...
    g_theta,h_theta,gx_theta,hx_theta,gxxc_theta,hxxc_theta,gxxxc_theta,hxxxc_theta);
end

if isempty(Wcols)
    Wcols=max(abs(J))';
    Wrows=max(abs([J*spdiags(1./Wcols,0,model.n_theta,model.n_theta)]'))';
end

scaledT=T./Wrows;

scaledJ=spdiags(1./Wrows,0,model.n_theta,model.n_theta)*J*spdiags(1./Wcols,0,model.n_theta,model.n_theta);


end

