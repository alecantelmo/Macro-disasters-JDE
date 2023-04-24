function [ncoeffs,R,J,iter]=tp_newton(coeffs,x0,model,params,eta,c0,nep,P,order,tolX,tolF,maxiter,varargin)
% The function solve Taylor projection by Newton method.
%
% © Copyright, Oren Levintal, June 13, 2016.

order=order(1);

% make sure all inputs are full
coeffs=full(coeffs);
x0=full(x0);
params=full(params);
eta=full(eta);
c0=full(c0);
P=full(P);

n_x=model.n_x;
n_x2=model.n_x2;
n_y=model.n_y;
n_f=model.n_f;
n_theta=model.n_theta;
U=model.U;
W=model.W;
unique2=model.unique2;
unique3=model.unique3;

% shift the center of the Polynomials to x0

c0_old=c0;

if ~isequal(c0_old,x0)
    coeffs=reshape(coeffs,n_f,model.n_b);
    GH0=coeffs(:,1);
    if order==1
        GH1=coeffs(:,2:1+n_x);
        [ GH0,GH1 ] = shiftpoly( x0,c0,[],GH0,GH1 );
        coeffs=[GH0,GH1];
        coeffs=coeffs(:);
        c0=x0;
    end
    if order==2
        GH1=coeffs(:,2:1+n_x);
        GH2=coeffs(:,2+n_x:1+n_x+unique2)*W{2};
        [ GH0,GH1,GH2 ] = shiftpoly( x0,c0,[],GH0,GH1,GH2 );
        coeffs=[GH0,GH1,GH2*U{2}];
        coeffs=coeffs(:);
        c0=x0;
    end
    if order==3
        GH1=coeffs(:,2:1+n_x);
        GH2=coeffs(:,2+n_x:1+n_x+unique2)*W{2};
        GH3=coeffs(:,2+n_x+unique2:1+n_x+unique2+unique3)*W{3};
        [ GH0,GH1,GH2,GH3 ] = shiftpoly( x0,c0,[],GH0,GH1,GH2,GH3 );
        coeffs=[GH0,GH1,GH2*U{2},GH3*U{3}];
        coeffs=coeffs(:);
        c0=x0;
    end
end

% precompute expressions that are independent of coeffs
[ g_theta,h_theta,gx_theta,hx_theta,gxxc_theta,hxxc_theta,gxxxc_theta,hxxxc_theta,model ] = precompute(x0,c0,model,order );


iter=1;

% start without weights.

Wcols=[]; Wrows=[]; coeffs=full(coeffs);

[R,J,Wcols,Wrows,model]=tpscale(coeffs,x0,model,params,eta,c0,nep,P,Wcols,Wrows,...
      g_theta,h_theta,gx_theta,hx_theta,gxxc_theta,hxxc_theta,gxxxc_theta,hxxxc_theta);  
normR=norm(full(R));

ncoeffs=coeffs-full(J\R)./Wcols;
d=ncoeffs-coeffs;
normd=norm(full(d));
optimality=norm(full(R'*J),inf);

disp(['iter    norm(R)   norm(step)    optimality measure'])

disp([' ' num2str(iter) '    ' num2str(normR,'%3.2e') '   ' num2str(normd,'%3.2e') '          ' num2str(optimality,'%3.2e')] )

% continue to convergence
while normd > tolX && normR>tolF && iter<maxiter 
    iter=iter+1;
    coeffs=ncoeffs;
    weighted_coeffs=full(Wcols.*coeffs);
    [R,J]=tpscale(weighted_coeffs,x0,model,params,eta,c0,nep,P,Wcols,Wrows,...
        g_theta,h_theta,gx_theta,hx_theta,gxxc_theta,hxxc_theta,gxxxc_theta,hxxxc_theta);  

    normR=norm(full(R));
    
    ncoeffs=coeffs-full(J\R)./Wcols;


    
    d=ncoeffs-coeffs;
    optimality=norm(full(R'*J),inf);
    normd=norm(full(d));

    disp([' ' num2str(iter) '    ' num2str(normR,'%3.2e') '   ' num2str(normd,'%3.2e') '          ' num2str(optimality,'%3.2e')] )
end

R=R.*Wrows;

% shift the center back to the original c0
if ~isequal(c0_old,x0)
    coeffs=reshape(coeffs,n_f,model.n_b);
    GH0=coeffs(:,1);
    if order==1
        GH1=coeffs(:,2:1+n_x);
        [ GH0,GH1 ] = shiftpoly( c0_old,c0,[],GH0,GH1 );
        coeffs=[GH0,GH1];
        coeffs=coeffs(:);
    end
    if order==2
        GH1=coeffs(:,2:1+n_x);
        GH2=coeffs(:,2+n_x:1+n_x+unique2)*W{2};
        [ GH0,GH1,GH2 ] = shiftpoly( c0_old,c0,[],GH0,GH1,GH2 );
        coeffs=[GH0,GH1,GH2*U{2}];
        coeffs=coeffs(:);
    end
    if order==3
        GH1=coeffs(:,2:1+n_x);
        GH2=coeffs(:,2+n_x:1+n_x+unique2)*W{2};
        GH3=coeffs(:,2+n_x+unique2:1+n_x+unique2+unique3)*W{3};
        [ GH0,GH1,GH2,GH3 ] = shiftpoly( c0_old,c0,[],GH0,GH1,GH2,GH3 );
        coeffs=[GH0,GH1,GH2*U{2},GH3*U{3}];
        coeffs=coeffs(:);
    end
end
ncoeffs=full(coeffs);
clearvars -global Wcols Wrows