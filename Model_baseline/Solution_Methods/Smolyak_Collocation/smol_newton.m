function [ ncoeffs,R,J,iter ] = smol_newton( coeffs,model,params,x_grid,cheb_c0,DELTA,nep,P,eta,maxload,tolX,tolF,maxiter )
% Newton method to solve the Smolyak coefficients.
%
% © Copyright, Jesus Fernandez-Villaverde and Oren Levintal, June 13, 2016.

iter=0;
ncoeffs=coeffs;
R=[];
J=[];
stopxrule=Inf;
stopfrule=Inf;

disp(['iter    norm(R)   norm(step)    optimality measure'])

while stopxrule > tolX && iter<maxiter && stopfrule>tolF 
    iter=iter+1;
    [R,J,model]=loop_cheb(coeffs,model,params,x_grid,cheb_c0,DELTA,nep,P,eta,maxload);
    
    ncoeffs=coeffs-full(J\R);
    d=ncoeffs-coeffs;
    step=norm(full(d));
    normR=norm(full(R));
    coeffs=ncoeffs;

    optimality=norm(R'*J,inf);

    disp([' ' num2str(iter) '    ' num2str(normR,'%3.2e') '   ' num2str(step,'%3.2e') '          ' num2str(optimality,'%3.2e')])
    stopxrule=step;
    stopfrule=normR;
end

end

