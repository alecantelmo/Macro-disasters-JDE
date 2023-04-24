function [ results,nu ] = fun_simulation( x0,shocks,coeffs,c0,N,n_y,n_x,n_f,model,params,eta,varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

T=size(shocks,2);
results=zeros(n_y+n_x,T);
results(n_y+1:end,1)=x0;
coeffs=reshape(coeffs,n_f,numel(coeffs)/n_f);

if ~isempty(varargin)
    polyspecs=varargin{1};
else
    polyspecs.type='power';
end

et=shocks;

if strcmp(polyspecs.type,'power')
    n_u=model.n_u;
    nu=zeros(n_u,T);
    for t=1:T
        x0=results(n_y+1:end,t);
        xpowers=[1;x0-c0];
        lastkron=x0-c0;
        for i=2:N
            lastkron=kron(lastkron,x0-c0);
            xpowers=[xpowers;model.W{i}*lastkron];
        end
        gh1=coeffs*xpowers;
        ny=gh1(1:n_y);
        Exp=[gh1(n_y+1:end);Phi_fun(x0,params)];
        results(1:n_y,t)=ny;
        if t<T
            results(n_y+1:end,t+1)=Exp+eta*et(:,t+1);
        end
        nv=[zeros(n_y,1);ny;Exp;x0];
        npreu=preu_fun(nv(model.preuvars),params); % all predetermined u
        nu(model.preurows,t)=npreu;
    end
elseif strcmp(polyspecs.type,'smol')
    addpath(polyspecs.Tcheb_fun_folder);
    n_u=model.n_u;
    nu=zeros(n_u,T);
    for t=1:T
        x0=results(n_y+1:end,t);
        Tcheb=Tcheb_fun((x0-c0)./polyspecs.DELTA,[]);
        gh1=coeffs*Tcheb;
        ny=gh1(1:n_y);
        Exp=[gh1(n_y+1:end);Phi_fun(x0,params)];
        results(1:n_y,t)=ny;
        if t<T
            results(n_y+1:end,t+1)=Exp+eta*et(:,t+1);
        end
        nv=[zeros(n_y,1);ny;Exp;x0];
        npreu=preu_fun(nv(model.preuvars),params); % all predetermined u
        nu(model.preurows,t)=npreu;
    end
    rmpath(polyspecs.Tcheb_fun_folder);
end

end

