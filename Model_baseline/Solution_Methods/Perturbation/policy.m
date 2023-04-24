function [ g,h ] = policy( nx,nyss,nxss,derivs,approx,varargin )
%[g,h]=policy(nx,nxss,derivs,approx) evaluates the policy functions g and h
% at state nx, given the steady state values nyss and nxss, the structure 
% derivs and the approximation order approx.
%
% © Copyright, Oren Levintal, June 13, 2016.

[n_y,n_x]=size(derivs.gx);


xhat=[nx(:)-nxss(:);1];

pf=[derivs.gx;derivs.hx]*xhat;

n_f=length(pf);

if isempty(varargin)

    if approx>=2
        xhat2=kron(xhat,xhat);
        pf=pf+[derivs.gxx;derivs.hxx]*xhat2/2;
    end
    if approx>=3
        xhat3=kron(xhat2,xhat);
        pf=pf+[derivs.gxxx;derivs.hxxx]*xhat3/6;
    end
    if approx>=4
        xhat4=kron(xhat3,xhat);
        pf=pf+[derivs.gxxxx;derivs.hxxxx]*xhat4/24;
    end
    if approx>=5
        xhat5=kron(xhat4,xhat);
        pf=pf+[derivs.gxxxxx;derivs.hxxxxx]*xhat5/120;
    end
else
    derivsc=varargin{1};
    if approx>=2
        pf=pf+innerkron(n_f,n_x,[derivsc.gxx;derivsc.hxx],xhat,xhat)/2;
    end
    if approx>=3
        pf=pf+innerkron(n_f,n_x,[derivsc.gxxx;derivsc.hxxx],xhat,xhat,xhat)/6;
    end
    if approx>=4
        pf=pf+innerkron(n_f,n_x,[derivsc.gxxxx;derivsc.hxxxx],xhat,xhat,xhat,xhat)/24;
    end
    if approx>=5
        pf=pf+innerkron(n_f,n_x,[derivsc.gxxxxx;derivsc.hxxxxx],xhat,xhat,xhat,xhat,xhat)/120;
    end
end
R=pf+[nyss;nxss];

g=R(1:n_y);
h=R(n_y+1:end);

end

