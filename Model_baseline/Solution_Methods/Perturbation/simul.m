function [yt,xt]=simul(x0,shocks,nyss,nxss,eta,derivs,approx,pruning,varargin)
% [yt,xt]=simul(x0,shocks,nyss,nxss,eta,derivs,approx,pruning) simulates
% the model from the initial state x0. shocks is a matrix with n_e rows
% and T columns, where n_e is the number of shocks (corresponds to the
% columns of eta), and T is the length of the simulation. The function
% returns yt and xt for T+1 periods. The first period is the initial state
% and the rest T periods correspond to the shocks. pruning=0 is a simple
% simulation without pruning. pruning=1 is a pruned simulation. The
% pruning algorithm follows Andreasen, Fernandez-Villaverde and
% Rubio-Ramirez (2013) "The Pruned State-Space System for Non-Linear DSGE Models:
% Theory and Empirical Applications".
%
% © Copyright, Oren Levintal, June 13, 2016.

if ~isempty(varargin)
    model=varargin{1};
    if approx>=2
        tempmat=model.UW.U2*model.UW.W2*model.UW.W2'*model.UW.U2';
        derivsc.gxx=sparse(derivs.gxx)*tempmat;
        derivsc.hxx=sparse(derivs.hxx)*tempmat;
    end
    if approx>=3
        tempmat=model.UW.U3*model.UW.W3*model.UW.W3'*model.UW.U3';
        derivsc.gxxx=sparse(derivs.gxxx)*tempmat;
        derivsc.hxxx=sparse(derivs.hxxx)*tempmat;
    end
    if approx>=4
        tempmat=model.UW.U4*model.UW.W4*model.UW.W4'*model.UW.U4';
        derivsc.gxxxx=sparse(derivs.gxxxx)*tempmat;
        derivsc.hxxxx=sparse(derivs.hxxxx)*tempmat;
    end
    if approx>=5
        tempmat=model.UW.U5*model.UW.W5*model.UW.W5'*model.UW.U5';
        derivsc.gxxxxx=sparse(derivs.gxxxxx)*tempmat;
        derivsc.hxxxxx=sparse(derivs.hxxxxx)*tempmat;
    end
end

T=size(shocks,2);
n_y=size(derivs.gx,1);
n_x=size(derivs.hx,1);
n_e=size(shocks,1);
shocks=[zeros(n_e,1),shocks,zeros(n_e,1)];
if pruning==0
    yt=zeros(n_y,T+2);
    xt=zeros(n_x,T+2);
    xt(:,1)=x0;
    for t=1:T+1
        nx=xt(:,t);
        if isempty(varargin)
            [g,h]=policy( nx,nyss,nxss,derivs,approx );
        else
            [g,h]=policy( nx,nyss,nxss,derivs,approx,derivsc );
        end
        yt(:,t)=g;
        xt(:,t+1)=h+eta*shocks(:,t+1);
    end
    xt=xt(:,1:T+1);
    yt=yt(:,1:T+1);
elseif pruning==1
    xt_f=zeros(n_x+1,T+2);
    yt=zeros(n_y,T+2);
    if approx>=2
        xt_s=zeros(n_x+1,T+2);
    end
    if approx>=3
        xt_rd=zeros(n_x+1,T+2);
    end
    if approx>=4
        xt_4th=zeros(n_x+1,T+2);
    end
    if approx>=5
        xt_5th=zeros(n_x+1,T+2);
    end

    xt_f(1:end-1,1)=x0-nxss;
    xt_f(end,:)=1;

    for t=1:T+1
        x_f=xt_f(:,t);
        xt_f(1:end-1,t+1)=derivs.hx*x_f+eta*shocks(:,t+1);
        if approx>=2
            x_s=xt_s(:,t);
            x_f2=kron(x_f,x_f);
            xt_s(1:end-1,t+1)=derivs.hx*x_s+derivs.hxx*x_f2/2;
        end
        if approx>=3
            x_rd=xt_rd(:,t);
            x_f3=kron(x_f2,x_f);
            x_f_x_s=kron(x_f,x_s);
            xt_rd(1:end-1,t+1)=derivs.hx*x_rd+derivs.hxx*(2*x_f_x_s)/2+derivs.hxxx*x_f3/6;
        end
        if approx>=4
            x_4th=xt_4th(:,t);
            x_f4=kron(x_f3,x_f);
            x_f2_x_s=kron(x_f,x_f_x_s);
            x_s2=kron(x_s,x_s);
            x_f_x_rd=kron(x_f,x_rd);
            xt_4th(1:end-1,t+1)=derivs.hx*x_4th+derivs.hxx*(2*x_f_x_rd+x_s2)/2 ...
                    +derivs.hxxx*(3*x_f2_x_s)/6 ...
                    +derivs.hxxxx*x_f4/24;
        end
        if approx>=5
            x_5th=xt_5th(:,t);
            x_f5=kron(x_f4,x_f);
            x_f3_x_s=kron(x_f,x_f2_x_s);
            x_f_x_s2=kron(x_f,x_s2);
            x_f2_x_rd=kron(x_f,x_f_x_rd);
            x_s_x_rd=kron(x_s,x_rd);
            x_f_x_4th=kron(x_f,x_4th);
            xt_5th(1:end-1,t+1)=derivs.hx*x_5th+derivs.hxx*(2*x_f_x_4th+2*x_s_x_rd)/2 ...
                +derivs.hxxx*(3*x_f2_x_rd+3*x_f_x_s2)/6 ...
                +derivs.hxxxx*(4*x_f3_x_s)/24 ...
                +derivs.hxxxxx*x_f5/120;
        end
        if approx==1
            yt(:,t)=derivs.gx*(x_f);
        elseif approx==2
            yt(:,t)=derivs.gx*(x_f+x_s)+derivs.gxx*(x_f2)/2;
        elseif approx==3
            yt(:,t)=derivs.gx*(x_f+x_s+x_rd)+derivs.gxx*(x_f2+2*x_f_x_s)/2 ...
                +derivs.gxxx*(x_f3)/6;
        elseif approx==4
            yt(:,t)=derivs.gx*(x_f+x_s+x_rd+x_4th)+derivs.gxx*(x_f2+2*x_f_x_s+2*x_f_x_rd+x_s2)/2 ...
                +derivs.gxxx*(x_f3+3*x_f2_x_s)/6 ...
                +derivs.gxxxx*x_f4/24;
        elseif approx==5
            yt(:,t)=derivs.gx*(x_f+x_s+x_rd+x_4th+x_5th)+derivs.gxx*(x_f2+2*x_f_x_s+2*x_f_x_rd+2*x_f_x_4th+x_s2+2*x_s_x_rd)/2 ...
                +derivs.gxxx*(x_f3+3*x_f2_x_s+3*x_f2_x_rd+3*x_f_x_s2)/6 ...
                +derivs.gxxxx*(x_f4+4*x_f3_x_s)/24 ...
                +derivs.gxxxxx*x_f5/120;
        end
    end
    yt=yt(:,1:T+1);
    xt=xt_f(:,1:T+1);
    if approx>=2
        xt=xt+xt_s(:,1:T+1);
    end
    if approx>=3
        xt=xt+xt_rd(:,1:T+1); 
    end
    if approx>=4
        xt=xt+xt_4th(:,1:T+1);
    end
    if approx>=5
        xt=xt+xt_5th(:,1:T+1);
    end
    yt=yt+repmat(nyss,1,T+1);
    xt=xt(1:end-1,:)+repmat(nxss,1,T+1);
end

