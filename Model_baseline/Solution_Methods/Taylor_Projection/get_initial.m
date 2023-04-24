function [ stoch_pert,nonstoch_pert ] = get_initial(model,order,derivs,nyss,nxss )
%This function transforms a pertrubation solution to an initial guess for Taylor
%Projection. Two types of guess are available: stoch_pert includes
%correction for volatility. nonstoch_pert does not include correction for
%volatility.
%
% © Copyright, Oren Levintal, June 13, 2016.

order=order(1);

n_x=model.n_x;
n_x1=model.n_x1;
n_y=model.n_y;

if order==1
    [GH0,GH1]=shiftpoly([nxss;1],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx;derivs.hx(1:n_x1,:)]);
    stoch_pert=[GH0,GH1];
    if nargout==2
        [GH0,GH1]=shiftpoly([nxss;0],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx;derivs.hx(1:n_x1,:)]);
        nonstoch_pert=[GH0,GH1];
    end
elseif order==2
    derivs.gxx=reshape(derivs.gxx,n_y,[]);
    derivs.hxx=reshape(derivs.hxx(1:n_x1,:),n_x1,[]);
    [GH0,GH1,GH2]=shiftpoly([nxss;1],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx;derivs.hx(1:n_x1,:)],...
        [derivs.gxx;derivs.hxx(1:n_x1,:)]/2);
    [U2,~]=create_UW(n_x,2);
    stoch_pert=[GH0,GH1,GH2*U2];
    if nargout==2
        [GH0,GH1,GH2]=shiftpoly([nxss;0],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx;derivs.hx(1:n_x1,:)],...
            [derivs.gxx;derivs.hxx(1:n_x1,:)]/2);
        nonstoch_pert=[GH0,GH1,GH2*U2];
    end
elseif order==3
    derivs.gxx=reshape(derivs.gxx,n_y,[]);
    derivs.hxx=reshape(derivs.hxx(1:n_x1,:),n_x1,[]);
    derivs.gxxx=reshape(derivs.gxxx,n_y,[]);
    derivs.hxxx=reshape(derivs.hxxx(1:n_x1,:),n_x1,[]);
    [GH0,GH1,GH2,GH3]=shiftpoly([nxss;1],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx;derivs.hx(1:n_x1,:)],...
        [derivs.gxx;derivs.hxx(1:n_x1,:)]/2,[derivs.gxxx;derivs.hxxx(1:n_x1,:)]/factorial(3));
    [U2,~]=create_UW(n_x,2);
    [U3,~]=create_UW(n_x,3);
    stoch_pert=[GH0,GH1,GH2*U2,GH3*U3];
    if nargout==2
        [GH0,GH1,GH2,GH3]=shiftpoly([nxss;0],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx;derivs.hx(1:n_x1,:)],...
            [derivs.gxx;derivs.hxx(1:n_x1,:)]/2,[derivs.gxxx;derivs.hxxx(1:n_x1,:)]/factorial(3));
        nonstoch_pert=[GH0,GH1,GH2*U2,GH3*U3];
    end
elseif order==4
    derivs.gxx=reshape(derivs.gxx,n_y,[]);
    derivs.hxx=reshape(derivs.hxx(1:n_x1,:),n_x1,[]);
    derivs.gxxx=reshape(derivs.gxxx,n_y,[]);
    derivs.hxxx=reshape(derivs.hxxx(1:n_x1,:),n_x1,[]);
    derivs.gxxxx=reshape(derivs.gxxxx,n_y,[]);
    derivs.hxxxx=reshape(derivs.hxxxx(1:n_x1,:),n_x1,[]);
    [GH0,GH1,GH2,GH3,GH4]=shiftpoly([nxss;1],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx;derivs.hx(1:n_x1,:)],...
        [derivs.gxx;derivs.hxx(1:n_x1,:)]/2,[derivs.gxxx;derivs.hxxx(1:n_x1,:)]/factorial(3),...
        [derivs.gxxxx;derivs.hxxxx(1:n_x1,:)]/factorial(4));
    [U2,~]=create_UW(n_x,2);
    [U3,~]=create_UW(n_x,3);
    [U4,~]=create_UW(n_x,4);
    stoch_pert=[GH0,GH1,GH2*U2,GH3*U3,GH4*U4];
    if nargout==2
        [GH0,GH1,GH2,GH3,GH4]=shiftpoly([nxss;0],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx;derivs.hx(1:n_x1,:)],...
            [derivs.gxx;derivs.hxx(1:n_x1,:)]/2,[derivs.gxxx;derivs.hxxx(1:n_x1,:)]/factorial(3),...
            [derivs.gxxxx;derivs.hxxxx(1:n_x1,:)]/factorial(4));
        nonstoch_pert=[GH0,GH1,GH2*U2,GH3*U3,GH4*U4];
    end
elseif order==5
    [GH0,GH1,GH2,GH3,GH4,GH5]=shiftpoly([nxss;1],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx;derivs.hx(1:n_x1,:)],...
        [derivs.gxx;derivs.hxx(1:n_x1,:)]/2,[derivs.gxxx;derivs.hxxx(1:n_x1,:)]/factorial(3),...
        [derivs.gxxxx;derivs.hxxxx(1:n_x1,:)]/factorial(4),[derivs.gxxxxx;derivs.hxxxxx(1:n_x1,:)]/factorial(5));
    [U2,~]=create_UW(n_x,2);
    [U3,~]=create_UW(n_x,3);
    [U4,~]=create_UW(n_x,4);
    [U5,~]=create_UW(n_x,5);
    stoch_pert=[GH0,GH1,GH2*U2,GH3*U3,GH4*U4,GH5*U5];
    if nargout==2
        [GH0,GH1,GH2,GH3,GH4,GH5]=shiftpoly([nxss;0],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx;derivs.hx(1:n_x1,:)],...
            [derivs.gxx;derivs.hxx(1:n_x1,:)]/2,[derivs.gxxx;derivs.hxxx(1:n_x1,:)]/factorial(3),...
            [derivs.gxxxx;derivs.hxxxx(1:n_x1,:)]/factorial(4),[derivs.gxxxxx;derivs.hxxxxx(1:n_x1,:)]/factorial(5));
        nonstoch_pert=[GH0,GH1,GH2*U2,GH3*U3,GH4*U4,GH5*U5];
    end
end
stoch_pert=stoch_pert(:);

if nargout==2
    nonstoch_pert=nonstoch_pert(:);
end

end

