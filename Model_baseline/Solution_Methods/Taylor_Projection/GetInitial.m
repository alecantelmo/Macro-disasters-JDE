function [ pert_coeffs,power_coeffs ] = GetInitial( order,derivs,choose,nyss,nxss,n_x,n_x1 )
%This function transforms a pertrubation solution to an initial guess for Taylor
%Projection. Two types of guess are available: pert_coeffs includes
%correction for volatility. power_coeffs does not include correction for
%volatility.
%
% © Copyright, Oren Levintal, June 13, 2016.


if order==1
    [GH0,GH1]=shiftpoly([nxss;1],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx(choose,:);derivs.hx(1:n_x1,:)]);
    pert_coeffs=[GH0,GH1];
    if nargout==2
        [GH0,GH1]=shiftpoly([nxss;0],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx(choose,:);derivs.hx(1:n_x1,:)]);
        power_coeffs=[GH0,GH1];
    end
elseif order==2
    [GH0,GH1,GH2]=shiftpoly([nxss;1],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx(choose,:);derivs.hx(1:n_x1,:)],...
        [derivs.gxx(choose,:);derivs.hxx(1:n_x1,:)]/2);
    [U2,~]=create_UW(n_x,2);
    pert_coeffs=[GH0,GH1,GH2*U2];
    if nargout==2
        [GH0,GH1,GH2]=shiftpoly([nxss;0],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx(choose,:);derivs.hx(1:n_x1,:)],...
            [derivs.gxx(choose,:);derivs.hxx(1:n_x1,:)]/2);
        power_coeffs=[GH0,GH1,GH2*U2];
    end
elseif order==3
    [GH0,GH1,GH2,GH3]=shiftpoly([nxss;1],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx(choose,:);derivs.hx(1:n_x1,:)],...
        [derivs.gxx(choose,:);derivs.hxx(1:n_x1,:)]/2,[derivs.gxxx(choose,:);derivs.hxxx(1:n_x1,:)]/factorial(3));
    [U2,~]=create_UW(n_x,2);
    [U3,~]=create_UW(n_x,3);
    pert_coeffs=[GH0,GH1,GH2*U2,GH3*U3];
    if nargout==2
        [GH0,GH1,GH2,GH3]=shiftpoly([nxss;0],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx(choose,:);derivs.hx(1:n_x1,:)],...
            [derivs.gxx(choose,:);derivs.hxx(1:n_x1,:)]/2,[derivs.gxxx(choose,:);derivs.hxxx(1:n_x1,:)]/factorial(3));
        power_coeffs=[GH0,GH1,GH2*U2,GH3*U3];
    end
elseif order==4
    [GH0,GH1,GH2,GH3,GH4]=shiftpoly([nxss;1],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx(choose,:);derivs.hx(1:n_x1,:)],...
        [derivs.gxx(choose,:);derivs.hxx(1:n_x1,:)]/2,[derivs.gxxx(choose,:);derivs.hxxx(1:n_x1,:)]/factorial(3),...
        [derivs.gxxxx(choose,:);derivs.hxxxx(1:n_x1,:)]/factorial(4));
    [U2,~]=create_UW(n_x,2);
    [U3,~]=create_UW(n_x,3);
    [U4,~]=create_UW(n_x,4);
    pert_coeffs=[GH0,GH1,GH2*U2,GH3*U3,GH4*U4];
    if nargout==2
        [GH0,GH1,GH2,GH3,GH4]=shiftpoly([nxss;0],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx(choose,:);derivs.hx(1:n_x1,:)],...
            [derivs.gxx(choose,:);derivs.hxx(1:n_x1,:)]/2,[derivs.gxxx(choose,:);derivs.hxxx(1:n_x1,:)]/factorial(3),...
            [derivs.gxxxx(choose,:);derivs.hxxxx(1:n_x1,:)]/factorial(4));
        power_coeffs=[GH0,GH1,GH2*U2,GH3*U3,GH4*U4];
    end
elseif order==5
    [GH0,GH1,GH2,GH3,GH4,GH5]=shiftpoly([nxss;1],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx(choose,:);derivs.hx(1:n_x1,:)],...
        [derivs.gxx(choose,:);derivs.hxx(1:n_x1,:)]/2,[derivs.gxxx(choose,:);derivs.hxxx(1:n_x1,:)]/factorial(3),...
        [derivs.gxxxx(choose,:);derivs.hxxxx(1:n_x1,:)]/factorial(4),[derivs.gxxxxx(choose,:);derivs.hxxxxx(1:n_x1,:)]/factorial(5));
    [U2,~]=create_UW(n_x,2);
    [U3,~]=create_UW(n_x,3);
    [U4,~]=create_UW(n_x,4);
    [U5,~]=create_UW(n_x,5);
    pert_coeffs=[GH0,GH1,GH2*U2,GH3*U3,GH4*U4,GH5*U5];
    if nargout==2
        [GH0,GH1,GH2,GH3,GH4,GH5]=shiftpoly([nxss;0],[nxss;0],1:n_x,[nyss;nxss(1:n_x1,:)],[derivs.gx(choose,:);derivs.hx(1:n_x1,:)],...
            [derivs.gxx(choose,:);derivs.hxx(1:n_x1,:)]/2,[derivs.gxxx(choose,:);derivs.hxxx(1:n_x1,:)]/factorial(3),...
            [derivs.gxxxx(choose,:);derivs.hxxxx(1:n_x1,:)]/factorial(4),[derivs.gxxxxx(choose,:);derivs.hxxxxx(1:n_x1,:)]/factorial(5));
        power_coeffs=[GH0,GH1,GH2*U2,GH3*U3,GH4*U4,GH5*U5];
    end
end
pert_coeffs=pert_coeffs(:);

if nargout==2
    power_coeffs=power_coeffs(:);
end

end

