%
% © Copyright, Oren Levintal, June 13, 2016.

OMEGA_x=[];
if approx>=3
    tempOMEGA = create_OMEGA( n_x,3 );
    OMEGA_x.OMEGA1=tempOMEGA.OMEGA1;
end

if approx>=4
    tempOMEGA = create_OMEGA( n_x,4 );
    OMEGA_x.OMEGA2=tempOMEGA.OMEGA2;
    OMEGA_x.OMEGA3=tempOMEGA.OMEGA3;
    OMEGA_x.OMEGA4=tempOMEGA.OMEGA4;
end
if approx>=5
    tempOMEGA = create_OMEGA( n_x,5 );
    OMEGA_x.OMEGA5=tempOMEGA.OMEGA5;
    OMEGA_x.OMEGA6=tempOMEGA.OMEGA6;
    OMEGA_x.OMEGA7=tempOMEGA.OMEGA7;
    OMEGA_x.OMEGA8=tempOMEGA.OMEGA8;
    OMEGA_x.OMEGA9=tempOMEGA.OMEGA9;
end
