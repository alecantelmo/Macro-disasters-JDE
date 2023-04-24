function y=solve_SS(x)
LOSS_iii=0;
LOSS=0;
RHOGRANT_iii=0;

parameters

%%% Variables to find with fsolve
tilyss=x(1);
lss=x(2);

%%% Known variables in steady state %%%
pibackss=PISS;
piss=PISS;
dss=MUD;
thetass=THETABAR;
loghatass=LAMBDA_A-(1-ALPHA)*dss*thetass;
loghatmuss=0; % assume no investment specific shock (or trend)
loghatzss=1/(1-ALPHA)*loghatass+ALPHA/(1-ALPHA)*loghatmuss;
hatzss=exp(loghatzss);
hatass=exp(loghatass);
hatmuss=exp(loghatmuss);
mss=piss/RSS;
BETA=mss*(hatzss^PSI);
tilqss=1;
tilrss=tilqss/(mss*exp(-dss*thetass))-tilqss*(1-DELTA);
taucss=TAUCSS;
pistarss=((1-THETA_P*(piss^CHI/piss)^(1-EPSILON))/(1-THETA_P))^(1/(1-EPSILON));
vpss=(1-THETA_P)*pistarss^(-EPSILON)/(1-THETA_P*(piss^CHI/piss)^(-EPSILON));
RSTARss=RSS;

%%% Steady state of variables in levels
tilg2ss=(pistarss*tilyss)/( 1-THETA_P*mss*((piss^CHI/piss)^(1-EPSILON)) *hatzss );
tilg1ss=(EPSILON-1)*tilg2ss/EPSILON;
mcss=tilg1ss*(1-THETA_P*mss*(piss^CHI/piss)^(-EPSILON)*hatzss)/tilyss;
tilxgss=XGY*tilyss;
GBAR=GY*tilyss;
Bbackss=BY*4*tilyss;
grantss=GRANTY*tilyss;
tilkgss=(tilxgss)/(hatzss*hatmuss*exp(dss*thetass)-(1-DELTAG));
tilkgstarss=tilkgss*hatzss*hatmuss*exp(dss*thetass);
tilkss=(tilrss/(ALPHA*(1-ALPHAG)*mcss*(tilkgss^(ALPHA*ALPHAG))*(lss^(1-ALPHA))))^(1/(ALPHA*(1-ALPHAG)-1));
tilkstarss=tilkss*hatzss*hatmuss*exp(dss*thetass);
tilxss=tilkstarss-(1-DELTA)*tilkss;
tilnxss=NXY*tilyss;
REMIT=-Bbackss*(1-RSTARss/(hatzss*piss))+tilnxss-grantss;
tilcss=tilyss-tilxss-GBAR-tilxgss+tilnxss;
LUMPSUM=-(1-RSTARss/(hatzss*piss))*Bbackss + GBAR +tilxgss - taucss*tilcss- grantss;
tiluss=tilcss*(1-lss)^NU;
ucss=(1-lss)^NU;
ulss=-NU*tilcss*(1-lss)^(NU-1);
tillambdass=(1-PSI)*(tiluss^(-PSI))*ucss/(1+taucss);
tilwss=-(1-PSI)*(tiluss^(-PSI))*ulss/tillambdass;
logbyss=log(BY);
tiluss_2=(1-LOSS)*tilcss*(1-lss)^NU;

%%%these are the equations for which fsolve will find the solution
y=[ tilyss - ((hatass/hatzss)*( (tilkgstarss*exp(-dss*thetass))^ALPHAG * (tilkstarss*exp(-dss*thetass))^(1-ALPHAG)  )^ALPHA * lss^(1-ALPHA)- PHI)/vpss;
    lss - ( ((1-ALPHA)*mcss*(tilkgss^ALPHAG * tilkss^(1-ALPHAG))^ALPHA)/tilwss )^(1/ALPHA);
];

end

