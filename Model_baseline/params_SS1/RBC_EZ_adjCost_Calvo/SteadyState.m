%%% Known variables in steady state %%%
parameters
LOSS=0;

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
zass=0;
LAMBDA_X=hatzss;
RSTARss=RSS;

%%% Initial values for fsolve
x0=[1,0.33];

[x,fval,exitflag] =fsolve(@solve_SS,x0,optimset('Display','on'));

%%% Variables to find with fsolve
tilyss=x(1);
lss=x(2);


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
tiluss_2=(1-LOSS)*tilcss*(1-lss)^NU;

%%% Logs of variables to save
logtilv_tilvss_ss=0;
auxvar1ss=(logtilv_tilvss_ss+loghatzss)*(1-GAMMA);
logtiluss=log(tiluss);

logtilv_tilvss_ss_2=0;
auxvar1ss_2=(logtilv_tilvss_ss_2+loghatzss)*(1-GAMMA);
logtiluss_2=log(tiluss_2);
tilv_2_ss=tiluss_2;
logtilv_2_ss=log(tilv_2_ss);

logtilcss=log(tilcss);
loglss=log(lss);
logtillambdass=log(tillambdass);
logtilwss=log(tilwss);
logucss=log(ucss);
lognegtilulss=log(NU)+logtilcss+(NU-1)*log(1-lss);
logmnomss=log(BETA)+logtillambdass-PSI*loghatzss+(PSI-GAMMA)*logtilv_tilvss_ss+(PSI-GAMMA)*loghatzss;
logmdenss=logtillambdass+(PSI-GAMMA)/(1-GAMMA)*auxvar1ss;
logtilrss=log(tilrss);
logtilqss=log(tilqss);
logtilxss=log(tilxss);
logtilkss=log(tilkss);
logtilyss=log(tilyss);
tildivss=tilyss-tilwss*lss-tilxss;
logtilqess=(logmnomss-logmdenss+loghatzss+log(tildivss))-log(1-exp(logmnomss-logmdenss+loghatzss));
logqfss=logmnomss-logmdenss;
Sss=0;
Sprimess=0;
logtilg1ss=log(tilg1ss);
logtilkstarss=log(tilkstarss);
loggss=log(GBAR);
logtaucss=log(taucss);
logBbackss=log(Bbackss);
SCALEPARAM=1-BETA*exp(auxvar1ss)^((1-PSI)/(1-GAMMA)); %SCALEPARAM is needed to state the utility function in ratios of steady state utility
LOGTILUSS=logtiluss;
LOGTILUSS_2=logtiluss_2;

logtilxgss=log(tilxgss);
logtilkgss=log(tilkgss);
logtilkgstarss=log(tilkgstarss);
logpibackss=log(piss);
logvpbackss=log(vpss);
logpiss=log(piss);
logpistarss=log(pistarss);
nlss=log(lss)-log(1-lss);
aux2ss=logpiss+logpistarss;
logthetass=log(thetass);
logtilkgstarbackss=logtilkgstarss;
e_xgss=0;
logmcss=log(mcss);
logbyss=log(BY);
logRSTARss=log(RSTARss);
loggrantss=log(grantss);
logtilnxss=log(tilnxss);

%ratios to match
Public_Investment_GDP=100*tilxgss/tilyss;
Private_Investment_GDP=100*tilxss/tilyss;
Revenue_GDP=100*taucss*tilcss/tilyss;
% PubK_share=100*tilkgss/(tilkgss+tilkss)

