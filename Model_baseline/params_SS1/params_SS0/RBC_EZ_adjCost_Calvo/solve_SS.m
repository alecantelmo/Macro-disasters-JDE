function [ residual,auxvar1ss,logtilv_tilvss_ss,logtiluss,logtilcss,loglss,logtillambdass,logtilwss,logucss,...
    lognegtilulss,logmnomss,logmdenss,logtilrss,logtilqss,logtilxss,logtilkss,logtilyss,logtilqess,logqfss,tildivss,...
    Sss,Sprimess,logtilg1ss,logtilkstarss, ...
    SCALEPARAM,LOGTILUSS] = solve_SS( tilwss,tilrss,vpss,dss,thetass,hatass,hatzss,hatmuss,piss,tilqss,...
    ALPHA,DELTA,THETA_P,EPSILON,CHI,NU,PHI,GAMMA,BETA,PSI)
%Auxiliary function to use with fsolve to find tilwss for the Calvo pricing
%model.

loghatzss=log(hatzss);
k_to_l=ALPHA/(1-ALPHA)*tilwss/tilrss;

kstar_to_l=k_to_l*exp(dss*thetass)*hatzss*hatmuss;

lss=(PHI+tilwss/NU*vpss)/(hatass/hatzss*kstar_to_l^ALPHA*exp(-dss*thetass*ALPHA)...
    +tilwss/NU*vpss-(kstar_to_l-(1-DELTA)*k_to_l)*vpss);

tilkstarss=kstar_to_l*lss;

tilkss=tilkstarss/hatzss/hatmuss*exp(-dss*thetass);
tilxss=tilkstarss-(1-DELTA)*tilkss;

tilcss=(1-lss)*tilwss/NU;
tilyss=tilcss+tilxss;

% logs
logtilcss=log(tilcss);

% continue

logtiluss=logtilcss+NU*log(1-lss);
ucss=(1-lss)^NU;
tilulss=-NU*tilcss*(1-lss)^(NU-1);

logthetass=log(thetass);
zass=0;

logtilv_tilvss_ss=0;
auxvar1ss=(logtilv_tilvss_ss+loghatzss)*(1-GAMMA);
SCALEPARAM=1-BETA*exp(auxvar1ss)^((1-PSI)/(1-GAMMA)); %SCALEPARAM is needed to state the utility function in ratios of steady state utility

loglss=log(lss);

logucss=log(ucss);
logtillambdass=log(1-PSI)-PSI*logtiluss+logucss;
lognegtilulss=log(NU)+logtilcss+(NU-1)*log(1-lss);


logmnomss=log(BETA)+logtillambdass-PSI*loghatzss+(PSI-GAMMA)*logtilv_tilvss_ss+(PSI-GAMMA)*loghatzss;
logmdenss=logtillambdass+(PSI-GAMMA)/(1-GAMMA)*auxvar1ss;

logtilwss=log(tilwss);
logtilrss=log(tilrss);
logtilqss=log(tilqss);
logtilxss=log(tilxss);
logtilkss=log(tilkss);
logtilkstarss=log(tilkstarss);
logtilyss=log(tilyss);

LOGTILUSS=logtiluss;

% logqeyss=(logmnomss-logmdenss+loghatzss)-log(1-exp(logmnomss-logmdenss+loghatzss));
tildivss=tilyss-tilwss*lss-tilxss;

logtilqess=(logmnomss-logmdenss+loghatzss+log(tildivss))-log(1-exp(logmnomss-logmdenss+loghatzss));

logqfss=logmnomss-logmdenss;

%%%%%% continue
mss=exp(logmnomss-logmdenss);

mcss=(1/(1-ALPHA))^(1-ALPHA)*(1/ALPHA)^ALPHA*tilwss^(1-ALPHA)*tilrss^ALPHA;

tilg1ss=mcss*tilyss/[1-THETA_P*mss*(piss^CHI/piss)^(-EPSILON)*hatzss];

tilg2ss=EPSILON*tilg1ss/(EPSILON-1);

pistarss=[tilg2ss-THETA_P*mss*(piss^CHI/piss)^(1-EPSILON)*tilg2ss*hatzss]/tilyss;

residual=-1+THETA_P*(piss^CHI/piss)^(1-EPSILON)+(1-THETA_P)*pistarss^(1-EPSILON);

Sss=0;
Sprimess=0;
logtilg1ss=log(tilg1ss);

end

