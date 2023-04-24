function [ residual,auxvar1ss,logtilv_tilvss_ss,logtiluss,logtilcss,loglss,logtillambdass,logtilwss,logucss,...
    lognegtilulss,logmnomss,logmdenss,logtilrss,logtilqss,logtilxss,logtilkss,logtilyss,logtilqess,logqfss,tildivss,...
    Sss,Sprimess,logtilkstarss,logtaucss,logBbackss,...
    logthetass,zass,nlss,loggss, tilwss,...
    SCALEPARAM,LOGTILUSS] = solve_SS(lss,tilrss,dss,thetass,hatass,hatzss,hatmuss,tilqss,...
    ALPHA,DELTA,NU,PHI,GAMMA,BETA,PSI,GBAR,GY,logtauc_ss,logBbackss_ss,BY,RHOTAUB,RHOTAU)


loghatzss=log(hatzss);

taucss=0.20;%GY-(1-RSS/(hatzss*piss))*BY;
logtaucss=log(taucss);
logtauc_ss=logtaucss;

kstar_to_l=(ALPHA*hatass*hatmuss*exp(-dss*thetass*(ALPHA-1))/tilrss)^(1/(1-ALPHA));
tilwss=(1-ALPHA)*hatass/hatzss*(kstar_to_l*exp(-dss*thetass))^ALPHA;

tilkstarss=kstar_to_l*lss;

tilkss=tilkstarss/hatzss/hatmuss*exp(-dss*thetass);
tilxss=tilkstarss-(1-DELTA)*tilkss;

tilcss=(1-lss)*tilwss/NU;
% tilyss=tilcss+tilxss;
tilyss=(tilcss+tilxss)/(1-GY);
GBAR=GY*tilyss;
loggss=log(GBAR);

Bbackss=BY*tilyss;
logBbackss=log(Bbackss);
logBbackss_ss=logBbackss;
%test


% tilrss-ALPHA*hatass*hatmuss*(tilkstarss*exp(-dss*thetass))^(ALPHA-1)*lss^(1-ALPHA)
% tilwss-(1-ALPHA)*hatass/hatzss*(tilkstarss*exp(-dss*thetass))^(ALPHA)*lss^(-ALPHA)

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
logtillambdass=log(1-PSI)-PSI*logtiluss+logucss-log(1+exp(logtaucss));
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

tildivss=tilyss-tilwss*lss-tilxss;

logtilqess=(logmnomss-logmdenss+loghatzss+log(tildivss))-log(1-exp(logmnomss-logmdenss+loghatzss));

logqfss=logmnomss-logmdenss;

nlss=log(lss)-log(1-lss);

%Define auxiliary variables for lss
ghelp=1/(1-GY);
whelp=tilwss/(NU*(1+taucss));  %Only this helping variable has to change to determine lss when the consumption tax is introduced 
dhelp=exp(-dss*thetass);
deltahelp=(1-DELTA)/hatzss/hatmuss*dhelp;


residual=lss-((PHI+ghelp*whelp)/(hatass/hatzss*(kstar_to_l*dhelp)^ALPHA-ghelp*kstar_to_l*(1-deltahelp)+ghelp*whelp));


% residual=-lss+(PHI+tilwss/NU+GBAR)/(hatass/hatzss*kstar_to_l^ALPHA*exp(-dss*thetass*ALPHA)...
%     +tilwss/NU-kstar_to_l*(1-(1-DELTA)/hatzss/hatmuss*exp(-dss*thetass)));

Sss=0;
Sprimess=0;


end

