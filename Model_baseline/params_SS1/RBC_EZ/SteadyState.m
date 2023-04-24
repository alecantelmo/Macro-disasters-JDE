
dss=MUD;
thetass=THETABAR;
loghatass=LAMBDA_A-(1-ALPHA)*dss*thetass;

loghatmuss=0; % assume no investment specific shock (or trend)
loghatzss=1/(1-ALPHA)*loghatass+ALPHA/(1-ALPHA)*loghatmuss;

hatzss=exp(loghatzss);

mss=BETA*hatzss^(-PSI);
tilqss=1;

tilrss=tilqss/(mss*exp(-dss*thetass))-tilqss*(1-DELTA);

hatass=exp(loghatass);
hatmuss=exp(loghatmuss);

lss=0.30; %initial value for fsolve

% now use fsolve (Alessandro)
OPTIONS = optimoptions('fsolve','tolF',1e-10);
[lss,GN,~,output,J]=fsolve(@(lss) solve_SS(lss,tilrss,dss,thetass,hatass,hatzss,hatmuss,tilqss,...
    ALPHA,DELTA,NU,PHI,GAMMA,BETA,PSI,GBAR,GY,logtauc_ss,logBbackss_ss,BY,RHOTAUB,RHOTAU),lss,OPTIONS);

[ ~,auxvar1ss,logtilv_tilvss_ss,logtiluss,logtilcss,loglss,logtillambdass,logtilwss,logucss,...
    lognegtilulss,logmnomss,logmdenss,logtilrss,logtilqss,logtilxss,logtilkss,logtilyss,logtilqess,logqfss,tildivss,...
    Sss,Sprimess,logtilkstarss,logtaucss,logBbackss,...
    logthetass,zass,nlss,loggss,tilwss, ...
    SCALEPARAM,LOGTILUSS] = solve_SS(lss,tilrss,dss,thetass,hatass,hatzss,hatmuss,tilqss,...
    ALPHA,DELTA,NU,PHI,GAMMA,BETA,PSI,GBAR,GY,logtauc_ss,logBbackss_ss,BY,RHOTAUB,RHOTAU);

%now comment all the remaining equations and copy them in solve_SS

% kstar_to_l=(ALPHA*hatass*hatmuss*exp(-dss*thetass*(ALPHA-1))/tilrss)^(1/(1-ALPHA));
% tilwss=(1-ALPHA)*hatass/hatzss*(kstar_to_l*exp(-dss*thetass))^ALPHA;

% lss=(PHI+tilwss/NU+GBAR)/(hatass/hatzss*kstar_to_l^ALPHA*exp(-dss*thetass*ALPHA)...
%     +tilwss/NU-kstar_to_l*(1-(1-DELTA)/hatzss/hatmuss*exp(-dss*thetass)));
% 
% tilkstarss=kstar_to_l*lss;
% 
% tilkss=tilkstarss/hatzss/hatmuss*exp(-dss*thetass);
% tilxss=tilkstarss-(1-DELTA)*tilkss;
% 
% tilcss=(1-lss)*tilwss/NU;
% % tilyss=tilcss+tilxss;
% tilyss=(tilcss+tilxss)/(1-GY);
% GBAR=GY*tilyss;
% loggss=log(GBAR);
% %test
% 
% % tilrss-ALPHA*hatass*hatmuss*(tilkstarss*exp(-dss*thetass))^(ALPHA-1)*lss^(1-ALPHA)
% % tilwss-(1-ALPHA)*hatass/hatzss*(tilkstarss*exp(-dss*thetass))^(ALPHA)*lss^(-ALPHA)
% 
% % logs
% logtilcss=log(tilcss);
% 
% % continue
% 
% logtiluss=logtilcss+NU*log(1-lss);
% ucss=(1-lss)^NU;
% tilulss=-NU*tilcss*(1-lss)^(NU-1);
% 
% logthetass=log(thetass);
% zass=0;
% 
% logtilv_tilvss_ss=0;
% auxvar1ss=(logtilv_tilvss_ss+loghatzss)*(1-GAMMA);
% SCALEPARAM=1-BETA*exp(auxvar1ss)^((1-PSI)/(1-GAMMA)); %SCALEPARAM is needed to state the utility function in ratios of steady state utility
% 
% loglss=log(lss);
% 
% logucss=log(ucss);
% logtillambdass=log(1-PSI)-PSI*logtiluss+logucss;
% lognegtilulss=log(NU)+logtilcss+(NU-1)*log(1-lss);
% 
% 
% logmnomss=log(BETA)+logtillambdass-PSI*loghatzss+(PSI-GAMMA)*logtilv_tilvss_ss+(PSI-GAMMA)*loghatzss;
% logmdenss=logtillambdass+(PSI-GAMMA)/(1-GAMMA)*auxvar1ss;
% 
% logtilwss=log(tilwss);
% logtilrss=log(tilrss);
% logtilqss=log(tilqss);
% logtilxss=log(tilxss);
% logtilkss=log(tilkss);
% logtilkstarss=log(tilkstarss);
% logtilyss=log(tilyss);
% 
% LOGTILUSS=logtiluss;
% 
% tildivss=tilyss-tilwss*lss-tilxss;
% 
% logtilqess=(logmnomss-logmdenss+loghatzss+log(tildivss))-log(1-exp(logmnomss-logmdenss+loghatzss));
% 
% logqfss=logmnomss-logmdenss;
% 
% nlss=log(lss)-log(1-lss);