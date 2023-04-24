pibackss=PISS;
piss=PISS;

pistarss=[(1-THETA_P*(piss^CHI/piss)^(1-EPSILON))/(1-THETA_P)]^(1/(1-EPSILON));
vpss=(1-THETA_P)*pistarss^(-EPSILON)/[1-THETA_P*(piss^CHI/piss)^(-EPSILON)];
dss=MUD;
thetass=THETABAR;
loghatass=LAMBDA_A-(1-ALPHA)*dss*thetass;
% loghatmuss=LAMBDA_MU;
loghatmuss=0; % assume no investment specific shock (or trend)
loghatzss=1/(1-ALPHA)*loghatass+ALPHA/(1-ALPHA)*loghatmuss;

hatzss=exp(loghatzss);
mss=BETA*hatzss^(-PSI);
tilqss=1;

tilrss=tilqss/(mss*exp(-dss*thetass))-tilqss*(1-DELTA);

hatass=exp(loghatass);
hatmuss=exp(loghatmuss);


% solve_SS( tilwss,tilrss,vpss,dss,thetass,hatass,hatzss,hatmuss,piss,tilqss,...
%     ALPHA,DELTA,THETA_P,EPSILON,CHI,NU,PHI,GAMMA,BETA,PSI);

OPTIONS = optimoptions('fsolve','tolF',1e-10);

[tilwss,R,~,output,J]=fsolve(@(tilwss) solve_SS(tilwss,tilrss,vpss,dss,thetass,hatass,hatzss,hatmuss,piss,tilqss,...
    ALPHA,DELTA,THETA_P,EPSILON,CHI,NU,PHI,GAMMA,BETA,PSI),tilwss,OPTIONS);


[ ~,auxvar1ss,logtilv_tilvss_ss,logtiluss,logtilcss,loglss,logtillambdass,logtilwss,logucss,...
    lognegtilulss,logmnomss,logmdenss,logtilrss,logtilqss,logtilxss,logtilkss,logtilyss,logtilqess,logqfss,tildivss,...
    Sss,Sprimess,logtilg1ss,logtilkstarss, ...
    SCALEPARAM,LOGTILUSS] = solve_SS( tilwss,tilrss,vpss,dss,thetass,hatass,hatzss,hatmuss,piss,tilqss,...
    ALPHA,DELTA,THETA_P,EPSILON,CHI,NU,PHI,GAMMA,BETA,PSI);

logpibackss=log(piss);
logvpbackss=log(vpss);
logpiss=log(piss);
logpistarss=log(pistarss);
lss=exp(loglss);
nlss=log(lss)-log(1-lss);
aux2ss=logpiss+logpistarss;
