syms d dp logtheta logthetap za zap loghata loghatap auxvar1 auxvar1p real
syms tilvp_tilvss tilv_tilvss hatzp hatz real
syms logtilv_tilvss logtilvp_tilvss real
syms tilu tilup tilc tilcp l lp real
syms logtilu logtilup logtilc logtilcp real
syms loguc logucp lognegtilul lognegtilulp real
syms logtillambda logtillambdap logtilw logtilwp real
syms logmden logmdenp logmnom logmnomp real
syms logtilq logtilqp logtilr logtilrp real
syms logtilx logtilxp logtilxback logtilxbackp real
syms logtilkstarback logtilkstarbackp real
syms logtilk logtilkp real
syms logl loglp real
syms logtily logtilyp real
syms logqey logqeyp real
syms logMp u17 real

syms MUD RHOTHETA THETABAR SDV_THETA SDV_ZA ALPHA LAMBDA_A GAMMA NU PSI BETA DELTA PHI SCALEPARAM LOGTILUSS KAPPA LAMBDA_X  real
syms logg loggp GBAR RHOG GY real 
syms logtauc logtaucp logBback logBbackp real
syms logtauc_ss logBbackss_ss BY RHOTAUB LUMPSUM real
syms logtaucback RHOTAU real

%add variables for public capital
syms e_xg e_xgp tilxg logtilxg tilxgp logtilxgp logtilkgstarback logtilkgstarbackp real
syms logtilkg logtilkgp logmcp logtilxgback logtilxgbackp real
  
%add parameters for public capital
syms SDV_XG DELTAG logtilxg_ss RHOXG logtilkg_ss ALPHAG XGY logtily_ss logtilkgstarback_ss RHOKG real

syms RHOTAUY RHOGB RHOGY real
syms logby logbyp real
syms RHOZA real
syms logRSTAR logRSTARp  ETA real
syms logtilnx logtilnxp tilnx tilnxp NXY real
syms loggrant loggrantp RHOGRANT RHOGRANT_ETA GRANTY loggrant_ss real% tilgrant tilgrantp real

% Exogenous state variables

x2p=[dp;
    logthetap;
    zap;

    e_xgp]; %Alessandro


Phi_fun=[MUD;
    (1-RHOTHETA)*log(THETABAR)+RHOTHETA*logtheta;
    RHOZA*za;
    0];

eta_mat=[0,0,0,0;
         0,0,0,0;
         0,0,0,0;
         0,0,0,0;
         0,0,0,0;
         0,0,0,0;
         0,0,0,0;
         0,0,0,0;
         0,0,0,0;
         0,0,0,0;
         1,0,0,0;
         0,SDV_THETA,0,0;
         0,0,SDV_ZA,0;
         0,0,0,SDV_XG];


     
% other variables that depend only on the exogenous state vars

syms loghata loghatap loghatz loghatzp real
loghata_=LAMBDA_A+za-(1-ALPHA)*d*exp(logtheta);
loghatap_=LAMBDA_A+zap-(1-ALPHA)*dp*exp(logthetap);

 % assume no investment specific shock
loghatmu=0; loghatmup=0;
%%%

loghatz_=1/(1-ALPHA)*loghata+ALPHA/(1-ALPHA)*loghatmu;
loghatzp_=1/(1-ALPHA)*loghatap+ALPHA/(1-ALPHA)*loghatmup;

% define auxiliary variable to bind l within [0,1]
syms nl nlp real

logl=nl-log(1+exp(nl));
loglp=nlp-log(1+exp(nlp));

syms  LOSS HATZ_ss auxvar1_2 auxvar1p_2 logtilv_tilvss_2 logtilvp_tilvss_2 logtilu_2 logtilup_2 LOGTILUSS_2 logtilv_2 logtilvp_2 tilv_2_ss real
syms u01p_2 u02_2 real
u01p_2_=logtilvp_tilvss_2+loghatzp;
u02_2_=exp(-auxvar1_2);
f37=-1+exp(u01p_2)^(1-GAMMA)*u02_2; % auxiliary variable


syms u1_2 u2_2 u1p_2 u2p_2 real

u1_2_=exp(u2_2)*SCALEPARAM+BETA*exp(auxvar1_2*((1-PSI)/(1-GAMMA)));
u1p_2_=exp(u2p_2)*SCALEPARAM+BETA*exp(auxvar1p_2*((1-PSI)/(1-GAMMA)));

u2_2_=(logtilu_2-LOGTILUSS_2)*(1-PSI);
u2p_2_=(logtilup_2-LOGTILUSS_2)*(1-PSI);

logtilv_tilvss_2_=1/(1-PSI)*log(u1_2);
logtilvp_tilvss_2_=1/(1-PSI)*log(u1p_2);

logtilu_2_=log(1-LOSS)+logtilc+NU*log(1-l);
logtilup_2_=log(1-LOSS)+logtilcp+NU*log(1-lp);

logtilv_2_=logtilv_tilvss_2+log(tilv_2_ss);
logtilvp_2_=logtilvp_tilvss_2+log(tilv_2_ss);

syms u01p u02 real
u01p_=logtilvp_tilvss+loghatzp;
u02_=exp(-auxvar1);
f0=-1+exp(u01p)^(1-GAMMA)*u02; % auxiliary variable

syms u1 u2 u1p u2p real

u1_=exp(u2)*SCALEPARAM+BETA*exp(auxvar1*((1-PSI)/(1-GAMMA)));
u1p_=exp(u2p)*SCALEPARAM+BETA*exp(auxvar1p*((1-PSI)/(1-GAMMA)));

u2_=(logtilu-LOGTILUSS)*(1-PSI);
u2p_=(logtilup-LOGTILUSS)*(1-PSI);

logtilv_tilvss_=1/(1-PSI)*log(u1);
logtilvp_tilvss_=1/(1-PSI)*log(u1p);

l_=exp(logl);
lp_=exp(loglp);

logtilu_=logtilc+NU*log(1-l);
logtilup_=logtilcp+NU*log(1-lp);

loguc_=NU*log(1-l);
logucp_=NU*log(1-lp);

logtillambda_=log(1-PSI)-PSI*logtilu+loguc-log(1+exp(logtaucback));
logtillambdap_=log(1-PSI)-PSI*logtilup+logucp-log(1+exp(logtauc));

% use f3, f4, f5, f6 to solve for logtilc

logtilc_=logtilw+log(1-l)-log(NU);
logtilcp_=logtilwp+log(1-lp)-log(NU);

logmnom_=log(BETA)+logtillambda-PSI*loghatz+(PSI-GAMMA)*logtilv_tilvss+(PSI-GAMMA)*loghatz;
logmnomp_=log(BETA)+logtillambdap-PSI*loghatzp+(PSI-GAMMA)*logtilvp_tilvss+(PSI-GAMMA)*loghatzp;

logmden_=logtillambda+(PSI-GAMMA)/(1-GAMMA)*auxvar1; %denominator of logmp

syms tilr tilrp tilq tilqp theta thetap real

tilrp_=exp(logtilrp);

theta_=exp(logtheta);
thetap_=exp(logthetap);


syms u03 u04 real
logMp_=logmnomp-logmden;
u03_=logMp-dp*thetap-loghatmup;
u04_=1/tilq;
f9=exp(u03)*(tilrp+tilqp*(1-DELTA))*u04-1;

syms S Sprime logtilxback Sp Sprimep logtilxbackp real
S_=KAPPA/2*(exp(logtilx-logtilxback+loghatz)-LAMBDA_X)^2;
Sprime_=KAPPA*(exp(logtilx-logtilxback+loghatz)-LAMBDA_X);
Sprimep_=KAPPA*(exp(logtilxp-logtilxbackp+loghatzp)-LAMBDA_X);

f9c=-logtilx+logtilxbackp;

f10=-1+exp(logtilq)*(1-S-Sprime*exp(logtilx-logtilxback+loghatz))+...
    exp(logmnomp-logmden+logtilqp)*Sprimep*exp(2*(logtilxp-logtilx+loghatzp));

tilq_=exp(logtilq);
tilqp_=exp(logtilqp);


% % use f13 and f14 to solve for tilx
syms u4 tilx logtilx tily logtily real
tilc_=exp(logtilc);
tilx_=exp(logtilx);
tilxg_=exp(logtilxg);
tily_=tilc+tilx+exp(logg)+tilxg-tilnx;%+tilgrant;
logtily_=log(tily);
tilnx_=exp(logtilnx);

% 
syms u4p tilyp logtilyp tilxp tilcp real
tilcp_=exp(logtilcp);
tilxp_=exp(logtilxp);
tilxgp_=exp(logtilxgp);
tilyp_=tilcp+tilxp+exp(loggp)+tilxgp-tilnxp;
logtilyp_=log(tilyp);
tilnxp_=exp(logtilnxp);



syms u05 real
u05_=exp(-logtilx);
f15=exp(logtilkstarbackp)*u05-(1-DELTA)*exp(logtilk)*u05-(1-S);

logtilk_=logtilkstarback-d*theta-loghatz-loghatmu;
logtilkp_=logtilkstarbackp-dp*thetap-loghatzp-loghatmup;

syms w wp real
syms tildiv tildivp logtilqe logtilqep real

w_=exp(logtilw);
wp_=exp(logtilwp);

u17_=logMp+loghatzp;

syms logqf logqfp real

% Calvo Pricing %%%%%%%%%%%%%
syms logmc logtilg1 logtilg1p logtilg2 logtilg2p logpi logpip logpistar logpistarp logpiback logpibackp real
syms logvpback logvpbackp real
syms THETA_P CHI EPSILON real

% FOC

f19=-1+[exp(logmc+logtily)+THETA_P*exp(logmnomp-logmden-EPSILON*(CHI*logpi-logpip)+logtilg1p+loghatzp)]*exp(-logtilg1);

logtilg2_=log(EPSILON)+logtilg1-log(EPSILON-1);
logtilg2p_=log(EPSILON)+logtilg1p-log(EPSILON-1);

syms aux2 aux2p pi_power pi_powerp real

pi_power_=THETA_P*exp(logpiback*(CHI*(1-EPSILON)))+(1-THETA_P)*exp(aux2*(1-EPSILON));
pi_powerp_=THETA_P*exp(logpibackp*(CHI*(1-EPSILON)))+(1-THETA_P)*exp(aux2p*(1-EPSILON));

logpi_=log(pi_power)/(1-EPSILON);
logpip_=log(pi_powerp)/(1-EPSILON);

logpistar_=aux2-logpi;
logpistarp_=aux2p-logpip;

f20=-1+[exp(logpistar+logtily)+THETA_P*exp(logmnomp-logmden+(1-EPSILON)*(CHI*logpi-logpip)...
    +logpistar-logpistarp+logtilg2p+loghatzp)]*exp(-logtilg2);

f22=-logpibackp+logpi;

% Price dispersion
f23=-1+[THETA_P*exp(-EPSILON*(CHI*logpiback-logpi)+logvpback)+(1-THETA_P)*exp(-EPSILON*logpistar)]*exp(-logvpbackp);

% Optimal factor composition
logtilr_=log(ALPHA)+log(1-ALPHAG)+logmc+ALPHA*(ALPHAG*logtilkg+(1-ALPHAG)*logtilk)-logtilk+(1-ALPHA)*logl;
logtilrp_=log(ALPHA)+log(1-ALPHAG)+logmcp+ALPHA*(ALPHAG*logtilkgp+(1-ALPHAG)*logtilkp)-logtilkp+(1-ALPHA)*loglp;

logmc_=logtilw-(log(1-ALPHA)+ALPHA*(ALPHAG*logtilkg+(1-ALPHAG)*logtilk)-ALPHA*logl);
logmcp_=logtilwp-(log(1-ALPHA)+ALPHA*(ALPHAG*logtilkgp+(1-ALPHAG)*logtilkp)-ALPHA*loglp);

% Aggregation 
syms u13a real
u13a_=exp(logtily+logvpbackp)+PHI;
f13=-log(u13a)+loghata-loghatz+(ALPHA)*(ALPHAG*(logtilkgstarback-d*exp(logtheta))+ (1-ALPHAG)*(logtilkstarback-d*exp(logtheta)) )+(1-ALPHA)*logl;

% End of Calvo Pricing %%%%%%%%%%%%%%%

% Simple Taylor rule
syms logR logRp RSS PISS GAMMA_PI real
logR_=log(RSS)+GAMMA_PI*(logpi-log(PISS));
f24=-1+exp(logmnomp-logmden+logR-logpip);


% Government budget
f25=-exp(logBbackp)+exp(logRSTAR)*exp(logBback)/(exp(loghatz)*exp(logpi))+exp(logg)+exp(logtilxg)-exp(logtaucback)*exp(logtilc)-LUMPSUM-exp(loggrantp);
f26=-(logtauc-logtauc_ss)+RHOTAU*(logtaucback-logtauc_ss)+RHOTAUB*(logBback-logBbackss_ss)+RHOTAUY*(logtily-logtily_ss);
f32=-(loggp-log(GBAR))+RHOG*(logg-log(GBAR))-RHOGB*(logBback-logBbackss_ss)-RHOGY*(logtily-logtily_ss);
logby_=logBback-log(4)-logtily;
logbyp_=logBbackp-log(4)-logtilyp;

logRSTAR_=log(RSS)+ETA*( exp(logBback) / exp(logBbackss_ss) -1);
logRSTARp_=log(RSS)+ETA*( exp(logBbackp) / exp(logBbackss_ss)-1);

f35=-(loggrantp-loggrant_ss)+RHOGRANT*(loggrant-loggrant_ss-loghatz+log(LAMBDA_X))+(1-RHOGRANT)*RHOGRANT_ETA*(   d*exp(logtheta) - MUD*THETABAR );

f36=exp(logBbackp)-exp(logtilnx)-exp(logRSTAR)*exp(logBback)/(exp(loghatz)*exp(logpi))+REMIT+exp(loggrantp);

% Add public capital block
f28=exp(logtilkgstarbackp)-(1-DELTAG)*exp(logtilkg)-exp(logtilxg);%-exp(loggrantp);
logtilkg_=logtilkgstarback-d*theta-loghatz-loghatmu;
logtilkgp_=logtilkgstarbackp-dp*thetap-loghatzp-loghatmup;
f30=-(logtilxgbackp-logtilxg_ss)+RHOXG*(logtilxgback-logtilxg_ss)+RHOKG*(d*exp(logtheta)-MUD*THETABAR)+e_xg;%might work for large rhokg
f31=-logtilxg+logtilxgbackp;


f=[f0;f9;f13;f15;f9c;f10;f19;f20;f22;f23;f24;f25;f26;f28;f30;f31;f32;f35;f36;f37];%
f=simplify(f);

x=[logtilkstarback;
    logtilxback;
    logpiback;
    logvpback;
    logBback;
    logtaucback;
    logtilkgstarback;
    logtilxgback;
    logg;
    loggrant;
    d;
    logtheta;
    za;
    e_xg];

xp=[logtilkstarbackp;
    logtilxbackp;
    logpibackp;
    logvpbackp;
    logBbackp;
    logtauc;
    logtilkgstarbackp;
    logtilxgbackp;
    loggp;
    loggrantp;
    dp;
    logthetap;
    zap;
    e_xgp];


y=[auxvar1,nl,logtilw,logtilx,logtilq,logtilg1,aux2,logtilxg,logtilnx,auxvar1_2];
yp=[auxvar1p,nlp,logtilwp,logtilxp,logtilqp,logtilg1p,aux2p,logtilxgp,logtilnxp,auxvar1p_2];

subsvars=[loghata;loghatap;loghatz;loghatzp;u1;u1p;u2;u2p;logtilv_tilvss;logtilvp_tilvss;l;lp;logtilu;logtilup;loguc;logucp;...
    logtillambda;logtillambdap;logtilc;logtilcp;logmnom;logmnomp;logmden;...
    theta;thetap;tilrp;tilq;tilqp;logtilr;logtilrp;...
    tilc;tilcp;tilx;tilxp;logtilk;logtilkp;u01p;u02;u03;u04;u05;logMp;u17;tily;logtily;tilyp;logtilyp;...
    w;wp;S;Sprime;Sprimep;logtilg2;logtilg2p;logpistar;logpistarp;...
    u13a;logR;pi_power;pi_powerp;logpi;logpip;tilxg;tilxgp;logtilkg;logtilkgp;logmc;logmcp;logby;logbyp;logRSTAR;logRSTARp;tilnx;tilnxp;...
    u1_2;u1p_2;u2_2;u2p_2;logtilv_tilvss_2;logtilvp_tilvss_2;logtilu_2;logtilup_2;u01p_2;u02_2;logtilv_2;logtilvp_2];


subsfuns=sym(zeros(size(subsvars)));
for i=1:length(subsfuns)
    subsfuns(i)=eval([char(subsvars(i)) '_']);
end

symparams=[MUD RHOTHETA THETABAR SDV_THETA SDV_ZA ALPHA LAMBDA_A GAMMA ...
    NU PSI BETA DELTA PHI SCALEPARAM LOGTILUSS KAPPA LAMBDA_X...
    THETA_P CHI EPSILON RSS PISS GAMMA_PI RHOG GY GBAR logtauc_ss logBbackss_ss BY RHOTAUB LUMPSUM RHOTAU...
    SDV_XG DELTAG logtilxg_ss RHOXG logtilkg_ss ALPHAG XGY logtily_ss logtilkgstarback_ss RHOKG RHOTAUY RHOGB RHOGY RHOZA  ETA REMIT NXY RHOGRANT RHOGRANT_ETA GRANTY loggrant_ss LOSS HATZ_ss LOGTILUSS_2 tilv_2_ss];%LOGTILWELFSS

