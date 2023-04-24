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

syms MUD RHOTHETA THETABAR SDV_THETA SDV_ZA ALPHA LAMBDA_A GAMMA NU PSI BETA DELTA PHI SCALEPARAM LOGTILUSS KAPPA LAMBDA_X real
syms logg loggp GBAR RHOG GY real 

syms logtauc logtaucp logBback logBbackp real
syms logtauc_ss logBbackss_ss BY RHOTAUB LUMPSUM real
syms logtaucback RHOTAU real

%add variables for public capital
syms e_xg e_xgp logtilxg logtilxgp logtilkgstarbackp logtilkg logtilkgp logtilkgstarback logtilxgback logtilxgbackp real

%add parameters for public capital
syms SDV_XG DELTAG logtilxg_ss RHOXG logtilkg_ss ALPHAG XGY logtily_ss logtilkgstarback_ss RHOKG real
syms RHOTAUY RHOGB RHOGY real
syms logby logbyp real
syms RHOZA real

syms logRSTAR logRSTARp  ETA real
syms logtilnx logtilnxp REMIT NXY real
syms loggrant loggrantp RHOGRANT RHOGRANT_ETA GRANTY loggrant_ss  real

% Exogenous state variables

x2p=[dp;
    logthetap;
    zap;
    e_xgp]; 


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
loghata=LAMBDA_A+za-(1-ALPHA)*d*exp(logtheta);
loghatap=LAMBDA_A+zap-(1-ALPHA)*dp*exp(logthetap);
loghatmu=0; loghatmup=0; % assume no investment specific shock
loghatz=1/(1-ALPHA)*loghata+ALPHA/(1-ALPHA)*loghatmu;
loghatzp=1/(1-ALPHA)*loghatap+ALPHA/(1-ALPHA)*loghatmup;

% define auxiliary variable nl to bind l within [0,1]
syms nl nlp real
l=exp(nl)/(1+exp(nl));
lp=exp(nlp)/(1+exp(nlp));
logl=nl-log(1+exp(nl));
loglp=nlp-log(1+exp(nlp));

syms  LOSS HATZ_ss auxvar1_2 auxvar1p_2  logtilv_tilvss_2 logtilvp_tilvss_2 logtilu_2 logtilup_2 LOGTILUSS_2 logtilv_2 logtilvp_2 tilv_2_ss real
f37=-1+exp(-auxvar1_2)*exp(logtilvp_tilvss_2+loghatzp)^(1-GAMMA); % auxiliary variable - VERY IMPORTANT to do it in log
f38=[-exp(logtilv_tilvss_2)^(1-PSI)+exp(logtilu_2-LOGTILUSS_2)^(1-PSI)*SCALEPARAM+BETA*exp(auxvar1_2*(1-PSI)/(1-GAMMA))]*exp(-logtilv_tilvss_2)^(1-PSI);
f39=-logtilu_2+log(1-LOSS)+logtilc+NU*log(1-l);
f40=-logtilv_2+logtilv_tilvss_2+log(tilv_2_ss);

f0=-1+exp(-auxvar1)*exp(logtilvp_tilvss+loghatzp)^(1-GAMMA); % auxiliary variable - VERY IMPORTANT to do it in log
f1=[-exp(logtilv_tilvss)^(1-PSI)+exp(logtilu-LOGTILUSS)^(1-PSI)*SCALEPARAM+BETA*exp(auxvar1*(1-PSI)/(1-GAMMA))]*exp(-logtilv_tilvss)^(1-PSI);
f2=-logtilu+logtilc+NU*log(1-l);

f3=-loguc+NU*log(1-l);
f4=-lognegtilul+log(NU)+logtilc+(NU-1)*log(1-l); %lognegtilul is the log of negative tilul, because tilul itself is negative
f5=log(1-PSI)-PSI*logtilu+lognegtilul-logtillambda-logtilw; % note the use again of lognegtilul (and not logtilul)
f6=log(1-PSI)-PSI*logtilu+loguc-log(1+exp(logtaucback))-logtillambda;

f7=-logmnom+log(BETA)+logtillambda-PSI*loghatz+(PSI-GAMMA)*logtilv_tilvss+(PSI-GAMMA)*loghatz; %nominator of logm


f8=-logmden+logtillambda+(PSI-GAMMA)/(1-GAMMA)*auxvar1; %denominator of logmp

f9=exp(logmnomp-logmden-dp*exp(logthetap)-loghatmup)*(exp(logtilrp)+exp(logtilqp)*(1-DELTA))*exp(-logtilq)-1;

syms S Sprime logtilxback Sp Sprimep logtilxbackp real
f9a=-S+KAPPA/2*(exp(logtilx-logtilxback+loghatz)-LAMBDA_X)^2;
f9b=-Sprime+KAPPA*(exp(logtilx-logtilxback+loghatz)-LAMBDA_X);
f9c=-logtilx+logtilxbackp;


f10=-1+exp(logtilq)*(1-S-Sprime*exp(logtilx-logtilxback+loghatz))+...
    exp(logmnomp-logmden+logtilqp)*Sprimep*exp(2*(logtilxp-logtilx+loghatzp));

f14=-exp(logtily)+exp(logtilc)+exp(logtilx)+exp(logg)+exp(logtilxg)-exp(logtilnx);

f15=exp(logtilkstarbackp)*exp(-logtilx)-(1-DELTA)*exp(logtilk)*exp(-logtilx)-(1-S);
f16=-logtilk+logtilkstarback-d*exp(logtheta)-loghatz-loghatmu;

f28=exp(logtilkgstarbackp)-(1-DELTAG)*exp(logtilkg)-exp(logtilxg);
f29=-logtilkg+logtilkgstarback-d*exp(logtheta)-loghatz-loghatmu;

f30=-(logtilxgbackp-logtilxg_ss)+RHOXG*(logtilxgback-logtilxg_ss)+RHOKG*(d*exp(logtheta)-MUD*THETABAR)+e_xg;%might work for large rhokg


f31=-logtilxg+logtilxgbackp;


% Calvo Pricing %%%%%%%%%%%%%
syms logmc logmcp logtilg1 logtilg1p logtilg2 logtilg2p logpi logpip logpistar logpistarp logpiback logpibackp real
syms logvpback logvpbackp real
syms THETA_P CHI EPSILON real

% FOC

f11=-logtilr+log(ALPHA)+log(1-ALPHAG)+logmc+ALPHA*(ALPHAG*logtilkg+(1-ALPHAG)*logtilk)-logtilk+(1-ALPHA)*logl;
f12=-logtilw+log(1-ALPHA)+logmc+ALPHA*(ALPHAG*logtilkg+(1-ALPHAG)*logtilk)-ALPHA*logl;

f19=[-exp(logtilg1)+exp(logmc+logtily)+THETA_P*exp(logmnomp-logmden-EPSILON*(CHI*logpi-logpip)+logtilg1p+loghatzp)]*exp(-logtilg1);

logtilg2=log(EPSILON)+logtilg1-log(EPSILON-1);
logtilg2p=log(EPSILON)+logtilg1p-log(EPSILON-1);

f20=[-exp(logtilg2)+exp(logpistar+logtily)+THETA_P*exp(logmnomp-logmden+(1-EPSILON)*(CHI*logpi-logpip)...
    +logpistar-logpistarp+logtilg2p+loghatzp)]*exp(-logtilg2);

f21=-1+THETA_P*exp((1-EPSILON)*(CHI*logpiback-logpi))+(1-THETA_P)*exp((1-EPSILON)*logpistar);
syms aux2 aux2p real
f21b=aux2-logpi-logpistar;

f22=-logpibackp+logpi;

% Price dispersion
f23=[-exp(logvpbackp)+THETA_P*exp(-EPSILON*(CHI*logpiback-logpi)+logvpback)+(1-THETA_P)*exp(-EPSILON*logpistar)]*exp(-logvpbackp);


% PRODUCTION FUNCTION 
f13=-log(exp(logtily+logvpbackp)+PHI)+loghata-loghatz+(ALPHA)*(ALPHAG*(logtilkgstarback-d*exp(logtheta))+ (1-ALPHAG)*(logtilkstarback-d*exp(logtheta)) )+(1-ALPHA)*logl;

% End of Calvo Pricing %%%%%%%%%%%%%%%

% Simple Taylor rule
syms logR logRp RSS PISS GAMMA_PI real
logR=log(RSS)+GAMMA_PI*(logpi-log(PISS));
f24=-1+exp(logmnomp-logmden+logR-logpip);

% Government budget
f25=-exp(logBbackp)+exp(logRSTAR)*exp(logBback)/(exp(loghatz)*exp(logpi))+exp(logg)+exp(logtilxg)-exp(logtaucback)*exp(logtilc)-LUMPSUM-exp(loggrantp);
f26=-(logtauc-logtauc_ss)+RHOTAU*(logtaucback-logtauc_ss)+RHOTAUB*(logBback-logBbackss_ss)+RHOTAUY*(logtily-logtily_ss);
f32=-(loggp-log(GBAR))+RHOG*(logg-log(GBAR))-RHOGB*(logBback-logBbackss_ss)-RHOGY*(logtily-logtily_ss);
f33=-logby+(logBback)-log(4)-logtily;


f34=-logRSTAR+log(RSS)+ETA* (exp(logBback)/ exp(logBbackss_ss)-1);

f35=-(loggrantp-loggrant_ss)+RHOGRANT*(loggrant-loggrant_ss-loghatz+log(LAMBDA_X))+(1-RHOGRANT)*RHOGRANT_ETA*(   d*exp(logtheta) - MUD*THETABAR );

f36=exp(logBbackp)-exp(logtilnx)-exp(logRSTAR)*exp(logBback)/(exp(loghatz)*exp(logpi))+REMIT+exp(loggrantp);


% End of model %%%%%%%%%%%%%%%%%%%%%%%

f=[f0;f1;f2;f3;f4;f5;f6;f7;f8;f9;f9a;f9b;f9c;f10;f11;f12;f13;f14;f15;f16;f19;f20;f21;f22;f23;f24;f21b;f25;f26;f28;f29;f30;f31;f32;f33;f34;f35;f36;f37;f38;f39;f40];%
f=simplify(f);

f=subs(f,[hatz,tilu,tilc,l],exp([loghatz,logtilu,logtilc,logl]));
f=subs(f,[hatzp,tilup,tilcp,lp],exp([loghatzp,logtilup,logtilcp,loglp]));

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
  


y=[auxvar1,logtilv_tilvss,logtilu,logtilc,nl,logtillambda,logtilw,loguc,...
    lognegtilul,logmnom,logmden,logtilr,logtilq,logtilx,logtilk,logtily,...
    S,Sprime,...
    logtilg1,logpi,logpistar,aux2,logtilkg,logmc,logtilxg,logby,logRSTAR,logtilnx,auxvar1_2,logtilv_tilvss_2,logtilu_2,logtilv_2];
yp=[auxvar1p,logtilvp_tilvss,logtilup,logtilcp,nlp,logtillambdap,logtilwp,logucp,...
    lognegtilulp,logmnomp,logmdenp,logtilrp,logtilqp,logtilxp,logtilkp,logtilyp,...
    Sp,Sprimep,...
    logtilg1p,logpip,logpistarp,aux2p,logtilkgp,logmcp,logtilxgp,logbyp,logRSTARp,logtilnxp,auxvar1p_2,logtilvp_tilvss_2,logtilup_2,logtilvp_2];

symparams=[MUD RHOTHETA THETABAR SDV_THETA SDV_ZA ALPHA LAMBDA_A GAMMA ...
    NU PSI BETA DELTA PHI SCALEPARAM LOGTILUSS KAPPA LAMBDA_X...
    THETA_P CHI EPSILON RSS PISS GAMMA_PI RHOG GY GBAR logtauc_ss logBbackss_ss BY RHOTAUB LUMPSUM RHOTAU...
    SDV_XG DELTAG logtilxg_ss RHOXG logtilkg_ss ALPHAG XGY logtily_ss logtilkgstarback_ss RHOKG RHOTAUY RHOGB RHOGY RHOZA  ETA REMIT NXY RHOGRANT RHOGRANT_ETA GRANTY loggrant_ss LOSS HATZ_ss  LOGTILUSS_2 tilv_2_ss];%LOGTILWELFSS
    
