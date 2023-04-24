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

syms LAMBDA_MU SDV_ZMU loghatmu loghatmup real
syms m mp SDV_M xi xip RHOXI SDV_XI real

% Exogenous state variables

x2p=[dp;
    logthetap;
    zap;
    loghatmup;
    mp;
    xip];

Phi_fun=[MUD;
    (1-RHOTHETA)*log(THETABAR)+RHOTHETA*logtheta;
    0;
    LAMBDA_MU;
    0;
    RHOXI*xi];

eta_mat=[0,0,0,0,0,0;
         0,0,0,0,0,0;
         0,0,0,0,0,0;
         0,0,0,0,0,0;
         0,0,0,0,0,0;
         0,0,0,0,0,0;
         1,0,0,0,0,0;
         0,SDV_THETA,0,0,0,0;
         0,0,SDV_ZA,0,0,0;
         0,0,0,SDV_ZMU,0,0;
         0,0,0,0,SDV_M,0;
         0,0,0,0,0,SDV_XI];

% other variables that depend only on the exogenous state vars
loghata=LAMBDA_A+za-(1-ALPHA)*d*exp(logtheta);
loghatap=LAMBDA_A+zap-(1-ALPHA)*dp*exp(logthetap);
% loghatmu=0; loghatmup=0; % assume no investment specific shock
loghatz=1/(1-ALPHA)*loghata+ALPHA/(1-ALPHA)*loghatmu;
loghatzp=1/(1-ALPHA)*loghatap+ALPHA/(1-ALPHA)*loghatmup;

% define auxiliary variable nl to bind l within [0,1]
syms nl nlp real
l=exp(nl)/(1+exp(nl));
lp=exp(nlp)/(1+exp(nlp));
logl=nl-log(1+exp(nl));
loglp=nlp-log(1+exp(nlp));


f0=-1+exp(-auxvar1)*exp(logtilvp_tilvss+loghatzp)^(1-GAMMA); % auxiliary variable - VERY IMPORTANT to do it in log
f1=[-exp(logtilv_tilvss)^(1-PSI)+exp(logtilu-LOGTILUSS)^(1-PSI)*SCALEPARAM+BETA*exp(auxvar1*(1-PSI)/(1-GAMMA))]*exp(-logtilv_tilvss)^(1-PSI);

f2=-logtilu+logtilc+NU*log(1-l)+xi;

f3=-loguc+NU*log(1-l)+xi;
f4=-lognegtilul+log(NU)+logtilc+(NU-1)*log(1-l)+xi; %lognegtilul is the log of negative tilul, because tilul itself is negative
f5=log(1-PSI)-PSI*logtilu+lognegtilul-logtillambda-logtilw; % note the use again of lognegtilul (and not logtilul)
f6=log(1-PSI)-PSI*logtilu+loguc-logtillambda;

f7=-logmnom+log(BETA)+logtillambda-PSI*loghatz+(PSI-GAMMA)*logtilv_tilvss+(PSI-GAMMA)*loghatz; %nominator of logm


f8=-logmden+logtillambda+(PSI-GAMMA)/(1-GAMMA)*auxvar1; %denominator of logmp

f9=exp(logmnomp-logmden-dp*exp(logthetap)-loghatmup)*(exp(logtilrp)+exp(logtilqp)*(1-DELTA))*exp(-logtilq)-1;

% S=0; Sprime=0; Sp=0; Sprimep=0; % no adj costs
syms S Sprime logtilxback Sp Sprimep logtilxbackp real
f9a=-S+KAPPA/2*(exp(logtilx-logtilxback+loghatz)-LAMBDA_X)^2;
f9b=-Sprime+KAPPA*(exp(logtilx-logtilxback+loghatz)-LAMBDA_X);
f9c=-logtilx+logtilxbackp;


f10=-1+exp(logtilq)*(1-S-Sprime*exp(logtilx-logtilxback+loghatz))+...
    exp(logmnomp-logmden+logtilqp)*Sprimep*exp(2*(logtilxp-logtilx+loghatzp));

% f11=-logtilr+log(ALPHA)+loghata+loghatmu+(ALPHA-1)*(logtilkstarback-d*exp(logtheta))+(1-ALPHA)*logl;
% f12=-logtilw+log(1-ALPHA)+loghata-loghatz+(ALPHA)*(logtilkstarback-d*exp(logtheta))+(-ALPHA)*logl;

% f13=-log(exp(logtily)+PHI)+loghata-loghatz+(ALPHA)*(logtilkstarback-d*exp(logtheta))+(1-ALPHA)*logl;

f14=[-exp(logtily)+exp(logtilc)+exp(logtilx)]*exp(-logtilx);

f15=exp(logtilkstarbackp)*exp(-logtilx)-(1-DELTA)*exp(logtilk)*exp(-logtilx)-(1-S);
f16=-logtilk+logtilkstarback-d*exp(logtheta)-loghatz-loghatmu;

syms tildiv tildivp logtilqe logtilqep real

f17=-1+exp(-logtilqe)*exp(logmnomp-logmden+loghatzp)*(tildivp+exp(logtilqep));
f17a=[-tildiv+exp(logtily)-exp(logtilw+logl)-exp(logtilx)]*exp(-logtilx); % dividends equal total output excluding labour income and investments


syms logqf logqfp real
f18=-1+exp(-logqf)*exp(logmnomp-logmden);

% Calvo Pricing %%%%%%%%%%%%%
syms logmc logtilg1 logtilg1p logtilg2 logtilg2p logpi logpip logpistar logpistarp logpiback logpibackp real
syms logvpback logvpbackp real
syms THETA_P CHI EPSILON real

% FOC

logmc=(ALPHA-1)*log(1-ALPHA)-ALPHA*log(ALPHA)+(1-ALPHA)*logtilw+ALPHA*logtilr;
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

% Optimal factor composition
f11=-(logtilk-logl)+log(ALPHA)-log(1-ALPHA)+logtilw-logtilr;

% Aggregation 
f13=-log(exp(logtily+logvpbackp)+PHI)+loghata-loghatz+(ALPHA)*(logtilkstarback-d*exp(logtheta))+(1-ALPHA)*logl;

% End of Calvo Pricing %%%%%%%%%%%%%%%

% Taylor rule with backwards output and interest rate smoothing
syms logR logRp RSS PISS GAMMA_PI real
syms GAMMA_Y LAMBDA_Y logtilyback logtilybackp real
syms logRback logRbackp GAMMA_R real
logR=log(RSS)+GAMMA_R*(logRback-log(RSS))+(1-GAMMA_R)*[GAMMA_PI*(logpi-log(PISS))+GAMMA_Y*(logtily-logtilyback+loghatz-LAMBDA_Y)]+m;
f24=-1+exp(logmnomp-logmden+logR-logpip);
f25=-logtilybackp+logtily;
f26=-logRbackp+logR;

% End of model %%%%%%%%%%%%%%%%%%%%%%%

f=[f0;f1;f2;f3;f4;f5;f6;f7;f8;f9;f9a;f9b;f9c;f10;f11;f13;f14;f15;f16;f17;f17a;f18;f19;f20;f21;f22;f23;f24;f25;f26;f21b];
f=simplify(f);

f=subs(f,[hatz,tilu,tilc,l],exp([loghatz,logtilu,logtilc,logl]));
f=subs(f,[hatzp,tilup,tilcp,lp],exp([loghatzp,logtilup,logtilcp,loglp]));

x=[logtilkstarback;
    logtilxback;
    logpiback;
    logvpback;
    logtilyback;
    logRback;
    d;
    logtheta;
    za;
    loghatmu;
    m;
    xi];

xp=[logtilkstarbackp;
    logtilxbackp;
    logpibackp;
    logvpbackp;
    logtilybackp;
    logRbackp;
    dp;
    logthetap;
    zap;
    loghatmup;
    mp;
    xip];

y=[auxvar1,logtilv_tilvss,logtilu,logtilc,nl,logtillambda,logtilw,loguc,...
    lognegtilul,logmnom,logmden,logtilr,logtilq,logtilx,logtilk,logtily,...
    logtilqe,logqf,tildiv,S,Sprime,...
    logtilg1,logpi,logpistar,aux2];
yp=[auxvar1p,logtilvp_tilvss,logtilup,logtilcp,nlp,logtillambdap,logtilwp,logucp,...
    lognegtilulp,logmnomp,logmdenp,logtilrp,logtilqp,logtilxp,logtilkp,logtilyp,...
    logtilqep,logqfp,tildivp,Sp,Sprimep,...
    logtilg1p,logpip,logpistarp,aux2p];

symparams=[MUD RHOTHETA THETABAR SDV_THETA SDV_ZA ALPHA LAMBDA_A GAMMA ...
    NU PSI BETA DELTA PHI SCALEPARAM LOGTILUSS KAPPA LAMBDA_X...
    THETA_P CHI EPSILON RSS PISS GAMMA_PI GAMMA_Y LAMBDA_Y GAMMA_R LAMBDA_MU SDV_ZMU SDV_M RHOXI];

