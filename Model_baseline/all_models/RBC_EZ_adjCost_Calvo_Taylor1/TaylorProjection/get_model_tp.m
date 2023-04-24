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


% Exogenous state variables

x2p=[dp;
    logthetap;
    zap];

Phi_fun=[MUD;
    (1-RHOTHETA)*log(THETABAR)+RHOTHETA*logtheta;
    0];

eta_mat=[0,0,0;
         0,0,0;
         0,0,0;
         0,0,0;
         0,0,0;
         1,0,0;
         0,SDV_THETA,0;
         0,0,SDV_ZA];
     
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

syms u01p u02 real
u01p_=logtilvp_tilvss+loghatzp;
u02_=exp(-auxvar1);
f0=-1+exp(u01p)^(1-GAMMA)*u02; % auxiliary variable
% f1=-exp(logtilv_tilvss)^(1-PSI)+exp(logtilu-LOGTILUSS)^(1-PSI)*SCALEPARAM+BETA*auxvar1^((1-PSI)/(1-GAMMA));

syms u1 u2 u1p u2p real

u1_=exp(u2)*SCALEPARAM+BETA*exp(auxvar1*((1-PSI)/(1-GAMMA)));
u1p_=exp(u2p)*SCALEPARAM+BETA*exp(auxvar1p*((1-PSI)/(1-GAMMA)));

u2_=(logtilu-LOGTILUSS)*(1-PSI);
u2p_=(logtilup-LOGTILUSS)*(1-PSI);

logtilv_tilvss_=1/(1-PSI)*log(u1);
logtilvp_tilvss_=1/(1-PSI)*log(u1p);

% f2=-logtilu+logtilc+NU*log(1-l);
l_=exp(logl);
lp_=exp(loglp);

logtilu_=logtilc+NU*log(1-l);
logtilup_=logtilcp+NU*log(1-lp);

% f3=-loguc+NU*log(1-l);
loguc_=NU*log(1-l);
logucp_=NU*log(1-lp);

% f4=-lognegtilul+log(NU)+logtilc+(NU-1)*log(1-l); %lognegtilul is the log of negative tilul, because tilul itself is negative
% f5=log(1-PSI)-PSI*logtilu+lognegtilul-logtillambda-logtilw; % note the use again of lognegtilul (and not logtilul)
% f6=log(1-PSI)-PSI*logtilu+loguc-logtillambda;
logtillambda_=log(1-PSI)-PSI*logtilu+loguc;
logtillambdap_=log(1-PSI)-PSI*logtilup+logucp;

% use f3, f4, f5, f6 to solve for logtilc

logtilc_=logtilw+log(1-l)-log(NU);
logtilcp_=logtilwp+log(1-lp)-log(NU);

% f7=-logmnom+log(BETA)+logtillambda-PSI*loghatz+(PSI-GAMMA)*logtilv_tilvss+(PSI-GAMMA)*loghatz; %nominator of logm
logmnom_=log(BETA)+logtillambda-PSI*loghatz+(PSI-GAMMA)*logtilv_tilvss+(PSI-GAMMA)*loghatz;
logmnomp_=log(BETA)+logtillambdap-PSI*loghatzp+(PSI-GAMMA)*logtilvp_tilvss+(PSI-GAMMA)*loghatzp;

% f8=-logmden+logtillambda+(PSI-GAMMA)/(1-GAMMA)*log(auxvar1); %denominator of logmp
logmden_=logtillambda+(PSI-GAMMA)/(1-GAMMA)*auxvar1; %denominator of logmp

syms tilr tilrp tilq tilqp theta thetap real

tilrp_=exp(logtilrp);

theta_=exp(logtheta);
thetap_=exp(logthetap);

% logtilq=0; logtilqp=0; % under no adj costs

syms u03 u04 real
logMp_=logmnomp-logmden;
u03_=logMp-dp*thetap-loghatmup;
u04_=1/tilq;
f9=exp(u03)*(tilrp+tilqp*(1-DELTA))*u04-1;

% S=0; Sprime=0; Sp=0; Sprimep=0; % no adj costs
syms S Sprime logtilxback Sp Sprimep logtilxbackp real
S_=KAPPA/2*(exp(logtilx-logtilxback+loghatz)-LAMBDA_X)^2;
Sprime_=KAPPA*(exp(logtilx-logtilxback+loghatz)-LAMBDA_X);
Sprimep_=KAPPA*(exp(logtilxp-logtilxbackp+loghatzp)-LAMBDA_X);

f9c=-logtilx+logtilxbackp;

f10=-1+exp(logtilq)*(1-S-Sprime*exp(logtilx-logtilxback+loghatz))+...
    exp(logmnomp-logmden+logtilqp)*Sprimep*exp(2*(logtilxp-logtilx+loghatzp));

tilq_=exp(logtilq);
tilqp_=exp(logtilqp);

% f11=-logtilr+log(ALPHA)+loghata+loghatmu+(ALPHA-1)*(logtilkstarback-d*exp(logtheta))+(1-ALPHA)*logl;
% logtilr_=log(ALPHA)+loghata+loghatmu+(ALPHA-1)*(logtilkstarback-d*theta)+(1-ALPHA)*logl;
% logtilrp_=log(ALPHA)+loghatap+loghatmup+(ALPHA-1)*(logtilkstarbackp-dp*thetap)+(1-ALPHA)*loglp;


% f12=-logtilw+log(1-ALPHA)+loghata-loghatz+(ALPHA)*(logtilkstarback-d*exp(logtheta))+(-ALPHA)*logl;
% logtilw_=log(1-ALPHA)+loghata-loghatz+(ALPHA)*(logtilkstarback-d*theta)+(-ALPHA)*logl;
% logtilwp_=log(1-ALPHA)+loghatap-loghatzp+(ALPHA)*(logtilkstarbackp-dp*thetap)+(-ALPHA)*loglp;

% f14=-1+[exp(logtilc)+exp(logtilx)]*exp(-logtily);

% % use f13 and f14 to solve for tilx
syms u4 tilx logtilx tily logtily real
% u4_=loghata-loghatz+(ALPHA)*(logtilkstarback-d*theta)+(1-ALPHA)*logl;
tilc_=exp(logtilc);
% tilx_=(exp(u4)-PHI)-tilc;
tilx_=exp(logtilx);
tily_=tilc+tilx;
% logtilx_=log(tilx);
logtily_=log(tily);
% 
syms u4p tilyp logtilyp tilxp tilcp real
% u4p_=loghatap-loghatzp+(ALPHA)*(logtilkstarbackp-dp*thetap)+(1-ALPHA)*loglp;
tilcp_=exp(logtilcp);
% tilxp_=exp(u4p)-PHI-tilcp;
tilxp_=exp(logtilxp);
% tilyp_=exp(u4p)-PHI;
tilyp_=tilcp+tilxp;
% logtilxp_=log(tilxp);
logtilyp_=log(tilyp);


syms u05 real
u05_=exp(-logtilx);
f15=exp(logtilkstarbackp)*u05-(1-DELTA)*exp(logtilk)*u05-(1-S);

% f16=-logtilk+logtilkstarback-d*exp(logtheta)-loghatz-loghatmu;
logtilk_=logtilkstarback-d*theta-loghatz-loghatmu;
logtilkp_=logtilkstarbackp-dp*thetap-loghatzp-loghatmup;

syms w wp real
syms tildiv tildivp logtilqe logtilqep real

w_=exp(logtilw);
wp_=exp(logtilwp);

tildiv_=tily-w*l-tilx;
tildivp_=tilyp-wp*lp-tilxp;
% logtildiv_=log(tildiv); % dividends
% logtildivp_=log(tildivp);

u17_=logMp+loghatzp;
f17=-1+exp(u17)*(tildivp+exp(logtilqep))*exp(-logtilqe);

syms logqf logqfp real
f18=-1+exp(logMp)*exp(-logqf);

% Calvo Pricing %%%%%%%%%%%%%
syms logmc logtilg1 logtilg1p logtilg2 logtilg2p logpi logpip logpistar logpistarp logpiback logpibackp real
syms logvpback logvpbackp real
syms THETA_P CHI EPSILON real

% FOC

logmc_=(ALPHA-1)*log(1-ALPHA)-ALPHA*log(ALPHA)+(1-ALPHA)*logtilw+ALPHA*logtilr;
f19=-1+[exp(logmc+logtily)+THETA_P*exp(logmnomp-logmden-EPSILON*(CHI*logpi-logpip)+logtilg1p+loghatzp)]*exp(-logtilg1);

logtilg2_=log(EPSILON)+logtilg1-log(EPSILON-1);
logtilg2p_=log(EPSILON)+logtilg1p-log(EPSILON-1);

% f21=-1+THETA_P*exp((1-EPSILON)*(CHI*logpiback-logpi))+(1-THETA_P)*exp((1-EPSILON)*logpistar); % use this to substitute out logpistar
% PROBLEM: logpistar can become a complex number if
% (1-THETA_P*exp((1-EPSILON)*(CHI*logpiback-logpi)))<0. See technical tips
% for the solution used here.
% logpistar_=[log(1-THETA_P*exp((1-EPSILON)*(CHI*logpiback-logpi)))-log(1-THETA_P)]/(1-EPSILON);
% logpistarp_=[log(1-THETA_P*exp((1-EPSILON)*(CHI*logpibackp-logpip)))-log(1-THETA_P)]/(1-EPSILON);

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
% f11=-(logtilk-logl)+log(ALPHA)-log(1-ALPHA)+logtilw-logtilr; %use this to
% substitute out logtilr (need to keep logtilw in)
logtilr_=-(logtilk-logl)+log(ALPHA)-log(1-ALPHA)+logtilw;
logtilrp_=-(logtilkp-loglp)+log(ALPHA)-log(1-ALPHA)+logtilwp;

% Aggregation 
syms u13a real
u13a_=exp(logtily+logvpbackp)+PHI;
f13=-log(u13a)+loghata-loghatz+(ALPHA)*(logtilkstarback-d*exp(logtheta))+(1-ALPHA)*logl;

% End of Calvo Pricing %%%%%%%%%%%%%%%

% Taylor rule with backwards output
syms logR logRp RSS PISS GAMMA_PI real
syms GAMMA_Y LAMBDA_Y logtilyback logtilybackp real

logR_=log(RSS)+GAMMA_PI*(logpi-log(PISS))+GAMMA_Y*(logtily-logtilyback+loghatz-LAMBDA_Y);
f24=-1+exp(logmnomp-logmden+logR-logpip);
f25=-logtilybackp+logtily;


f=[f0;f9;f13;f15;f17;f18;f9c;f10;f19;f20;f22;f23;f24;f25];
f=simplify(f);


x=[logtilkstarback;
    logtilxback;
    logpiback;
    logvpback;
    logtilyback;
    d;
    logtheta;
    za];

xp=[logtilkstarbackp;
    logtilxbackp;
    logpibackp;
    logvpbackp;
    logtilybackp;
    dp;
    logthetap;
    zap];

y=[auxvar1,nl,logtilw,logtilx,logtilqe,logqf,logtilq,logtilg1,aux2];
yp=[auxvar1p,nlp,logtilwp,logtilxp,logtilqep,logqfp,logtilqp,logtilg1p,aux2p];


subsvars=[loghata;loghatap;loghatz;loghatzp;u1;u1p;u2;u2p;logtilv_tilvss;logtilvp_tilvss;l;lp;logtilu;logtilup;loguc;logucp;...
    logtillambda;logtillambdap;logtilc;logtilcp;logmnom;logmnomp;logmden;...
    theta;thetap;tilrp;tilq;tilqp;logtilr;logtilrp;...
    tilc;tilcp;tilx;tilxp;logtilk;logtilkp;u01p;u02;u03;u04;u05;logMp;u17;tily;logtily;tilyp;logtilyp;...
    w;wp;tildiv;tildivp;S;Sprime;Sprimep;logmc;logtilg2;logtilg2p;logpistar;logpistarp;...
    u13a;logR;pi_power;pi_powerp;logpi;logpip];


subsfuns=sym(zeros(size(subsvars)));
for i=1:length(subsfuns)
    subsfuns(i)=eval([char(subsvars(i)) '_']);
end

symparams=[MUD RHOTHETA THETABAR SDV_THETA SDV_ZA ALPHA LAMBDA_A GAMMA ...
    NU PSI BETA DELTA PHI SCALEPARAM LOGTILUSS KAPPA LAMBDA_X...
    THETA_P CHI EPSILON RSS PISS GAMMA_PI GAMMA_Y LAMBDA_Y];
