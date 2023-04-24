tset 
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
syms logg loggp g gp GBAR RHOG GY real % add government spending 

syms logtauc logtaucp logBback logBbackp real% add taxes and debt variables (Alessandro)
syms logtauc_ss logBbackss_ss BY RHOTAUB LUMPSUM real% add taxes and debt parameters (Alessandro)
syms logtaucback RHOTAU real

% Exogenous state variables

x2p=[dp;
    logthetap;
    zap;
    loggp]; %Alessandro


Phi_fun=[MUD;
    (1-RHOTHETA)*log(THETABAR)+RHOTHETA*logtheta;
    0;
    (1-RHOG)*log(GBAR)+RHOG*logg]; %(Alessandro)


eta_mat=[0,0,0,0;
         0,0,0,0; 
         0,0,0,0;
         0,0,0,0;
         1,0,0,0;
         0,SDV_THETA,0,0;
         0,0,SDV_ZA,0;
         0,0,0,0]; %(Alessandro)


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


f0=-1+exp(-auxvar1)*exp(logtilvp_tilvss+loghatzp)^(1-GAMMA); % auxiliary variable - VERY IMPORTANT to do it in log
f1=[-exp(logtilv_tilvss)^(1-PSI)+exp(logtilu-LOGTILUSS)^(1-PSI)*SCALEPARAM+BETA*exp(auxvar1*(1-PSI)/(1-GAMMA))]*exp(-logtilv_tilvss)^(1-PSI);

f2=-logtilu+logtilc+NU*log(1-l);

f3=-loguc+NU*log(1-l);
f4=-lognegtilul+log(NU)+logtilc+(NU-1)*log(1-l); %lognegtilul is the log of negative tilul, because tilul itself is negative
f5=log(1-PSI)-PSI*logtilu+lognegtilul-logtillambda-logtilw; % note the use again of lognegtilul (and not logtilul)
f6=log(1-PSI)-PSI*logtilu+loguc-logtillambda;

f7=-logmnom+log(BETA)+logtillambda-PSI*loghatz+(PSI-GAMMA)*logtilv_tilvss+(PSI-GAMMA)*loghatz; %nominator of logm


f8=-logmden+logtillambda+(PSI-GAMMA)/(1-GAMMA)*auxvar1; %denominator of logmp

f9=exp(logmnomp-logmden-dp*exp(logthetap)-loghatmup)*(exp(logtilrp)+exp(logtilqp)*(1-DELTA))*exp(-logtilq)-1;


syms S Sprime logtilxback Sp Sprimep logtilxbackp real
f9a=-S+KAPPA/2*(exp(logtilx-logtilxback+loghatz)-LAMBDA_X)^2;
f9b=-Sprime+KAPPA*(exp(logtilx-logtilxback+loghatz)-LAMBDA_X);
f9c=-logtilx+logtilxbackp;


f10=-1+exp(logtilq)*(1-S-Sprime*exp(logtilx-logtilxback+loghatz))+...
    exp(logmnomp-logmden+logtilqp)*Sprimep*exp(2*(logtilxp-logtilx+loghatzp));

f11=-logtilr+log(ALPHA)+loghata+loghatmu+(ALPHA-1)*(logtilkstarback-d*exp(logtheta))+(1-ALPHA)*logl;
f12=-logtilw+log(1-ALPHA)+loghata-loghatz+(ALPHA)*(logtilkstarback-d*exp(logtheta))+(-ALPHA)*logl;

f13=-log(exp(logtily)+PHI)+loghata-loghatz+(ALPHA)*(logtilkstarback-d*exp(logtheta))+(1-ALPHA)*logl;

% f14=[-exp(logtily)+exp(logtilc)+exp(logtilx)]*exp(-logtilx);
f14=[-exp(logtily)+exp(logtilc)+exp(logtilx)+exp(logg)]*exp(-logtilx); %Alessandro added gov. spending

f15=exp(logtilkstarbackp)*exp(-logtilx)-(1-DELTA)*exp(logtilk)*exp(-logtilx)-(1-S);
f16=-logtilk+logtilkstarback-d*exp(logtheta)-loghatz-loghatmu;

syms tildiv tildivp logtilqe logtilqep real

f17=-1+exp(-logtilqe)*exp(logmnomp-logmden+loghatzp)*(tildivp+exp(logtilqep));
f17a=[-tildiv+exp(logtily)-exp(logtilw+logl)-exp(logtilx)]*exp(-logtilx); % dividends equal total output excluding labour income and investments


syms logqf logqfp real
f18=-1+exp(-logqf)*exp(logmnomp-logmden);

% Government budget
f25=-exp(logBbackp)+exp(logR)*exp(logBback)/exp(loghatz)+exp(logg)-exp(logtaucback)*exp(logtilc)-LUMPSUM;
f26=-(logtauc-logtauc_ss)+RHOTAU*(logtaucback-logtauc_ss)+(1-RHOTAU)*RHOTAUB*(logBback-logBbackss_ss);

%Euler Equation
syms logR logRp real
f20=-1+exp(logmnomp-logmden+logR);


f=[f0;f1;f2;f3;f4;f5;f6;f7;f8;f9;f9a;f9b;f9c;f10;f11;f12;f13;f14;f15;f16;f17;f17a;f18;f25;f26;f20];
f=simplify(f);

f=subs(f,[hatz,tilu,tilc,l],exp([loghatz,logtilu,logtilc,logl]));
f=subs(f,[hatzp,tilup,tilcp,lp],exp([loghatzp,logtilup,logtilcp,loglp]));

x=[logtilkstarback;
    logtilxback;
    logBback;
    logtaucback;
    d;
    logtheta;
    za;
    logg];  %Alessandro

xp=[logtilkstarbackp;
    logtilxbackp;
    logBbackp;
    logtauc;
    dp;
    logthetap;
    zap;
    loggp];  %Alessandro

y=[auxvar1,logtilv_tilvss,logtilu,logtilc,nl,logtillambda,logtilw,loguc,...
    lognegtilul,logmnom,logmden,logtilr,logtilq,logtilx,logtilk,logtily,...
    logtilqe,logqf,tildiv,S,Sprime,logR];
yp=[auxvar1p,logtilvp_tilvss,logtilup,logtilcp,nlp,logtillambdap,logtilwp,logucp,...
    lognegtilulp,logmnomp,logmdenp,logtilrp,logtilqp,logtilxp,logtilkp,logtilyp,...
    logtilqep,logqfp,tildivp,Sp,Sprimep,logRp];

symparams=[MUD RHOTHETA THETABAR SDV_THETA SDV_ZA ALPHA LAMBDA_A GAMMA ...
    NU PSI BETA DELTA PHI SCALEPARAM LOGTILUSS KAPPA LAMBDA_X RHOG GY GBAR logtauc_ss logBbackss_ss BY RHOTAUB LUMPSUM RHOTAU];
