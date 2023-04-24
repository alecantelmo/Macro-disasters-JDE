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

syms MUD RHOTHETA THETABAR SDV_THETA SDV_ZA ALPHA LAMBDA_A GAMMA NU PSI BETA DELTA PHI SCALEPARAM LOGTILUSS real
syms logg loggp  GBAR RHOG GY real % add government spending 
syms logtauc logtaucp logBback logBbackp real% add taxes and debt variables (Alessandro)
syms logtauc_ss logBbackss_ss BY RHOTAUB LUMPSUM real% add taxes and debt parameters (Alessandro)
syms logtaucback RHOTAU real


%Equation numbers refer to the online appendix of the paper (Alessandro, Feb'18)

% Exogenous state variables

%x2p is the second block of the vector x at time t+1, which contains the exogenous state variables  (Alessandro, Feb'18)
x2p=[dp;        %disaster variable
    logthetap;  %disaster size 
    zap;        %Gaussian component of TFP shock
    loggp];       

%Phi_fun is the lower block of the matrix h: it contains the expected
%value of future exogenous state variables as a function of current state
%variables (Alessandro, Feb'18)
Phi_fun=[MUD; %disaster probability (Alessandro, Feb'18)
    (1-RHOTHETA)*log(THETABAR)+RHOTHETA*logtheta; %rhs of logthetadp (Alessandro, Feb'18)
    0;
    (1-RHOG)*log(GBAR)+RHOG*logg];

% matrix of shocks (Alessandro, Feb'18)
eta_mat=[0,0,0,0;
         0,0,0,0;
         0,0,0,0;
         1,0,0,0;
         0,SDV_THETA,0,0;
         0,0,SDV_ZA,0;
         0,0,0,0]; %(Alessandro)

% other variables that depend only on the exogenous state vars
% equations page 11 (Alessandro, Feb'18)   
syms loghata loghatap loghatz loghatzp real 
loghata_=LAMBDA_A+za-(1-ALPHA)*d*exp(logtheta);                        % TFP growth  (Alessandro, Feb'18)        
loghatap_=LAMBDA_A+zap-(1-ALPHA)*dp*exp(logthetap);                    % TFP growth t+1 (Alessandro, Feb'18)   
loghatmu=0; loghatmup=0; % assume no investment specific shock
loghatz_=1/(1-ALPHA)*loghata+ALPHA/(1-ALPHA)*loghatmu;                 % TFP factor in fixed costs growth (Alessandro, Feb'18)  
loghatzp_=1/(1-ALPHA)*loghatap+ALPHA/(1-ALPHA)*loghatmup;              % TFP factor in fixed costs growth t+1 (Alessandro, Feb'18) 

% define auxiliary variable nl to bind l within [0,1]
syms nl nlp real

logl=nl-log(1+exp(nl));
loglp=nlp-log(1+exp(nlp));

%the next blockes define the components of equation 11 (Alessandro, Feb'18) 
syms u01p u02 real
u01p_=logtilvp_tilvss+loghatzp; %last term on rhs (Alessandro, Feb'18) 
u02_=exp(-auxvar1);             %auxvar1 is Vtilde(+1)/Vtildess (Alessandro, Feb'18)
f0=-1+exp(u01p)^(1-GAMMA)*u02; % auxiliary variable
% f1=-exp(logtilv_tilvss)^(1-PSI)+exp(logtilu-LOGTILUSS)^(1-PSI)*SCALEPARAM+BETA*auxvar1^((1-PSI)/(1-GAMMA));

syms u1 u2 u1p u2p real

u1_=exp(u2)*SCALEPARAM+BETA*exp(auxvar1*((1-PSI)/(1-GAMMA)));       % rhs of eq 11 (Alessandro, Feb'18)
u1p_=exp(u2p)*SCALEPARAM+BETA*exp(auxvar1p*((1-PSI)/(1-GAMMA)));    %rhs of eq 11 t+1 (Alessandro, Feb'18)

u2_=(logtilu-LOGTILUSS)*(1-PSI); %first term on rhs of eq 11 (Alessandro, Feb'18)
u2p_=(logtilup-LOGTILUSS)*(1-PSI);  %first term on rhs of eq 11 t+1 (Alessandro, Feb'18)

logtilv_tilvss_=1/(1-PSI)*log(u1);   %eq 11 (Alessandro, Feb'18)
logtilvp_tilvss_=1/(1-PSI)*log(u1p); %eq 11 t+1 (Alessandro, Feb'18)

% f2=-logtilu+logtilc+NU*log(1-l);
l_=exp(logl);                       % labor (Alessandro, Feb'18)
lp_=exp(loglp);                     % labor t+1  (Alessandro, Feb'18)

logtilu_=logtilc+NU*log(1-l);       % eq 12: utility function (Alessandro, Feb'18)
logtilup_=logtilcp+NU*log(1-lp);    % eq 12: utility function t+1 (Alessandro, Feb'18)

% f3=-loguc+NU*log(1-l);
loguc_=NU*log(1-l);                 % eq 13: marginal utility of consumption (Alessandro, Feb'18)
logucp_=NU*log(1-lp);               % eq 13: marginal utility of consumption t+1 (Alessandro, Feb'18)

% f4=-lognegtilul+log(NU)+logtilc+(NU-1)*log(1-l); %lognegtilul is the log of negative tilul, because tilul itself is negative
% f5=log(1-PSI)-PSI*logtilu+lognegtilul-logtillambda-logtilw; % note the use again of lognegtilul (and not logtilul)
% f6=log(1-PSI)-PSI*logtilu+loguc-logtillambda;
logtillambda_=log(1-PSI)-PSI*logtilu+loguc-log(1+exp(logtaucback));
logtillambdap_=log(1-PSI)-PSI*logtilup+logucp-log(1+exp(logtauc));

% use f3, f4, f5, f6 to solve for logtilc

logtilc_=logtilw+log(1-l)-log(NU);                  % eq. 15 labor supply schedule (Alessandro, Feb'18)
logtilcp_=logtilwp+log(1-lp)-log(NU);               % eq. 15 labor supply schedule t+1 (Alessandro, Feb'18)

% f7=-logmnom+log(BETA)+logtillambda-PSI*loghatz+(PSI-GAMMA)*logtilv_tilvss+(PSI-GAMMA)*loghatz; %nominator of logm
logmnom_=log(BETA)+logtillambda-PSI*loghatz+(PSI-GAMMA)*logtilv_tilvss+(PSI-GAMMA)*loghatz;       % eq 17: numerator of stocastic discount factor at time t (Alessandro, Feb'18)
logmnomp_=log(BETA)+logtillambdap-PSI*loghatzp+(PSI-GAMMA)*logtilvp_tilvss+(PSI-GAMMA)*loghatzp;  % eq 17: numerator of stocastic discount factor at time t+1 (Alessandro, Feb'18)

% f8=-logmden+logtillambda+(PSI-GAMMA)/(1-GAMMA)*log(auxvar1); %denominator of logmp
logmden_=logtillambda+(PSI-GAMMA)/(1-GAMMA)*auxvar1; %denominator of logmp %% eq 17: denominator of stocastic discount factor (Alessandro, Feb'18)

syms tilr tilrp tilq tilqp theta thetap real

tilrp_=exp(logtilrp);

theta_=exp(logtheta);
thetap_=exp(logthetap);

logtilq=0; logtilqp=0; % under no adj costs

syms u03 u04 real
logMp_=logmnomp-logmden;                    % eq 17 stocastic discount factor (Alessandro, Feb'18)
u03_=logMp-dp*thetap-loghatmup;             % eq 18 stocastic discount factor * disaster variables / investment shock growth (Alessandro, Feb'18)
u04_=1/tilq;                                % 1/q (Alessandro, Feb'18)
f9=exp(u03)*(tilrp+tilqp*(1-DELTA))*u04-1;  % eq 18 (Alessandro, Feb'18)

S=0; Sprime=0; Sp=0; Sprimep=0; % no adj costs

% f10=-1+exp(logtilq)*(1-S-S'*exp(logtilx-logtilxback+loghatz))+exp(logmnomp-logmden+logtilqp)*Sprimep*exp(2*(logtilxp-logtilx+loghatzp));

tilq_=exp(logtilq);
tilqp_=exp(logtilqp);

% f11=-logtilr+log(ALPHA)+loghata+loghatmu+(ALPHA-1)*(logtilkstarback-d*exp(logtheta))+(1-ALPHA)*logl;
logtilr_=log(ALPHA)+loghata+loghatmu+(ALPHA-1)*(logtilkstarback-d*theta)+(1-ALPHA)*logl;            % eq 36 marginal product of capital (Alessandro, Feb'18)
logtilrp_=log(ALPHA)+loghatap+loghatmup+(ALPHA-1)*(logtilkstarbackp-dp*thetap)+(1-ALPHA)*loglp;     % eq 36 marginal product of capital t+1 (Alessandro, Feb'18)


% f12=-logtilw+log(1-ALPHA)+loghata-loghatz+(ALPHA)*(logtilkstarback-d*exp(logtheta))+(-ALPHA)*logl;
logtilw_=log(1-ALPHA)+loghata-loghatz+(ALPHA)*(logtilkstarback-d*theta)+(-ALPHA)*logl;              % eq 37 marginal product of labor (Alessandro, Feb'18)
logtilwp_=log(1-ALPHA)+loghatap-loghatzp+(ALPHA)*(logtilkstarbackp-dp*thetap)+(-ALPHA)*loglp;       % eq 37 marginal product of labor t+1 (Alessandro, Feb'18)

% f13=-log(exp(logtily)+PHI)+loghata-loghatz+(ALPHA)*(logtilkstarback-d*exp(logtheta))+(1-ALPHA)*logl;
% f14=-exp(logtily)+exp(logtilc)+exp(logtilx);

% use f13 and f14 to solve for tilx
syms u4 tilx logtilx tily logtily auxg auxgp real
u4_=loghata-loghatz+(ALPHA)*(logtilkstarback-d*theta)+(1-ALPHA)*logl;           % eq 38 production function without fixed costs parameter (Alessandro, Feb'18)
tilc_=exp(logtilc);
% tilx_=exp(u4)-PHI-tilc;                                                         % eq 20 investment from resource constraint (Alessandro, Feb'18)
tilx_=exp(u4)-PHI-tilc-exp(logg);                                                         % eq 20 investment from resource constraint (Alessandro, Feb'18)
% tily_=tilc+tilx;                                                                % eq 20 resource constraint (Alessandro, Feb'18)
tily_=tilc+tilx+exp(logg);                                                                % eq 20 resource constraint (Alessandro, Feb'18)
logtilx_=log(tilx);
logtily_=log(tily);
auxg_=logg;
auxgp_=loggp;

syms u4p tilyp logtilyp tilxp tilcp real % same block as above at t+1 (Alessandro, Feb'18)
u4p_=loghatap-loghatzp+(ALPHA)*(logtilkstarbackp-dp*thetap)+(1-ALPHA)*loglp;    
tilcp_=exp(logtilcp); 
% tilxp_=exp(u4p)-PHI-tilcp;
tilxp_=exp(u4p)-PHI-tilcp-exp(loggp);
tilyp_=exp(u4p)-PHI;
logtilyp_=log(tilyp);


syms u05 real
u05_=exp(-logtilx);
f15=exp(logtilkstarbackp)*u05-(1-DELTA)*exp(logtilk)*u05-(1-S);             % eq 21 law of motion of capital (Alessandro, Feb'18)

% f16=-logtilk+logtilkstarback-d*exp(logtheta)-loghatz-loghatmu;
logtilk_=logtilkstarback-d*theta-loghatz-loghatmu;                          % eq 22 capital net of disaster (Alessandro, Feb'18)


%asset prices block (Alessandro, Feb'18)
syms logtildiv logtildivp w wp tildiv tildivp real
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

%end of asset price block (Alessandro, Feb'18)

% Government budget
f25=-exp(logBbackp)+exp(logR)*exp(logBback)/(exp(loghatz))+exp(logg)-exp(logtaucback)*exp(logtilc)-LUMPSUM;
f26=-(logtauc-logtauc_ss)+RHOTAU*(logtaucback-logtauc_ss)+(1-RHOTAU)*RHOTAUB*(logBback-logBbackss_ss);

% Euler Equation
syms logR logRp real
f20=-1+exp(logmnomp-logmden+logR);


f=[f0;f9;f15;f17;f18;f25;f26;f20];
f=simplify(f);

x=[logtilkstarback;
    logBback;
    logtaucback;
    d;
    logtheta;
    za;
    logg];  %Alessandro

xp=[logtilkstarbackp;
    logBbackp;
    logtauc;
    dp;
    logthetap;
    zap;
    loggp];  %Alessandro

y=[auxvar1,nl,logtilqe,logqf,logR];      %vector of control variables at time t (Alessandro, Feb'18)
yp=[auxvar1p,nlp,logtilqep,logqfplogR]; %vector of control variables at time t+1 (Alessandro, Feb'18)


subsvars=[loghata;loghatap;loghatz;loghatzp;u1;u1p;u2;u2p;logtilv_tilvss;logtilvp_tilvss;l;lp;logtilu;logtilup;loguc;logucp;...
    logtillambda;logtillambdap;logtilc;logtilcp;logmnom;logmnomp;logmden;...
    theta;thetap;tilrp;tilq;tilqp;logtilr;logtilrp;logtilw;logtilwp;...
    u4;tilc;tilcp;tilx;tilxp;logtilx;logtilk;u01p;u02;u03;u04;u05;logMp;u17;u4p;tily;logtily;tilyp;logtilyp;...
    w;wp;tildiv;tildivp;auxg;auxgp];

subsfuns=sym(zeros(size(subsvars)));
for i=1:length(subsfuns)
    subsfuns(i)=eval([char(subsvars(i)) '_']);
end

symparams=[MUD RHOTHETA THETABAR SDV_THETA SDV_ZA ALPHA LAMBDA_A GAMMA NU PSI BETA DELTA PHI SCALEPARAM LOGTILUSS RHOG GY GBAR logtauc_ss logBbackss_ss BY RHOTAUB LUMPSUM RHOTAU]; %RHOG GY added by Alessandro