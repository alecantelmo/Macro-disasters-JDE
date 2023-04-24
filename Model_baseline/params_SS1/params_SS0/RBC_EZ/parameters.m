% Parameter values

NU=2.33; % Leisure preference
GAMMA=3.8;  % Risk aversion
hatPSI=.5; % inverse IES
PSI=1-(1-hatPSI)/(1+NU); % Epstein-Zin
% PSI=GAMMA; hatPSI=1-(1-PSI)*(1+NU); % Expected Utility Assumption

% From Fernandez-Villaverde et al (Journal of Econometrics 2015)
LAMBDA_A=.0028;
LAMBDA_MU=0;
BETA=.99; 
ALPHA=0.21;
DELTA=.025;
PHI=0;

SDV_ZA=.01 ;  % from Gourio (2012)

% Disaster Probability: from the Working paper version of Gourio (2010),
% which is based on Barro (2006)
prob_disaster=1-(1-.017)^.25; % disaster probability (quareterly)

MUD=prob_disaster;
THETABAR=-log(1-.00001); % mean disaster impact
SDV_THETA=.00001; % standrad deviation of disaster impact.
RHOTHETA=0.9; % persistence of disaster effect

%New parameters (Alessandro)
GBAR=0.1125;  %This needs to be calibrated to initialize the steady state
RHOG=0.80;
GY=0.20;

% Define eta as in Schmitt-Grohe and Uribe (2004). 
eta_mat=[0,0,0,0;
         1,0,0,0;
         0,SDV_THETA,0,0;
         0,0,SDV_ZA,0;
         0,0,0,0]; %(Alessandro)

