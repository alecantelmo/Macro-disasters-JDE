home2=pwd;

% Select calibration of natural disasters
do_disaster_prone=1;
do_non_disaster_prone=0;
do_cc_prob=0;
do_cc_damages=0;
do_cc_both=0;
do_calibration_Haiti=0;
do_no_disaster_scenario=0;
do_disaster_prone_alternative=0;
do_non_disaster_prone_alternative=0;

% Parameter values

NU=1.59; % Leisure preference
GAMMA=3.8; % Risk aversion
hatPSI=.5; % inverse IES
PSI=1-(1-hatPSI)/(1+NU); % Epstein-Zin

LAMBDA_A=.0035;
LAMBDA_MU=0;


ALPHA=0.32;
DELTA=.025;
PHI=0;

RHOZA=0.5;
SDV_ZA=0.025 ;

% Disasters calibration
if do_disaster_prone==1
    prob_disaster=1-(1-.162)^.25; 
    MUD=prob_disaster;
    THETABAR=-log(1-.0665); 
    SDV_THETA=0.127;
    RSS= 1.0184;
elseif do_non_disaster_prone==1
    prob_disaster=1-(1-.0028)^.25; 
    MUD=prob_disaster;
    THETABAR=-log(1-.0052);
    SDV_THETA=0.017;
    RSS= 1.0213;
elseif do_cc_prob==1
    prob_disaster=1-(1-.219)^.25; 
    MUD=prob_disaster;
    THETABAR=-log(1-.0665); 
    SDV_THETA=0.127;
    RSS= 1.0184;
elseif do_cc_damages==1
    prob_disaster=1-(1-.162)^.25; 
    MUD=prob_disaster;
    THETABAR=-log(1-.121);
    SDV_THETA=0.127;
    RSS= 1.0184;
elseif do_cc_both==1
    prob_disaster=1-(1-.219)^.25; 
    MUD=prob_disaster;
    THETABAR=-log(1-.121);
    SDV_THETA=0.127;
    RSS= 1.0184;
elseif do_calibration_Haiti==1
    prob_disaster=1-(1-.00001)^.25;
    MUD=prob_disaster;
    THETABAR=-log(1-.222);
    SDV_THETA=0.0127;
    RSS= 1.0184;
elseif do_no_disaster_scenario==1
    prob_disaster=1-(1-.00000001)^.25;
    MUD=prob_disaster;
    THETABAR=-log(1-.00000001);
    SDV_THETA=0;
    RSS= 1.0213;
elseif do_disaster_prone_alternative==1
    prob_disaster=1-(1-.103)^.25; 
    MUD=prob_disaster;
    THETABAR=-log(1-.107); 
    SDV_THETA=0.127;
    RSS= 1.0184;
elseif do_non_disaster_prone_alternative==1
    prob_disaster=1-(1-.093)^.25; 
    MUD=prob_disaster;
    THETABAR=-log(1-.027);
    SDV_THETA=0.127;
    RSS= 1.0184;

end

RHOTHETA=0.9; % persistence of disaster impact




RHOG=0.99;
RHOGB=0;
RHOGY=0;

GY=0.16;
XGY=0.07;

BY=0.580;
RHOTAUB=0.225;
RHOTAU=0.90;
TAUCSS=0.21;
RHOTAUY=0;

KAPPA=12; 


THETA_P=0;
CHI=0;
EPSILON=100; 

PISS=1;
GAMMA_PI=1.3;



DELTAG=0.0125;
ALPHAG=0.22;
RHOXG=0.95;
RHOKG=1.5;
SDV_XG=0;

ETA=0;

RHOGRANT=0;
RHOGRANT_ETA=0;
GRANTY=0.00001;

REMIT=0;
NXY=0.12;

LOSS=0;



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
     
cd(home2);