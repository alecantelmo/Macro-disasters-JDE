home=pwd;
cd ..\;
cd RBC_EZ

parameters;

KAPPA=9.5; %from Fernandez-Vaillaverde et al (2015)

SteadyState; % to get hatzss;
LAMBDA_X=hatzss;

eta_mat=[0,0,0,0;
         0,0,0,0; 
         0,0,0,0;
         0,0,0,0;
         1,0,0,0;
         0,SDV_THETA,0,0;
         0,0,SDV_ZA,0;
         0,0,0,0]; %(Alessandro)

cd(home);