home3=pwd;

cd ..\;
cd RBC_EZ_adjCost_Calvo

parameters;

GAMMA_Y=exp(-1.4034); %from Fernandez-Vaillaverde et al (2015)
LAMBDA_Y=loghatzss;

eta_mat=[0,0,0;
         0,0,0;
         0,0,0;
         0,0,0;
         0,0,0;
         1,0,0;
         0,SDV_THETA,0;
         0,0,SDV_ZA];


cd(home3);