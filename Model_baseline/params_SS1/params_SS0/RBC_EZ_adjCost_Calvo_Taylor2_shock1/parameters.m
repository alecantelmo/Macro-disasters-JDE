home5=pwd;

cd ..\;
cd RBC_EZ_adjCost_Calvo_Taylor2

parameters;

LAMBDA_MU=0; % need to use zero trend, because i assumed that in the previous models. To make nonzero trend, change the previous models.
SDV_ZMU=exp(-6.0283); % from Fernandez-Villaverde et al (2015)

eta_mat=[0,0,0,0;
         0,0,0,0;
         0,0,0,0;
         0,0,0,0;
         0,0,0,0;
         0,0,0,0;
         1,0,0,0;
         0,SDV_THETA,0,0;
         0,0,SDV_ZA,0;
         0,0,0,SDV_ZMU];


cd(home5);