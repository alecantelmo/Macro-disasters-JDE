home7=pwd;

cd ..\;
cd RBC_EZ_adjCost_Calvo_Taylor2_shock2

parameters;

RHOXI=0.1182; % from Fernandez-Villaverde et al (2015)
SDV_XI=exp(-1.9834); % from Fernandez-Villaverde et al (2015)


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


cd(home7);