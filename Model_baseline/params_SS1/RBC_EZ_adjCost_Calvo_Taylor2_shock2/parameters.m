home6=pwd;

cd ..\;
cd RBC_EZ_adjCost_Calvo_Taylor2_shock1

parameters;

SDV_M=exp(-6); % from Fernandez-Villaverde et al (2015)

eta_mat=[0,0,0,0,0;
         0,0,0,0,0;
         0,0,0,0,0;
         0,0,0,0,0;
         0,0,0,0,0;
         0,0,0,0,0;
         1,0,0,0,0;
         0,SDV_THETA,0,0,0;
         0,0,SDV_ZA,0,0;
         0,0,0,SDV_ZMU,0;
         0,0,0,0,SDV_M];


cd(home6);