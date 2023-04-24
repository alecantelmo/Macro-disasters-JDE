home2=pwd;

cd ..\;
cd RBC_EZ_adjCost

parameters;

THETA_P=0.8139;  %from Fernandez-Vaillaverde et al (2015)
CHI=0.6186; %from Fernandez-Vaillaverde et al (2015)
EPSILON=10; %from Fernandez-Vaillaverde et al (2015)

PISS=1.02^0.25; % Inflation target of 2 percent annually
GAMMA_PI=1.3; % Inflation coefficient in Taylor rule.

SteadyState;
RSS=PISS/exp(logmnomss-logmdenss);

eta_mat=[0,0,0;
         0,0,0;
         0,0,0;
         0,0,0;
         1,0,0;
         0,SDV_THETA,0;
         0,0,SDV_ZA];


cd(home2);