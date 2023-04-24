% This file presents the results

% clear,
clc, close all
load('models')
main_folder=pwd;

all_versions=[1]; % do disaster parameterizaion

names={'tp3'};

% only cheap solutions

if exist('options.mat','file')~=0
    load('options');
else
    do_all_solutions=1;
end

if do_all_solutions==0
    names={'pert1','pert2','pert3','tp1','tp2','tp3','smol1','smol2'};
end

for veri=all_versions

    version=['ver' num2str(veri)];


    n_models=length(models);
    n_methods=length(names);
    smol3_grid_expansion=zeros(n_models,1);

    maxEE=zeros(n_models,n_methods);
    meanEE=zeros(n_models,n_methods);
    runtime=zeros(n_models,n_methods);
    tily_mean=zeros(n_models,n_methods);
    tilc_mean=zeros(n_models,n_methods);
    tilx_mean=zeros(n_models,n_methods);
    consumption_fall=zeros(n_models,n_methods);
    
    tilw_mean=zeros(n_models,n_methods);
    
    tauc_mean=zeros(n_models,n_methods);
    B_mean=zeros(n_models,n_methods);
    R_mean=zeros(n_models,n_methods);
    

    for modeli=1:length(models)

        load(['performance\' version '\' models{modeli} '\accuracy_and_runtime'],'pert*','tp*','smol*','shocks','T')

        % Euler Errors - max and mean
        t_start=100; t_end=1000; % take a sample from t_start to t_end

        for j=1:n_methods
           eval(['solution=' names{j} '_r;'])
           maxEE(modeli,j)=max(log10(max(abs(solution(:,t_start:t_end)))));
           meanEE(modeli,j)=mean(log10(max(abs(solution(:,t_start:t_end))))); % mean across simulation points of max error across equations.
        end

        % runtime
        for i=1:length(names)
            eval(['runtime(modeli,i)=' names{i} '_time;']);
        end

        % Simulation moments

        load([main_folder '\performance\' version '\' models{modeli} '\simulation'],'results_*','x_tp*','x_pert*','x_smol*','nu_*')

        for i=1:length(names)
            load([main_folder '\files_for_TaylorProjection\' models{modeli} '\order3\model'],'subsvars','y','x')
            eval(['results=results_' names{i} '(:,t_start:t_end);'])
            eval(['x_simul=x_' names{i} '(:,t_start:t_end);'])
            eval(['nu=nu_' names{i} '(:,t_start:t_end);'])

            loghatz=nu(logical(subsvars==sym('loghatz')),:);
            logtheta=x_simul(logical(x==sym('logtheta')),:);
            logtilxback=x_simul(logical(x==sym('logtilxback')),:);
            logtilxgback=x_simul(logical(x==sym('logtilxgback')),:);
            grant=exp(x_simul(logical(x==sym('loggrant')),:));
            logtilw=nu(logical(subsvars==sym('logtilw')),:);  
            tilw=exp(logtilw);
            tilv_tilvss=exp(nu(logical(subsvars==sym('logtilv_tilvss')),:));
            tilv_tilvss_2=exp(nu(logical(subsvars==sym('logtilv_tilvss_2')),:));
            tilv_2=exp(nu(logical(subsvars==sym('logtilv_2')),:));
            tilB=x_simul(logical(x==sym('logBback')),:);
            tilB=exp(tilB);
            logtauc=x_simul(logical(x==sym('logtaucback')),:);
            tauc=exp(logtauc);
            R=exp(nu(logical(subsvars==sym('logR')),:));
            logpi=nu(logical(subsvars==sym('logpi')),:);
            pilevel=exp(logpi);
            RSTAR=exp(nu(logical(subsvars==sym('logRSTAR')),:));
            tily=nu(logical(subsvars==sym('tily')),:);
            tilc=nu(logical(subsvars==sym('tilc')),:);
            tilx=nu(logical(subsvars==sym('tilx')),:);
            logtilk=nu(logical(subsvars==sym('logtilk')),:);
            tilxg=nu(logical(subsvars==sym('tilxg')),:);
            tilk=exp(logtilk);
            logtilkg=nu(logical(subsvars==sym('logtilkg')),:);
            tilkg=exp(logtilkg);
            tilnx=nu(logical(subsvars==sym('tilnx')),:);
            logby=nu(logical(subsvars==sym('logby')),:);
            by=exp(logby);

             %variables with trend 
             logtily=nu(logical(subsvars==sym('logtily')),:);
             logz=cumsum(loghatz);
             logy=logtily+logz;
             ytrend=exp(logy);
             
             logtilx=log(tilx);
             xtrend=exp(logtilx+logz)';
             logtilc=log(tilc);
             ctrend=exp(logtilc+logz)';
             granttrend=exp(log(grant)+logz);
             
             hatz=exp(loghatz);
                      
             xgtrend=tilxg.*hatz;
             btrend=tilB.*hatz; 
             
             

             %
            tily_mean(modeli,i)=mean(tily);
            tilc_mean(modeli,i)=mean(tilc);
            tilx_mean(modeli,i)=mean(tilx);
            tilk_mean(modeli,i)=mean(tilk);
            tilxg_mean(modeli,i)=mean(tilxg);
            tauc_mean(modeli,i)=mean(tauc);
            B_mean(modeli,i)=mean(tilB);
            R_mean(modeli,i)=mean(R);
            pi_mean(modeli,i)=mean(pilevel);
            by_mean(modeli,i)=mean(by);
            RSTAR_mean(modeli,i)=mean(RSTAR);
            tilnx_mean(modeli,i)=mean(tilnx);
            tilv_tilvss_mean(modeli,i)=mean(tilv_tilvss);
            tilv_tilvss_2_mean(modeli,i)=mean(tilv_tilvss_2);
            tilv_2_mean(modeli,i)=mean(tilv_2);
            tily_std(modeli,i)=std(tily);
            tilc_std(modeli,i)=std(tilc);
            tilx_std(modeli,i)=std(tilx);
            tilk_std(modeli,i)=std(tilk);
            ytrend_mean(modeli,i)=mean(ytrend);
            ctrend_mean(modeli,i)=mean(ctrend);
            xtrend_mean(modeli,i)=mean(xtrend);
            xgtrend_mean(modeli,i)=mean(xgtrend);
            
            % Construct GDP Growth
            for jj=1:length(ytrend)-4
                ygrowth(1,jj)=100*(ytrend(1,jj+4)/ytrend(1,jj)-1);
            end
            ygrowth_mean(modeli,i)=mean(ygrowth);
            
            ytrend=exp(logy)';
        end
        save(['present_' version '_' models{modeli}])
    end
    % present tables
    if veri==0
        disp('no-disaster parametrization')
        disp(names)
        meanEE
        disp(names)
        maxEE
    elseif veri==1
        disp('disaster parametrization')
        disp(names)
        meanEE
        disp(names)
        maxEE
        disp(names)
        tily_mean
        disp(names)
        tilc_mean
        disp(names)
        tilx_mean
        disp(names)
        tilxg_mean
        disp(names)
        B_mean
        disp(names)
        by_mean
        disp(names)
        tauc_mean
        disp(names)
        RSTAR_mean
        disp(names)
        ygrowth_mean
        disp(names) 
        tily_std
        disp(names) 
        tilc_std
        disp(names) 
        tilx_std
        disp(names)
        tilnx_mean
%         disp(names)
        ytrend_mean
         disp(names)
        ctrend_mean
         disp(names)
        xtrend_mean
         disp(names)
        xgtrend_mean
        runtime
    end
    
    welfare=tilv_2_mean;
vector_means=[tily_mean tilc_mean tilx_mean by_mean welfare ygrowth_mean]'
divergence_vars=[ytrend ctrend xtrend];


disaster_dummy=x_tp3(11,100:end);
damages=exp(x_tp3(12,100:end));

for ii=1:size(tily,2)
    
    pubK_gdp(1,ii)=100*tilkg(1,ii)/tily(1,ii);
    pubINV_gdp(1,ii)=100*tilxg(1,ii)/tily(1,ii);
end

for ii=1:size(tily,2)-1
grant_dis(1,ii)=(disaster_dummy(1,ii)*grant(1,ii+1));
end

grant_dis=grant_dis';

end



