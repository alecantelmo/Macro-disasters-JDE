% This file produces impulse response functions of a disaster shock.

clear,clc
load('models') 
main_folder=pwd;
addpath(main_folder)

all_versions={'ver1'}; 
models=models(1); 

for veri=1:length(all_versions)
    version=all_versions{veri};
    load(version)

    for modeli=1:length(models)
        % Take the model equations for accuracy tests from this path:
        addpath([main_folder '\files_for_TaylorProjection\' models{modeli} '\order3\files']);

        load([main_folder '\all_models\' models{modeli} '\TaylorProjection\choose']);

        % load some parameters
        load([main_folder '\' pert_folder '\' models{modeli} '\order3\perturbation'],'eta_mat','nyss','nxss','params','P1','nep1','n_e','prob_disaster','MUD')
        nep=nep1;P=P1;
        load([main_folder '\files_for_TaylorProjection\' models{modeli} '\order3\model'],'model','y','x','subsvars','symparams');
        model0=model;

        n_y=model0.n_y;
        n_x1=model0.n_x1;
        n_x=model0.n_x;
        n_f=model0.n_f;
        for k=2:5
            [U,W]=create_UW(n_x,k);
            model.U{k}=U;
            model.W{k}=W;
        end

        c0=nxss;
        eta=eta_mat;
        
        % load solutions
        load([main_folder '\performance\' version '\' models{modeli} '\accuracy_and_runtime'],'pert*','tp*','smol*','shocks','T')

        %Impulse response functions
        impulse=1;

        if impulse==1
            ir_shocks=[zeros(size(shocks,1),120)];
            ir_shocks(1,:)=-MUD;
            ir_shocks(1,100)=1-MUD;

            % Taylor projection
            poly.type='power';
           
            coeffs=tp3(:); order=3;
            [ results,nu ] = fun_simulation( nxss,ir_shocks,coeffs,c0,order,n_y,n_x,n_f,model,params,eta );
            ir_tp3=results;
            nu_tp3=nu;
            tp3_EE=eval_eqs(coeffs,results(n_y+1:end,100),model,params,eta,c0,nep,P,order,poly);


            save('ir','ir_*','nu_*')
        else
            load('ir')
        end

        tvec=99:130; T=length(tvec);
          names={'tp3'};
        for i=1:length(names)
            disp(names{i})
            eval(['results=ir_' names{i} ';'])
            eval(['nu=nu_' names{i} ';'])
            
    
            
            
            tildiv=nu(logical(subsvars==sym('tildiv')),:);
            loghatz=nu(logical(subsvars==sym('loghatz')),:);
            AnnTFPgrowth=loghatz.*400;
            logz=cumsum(loghatz);
            logtily=nu(logical(subsvars==sym('logtily')),:);
            logtilc=nu(logical(subsvars==sym('logtilc')),:);
            logtilx=results(logical(y==sym('logtilx')),:);
            logy=logtily+logz;
            logx=logtilx+logz;
            logc=logtilc+logz;
            dlogy=logtily(2:end)-logtily(1:end-1)+loghatz(2:end);
            Growth=exp(4*dlogy)-1;
            logtilB=nu(logical(x==sym('logBback')),:);
            logB=logtilB.*logz;
            logtaucback=results(logical(x==sym('logtaucback')),:);
            logR=nu(logical(subsvars==sym('logR')),:);
            logpi=nu(logical(subsvars==sym('logpi')),:);
            logtilxg=results(logical(y==sym('logtilxg')),:);
            logxg=logtilxg+logz;
            logtilkg=nu(logical(subsvars==sym('logtilkg')),:);
            logtilk=nu(logical(subsvars==sym('logtilk')),:);
            logk=logtilk+logz;
            logkg=logtilkg+logz;
            logby=nu(logical(subsvars==sym('logby')),:);
            logRreal=logR-logpi;
            logRSTAR=nu(logical(subsvars==sym('logRSTAR')),:);
            logRSTARreal=logRSTAR;%-logpi;
            xy_ratio=logtilx-logtily;
            xgy_ratio=logtilxg-logtily; 
            by=exp(logby);
            ky_ratio=logtilk-logtily;
            kgy_ratio=logtilkg-logtily;
            cy_ratio=logtilc-logtily;

            eval(['Growth_' names{i} '=Growth(tvec(1):end);'])
            eval(['logtily_' names{i} '=logtily(tvec(1):end);'])
            eval(['logy_' names{i} '=logy(tvec(1):end);'])
            eval(['logtilx_' names{i} '=logtilx(tvec(1):end);'])
            eval(['logx_' names{i} '=logx(tvec(1):end);'])
            eval(['logtilc_' names{i} '=logtilc(tvec(1):end);'])
            eval(['logc_' names{i} '=logc(tvec(1):end);'])
            eval(['logB_' names{i} '=logB(tvec(1):end);'])
            eval(['logtaucback_' names{i} '=logtaucback(tvec(2):tvec(21)+1);'])
            eval(['logtilB_' names{i} '=logtilB(tvec(1):end);'])
            eval(['logR_' names{i} '=logR(tvec(1):end);'])
            eval(['logpi_' names{i} '=logpi(tvec(1):end);'])
            eval(['by_' names{i} '=by(tvec(1):end);'])
            eval(['logtilxg_' names{i} '=logtilxg(tvec(1):end);'])
            eval(['logxg_' names{i} '=logxg(tvec(1):end);'])
            eval(['logtilkg_' names{i} '=logtilkg(tvec(1):end);'])
            eval(['logtilk_' names{i} '=logtilk(tvec(1):end);'])
            eval(['xy_ratio_' names{i} '=xy_ratio(tvec(1):end);'])
            eval(['xgy_ratio_' names{i} '=xgy_ratio(tvec(1):end);'])
            eval(['logby_' names{i} '=logby(tvec(1):end);'])
            eval(['logRreal_' names{i} '=logRreal(tvec(1):end);'])
            eval(['logRSTARreal_' names{i} '=logRSTARreal(tvec(1):end);'])
            eval(['ky_ratio_' names{i} '=ky_ratio(tvec(1):end);'])
            eval(['kgy_ratio_' names{i} '=kgy_ratio(tvec(1):end);'])
            eval(['AnnTFPgrowth_' names{i} '=AnnTFPgrowth(tvec(1):end);'])
            eval(['logz_' names{i} '=logz(tvec(1):end);'])
            eval(['cy_ratio_' names{i} '=cy_ratio(tvec(1):end);'])
            eval(['TFPgrowth_' names{i} '=loghatz(tvec(1):end);'])
%           
            
        end
        close all
        time=0:length(tvec);
        
        % construct variables in deviations from steady state
        logtily_dev_tp3=zeros(size(logtily_tp3));
        logtilc_dev_tp3=zeros(size(logtilc_tp3));
        logtilx_dev_tp3=zeros(size(logtilx_tp3));
        logy_dev_tp3=zeros(size(logy_tp3));
        logc_dev_tp3=zeros(size(logc_tp3));
        logx_dev_tp3=zeros(size(logx_tp3));
        logtaucback_dev_tp3=zeros(size(logtaucback_tp3));
        logtilB_dev_tp3=zeros(size(logtilB_tp3));
        logB_dev_tp3=zeros(size(logB_tp3));
        logR_dev_tp3=zeros(size(logR_tp3));
        logpi_dev_tp3=zeros(size(logpi_tp3));
        by_dev_tp3=zeros(size(by_tp3));
        logtilxg_dev_tp3=zeros(size(logtilxg_tp3));
        logxg_dev_tp3=zeros(size(logxg_tp3));
        xy_dev_tp3=zeros(size(xy_ratio_tp3));
        xgy_dev_tp3=zeros(size(xgy_ratio_tp3));
        logby_dev_tp3=zeros(size(logby_tp3));
        logRreal_dev_tp3=zeros(size(logRreal_tp3));
        logRSTARreal_dev_tp3=zeros(size(logRSTARreal_tp3));
        ky_dev_tp3=zeros(size(ky_ratio_tp3));
        kgy_dev_tp3=zeros(size(kgy_ratio_tp3));
        logtilk_dev_tp3=zeros(size(logtilk_tp3));
        logtilkg_dev_tp3=zeros(size(logtilkg_tp3));
        cy_dev_tp3=zeros(size(cy_ratio_tp3));

        
        for iii=1:size(logtily_tp3,2)
        logtily_dev_tp3(1,iii)=(logtily_tp3(1,iii)-logtily_tp3(1,1))*100;
        logtilc_dev_tp3(1,iii)=(logtilc_tp3(1,iii)-logtilc_tp3(1,1))*100;
        logtilx_dev_tp3(1,iii)=(logtilx_tp3(1,iii)-logtilx_tp3(1,1))*100;
        logy_dev_tp3(1,iii)=(logy_tp3(1,iii)-logy_tp3(1,1))*100;
        logc_dev_tp3(1,iii)=(logc_tp3(1,iii)-logc_tp3(1,1))*100;
        logx_dev_tp3(1,iii)=(logx_tp3(1,iii)-logx_tp3(1,1))*100;
        logR_dev_tp3(1,iii)=(logR_tp3(1,iii)-logR_tp3(1,1))*100;
        logpi_dev_tp3(1,iii)=(logpi_tp3(1,iii)-logpi_tp3(1,1))*100;
        by_dev_tp3(1,iii)=(by_tp3(1,iii)-by_tp3(1,1))*400;
        logtilxg_dev_tp3(1,iii)=(logtilxg_tp3(1,iii)-logtilxg_tp3(1,1))*100;
        logxg_dev_tp3(1,iii)=(logxg_tp3(1,iii)-logxg_tp3(1,1))*100;
        xy_dev_tp3(1,iii)=(xy_ratio_tp3(1,iii)-xy_ratio_tp3(1,1))*100;
        xgy_dev_tp3(1,iii)=(xgy_ratio_tp3(1,iii)-xgy_ratio_tp3(1,1))*100;
        logRreal_dev_tp3(1,iii)=(logRreal_tp3(1,iii)-logRreal_tp3(1,1))*100;
        logRSTARreal_dev_tp3(1,iii)=(logRSTARreal_tp3(1,iii)-logRSTARreal_tp3(1,1))*100;
        ky_dev_tp3(1,iii)=(ky_ratio_tp3(1,iii)-ky_ratio_tp3(1,1))*100;
        kgy_dev_tp3(1,iii)=(kgy_ratio_tp3(1,iii)-kgy_ratio_tp3(1,1))*100;
        logtilk_dev_tp3(1,iii)=(logtilk_tp3(1,iii)-logtilk_tp3(1,1))*100;
        logtilkg_dev_tp3(1,iii)=(logtilkg_tp3(1,iii)-logtilkg_tp3(1,1))*100;
        cy_dev_tp3(1,iii)=(cy_ratio_tp3(1,iii)-cy_ratio_tp3(1,1))*100;
        end
        
        for iii=1:size(logtaucback_tp3,2)
        logtaucback_dev_tp3(1,iii)=(logtaucback_tp3(1,iii)-logtaucback_tp3(1,1))*100;
        end
        
        for iii=1:size(logtily_tp3,2)
        logtilB_dev_tp3(1,iii)=(logtilB_tp3(1,iii)-logtilB_tp3(1,2))*100;
        logB_dev_tp3(1,iii)=(logB_tp3(1,iii)-logB_tp3(1,2))*100;
        
        logby_dev_tp3(1,iii)=400*(logby_tp3(1,iii)-logby_tp3(1,2));
        end


        load([main_folder '\performance\' version '\' models{modeli} '\accuracy_and_runtime'],'pert*','tp*','smol*','shocks','T')

    end
end

