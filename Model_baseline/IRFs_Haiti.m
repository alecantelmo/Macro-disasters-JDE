
by_tp3_high=exp(logby_tp3);
logy_dev_tp3_high=logy_dev_tp3;
cy_dev_tp3_high=cy_dev_tp3;
xy_dev_tp3_high=xy_dev_tp3;
ky_dev_tp3_high=ky_dev_tp3;
xgy_dev_tp3_high=xgy_dev_tp3;
kgy_dev_tp3_high=kgy_dev_tp3;
logtaucback_dev_tp3_high=logtaucback_dev_tp3;
TFPgrowth_tp3_high=TFPgrowth_tp3.*100;


load baseline_replic %%% load baseline results

 by_tp3=exp(logby_tp3);
 
 for iii=1:size(logtily_tp3,2)
     by_dev_tp3(1,iii)=(by_tp3(1,iii)-by_tp3(1,1))*100;
     by_dev_tp3_high(1,iii)=(by_tp3_high(1,iii)-by_tp3_high(1,1))*100;
 end 

       time=0:length(tvec);
       
subplot(3,3,1);
hh1=plot(time(1,1:21),logy_dev_tp3(1,1:21),'LineWidth',2);
hold on
hh2=plot(time(1,1:21),logy_dev_tp3_high(1,1:21),':r','LineWidth',3);
title('Output','Fontsize',20)
xlim([0 20])
axis tight
       
subplot(3,3,2)
 plot(time(1,1:21),cy_dev_tp3(1,1:21),'LineWidth',2);
hold on
plot(time(1,1:21),cy_dev_tp3_high(1,1:21),':r','LineWidth',3);
title('Consumption to GDP','Fontsize',20)
xlim([0 20])
axis tight

subplot(3,3,3)
plot(time(1,1:21),xy_dev_tp3(1,1:21),'LineWidth',2);
hold on
plot(time(1,1:21),xy_dev_tp3_high(1,1:21),':r','LineWidth',3);
title('Private Investment to GDP','Fontsize',20)
xlim([0 20])
axis tight

subplot(3,3,4)
plot(time(1,1:21),ky_dev_tp3(1,1:21),'LineWidth',2);
hold on
plot(time(1,1:21),ky_dev_tp3_high(1,1:21),':r','LineWidth',3);
title('Private Capital to GDP','Fontsize',20)
xlim([0 20])
axis tight
        
subplot(3,3,5)
plot(time(1,1:21),xgy_dev_tp3(1,1:21),'LineWidth',2);
hold on
plot(time(1,1:21),xgy_dev_tp3_high(1,1:21),':r','LineWidth',3);
title('Public Investment to GDP','Fontsize',20)
xlim([0 20])
axis tight

subplot(3,3,6)
plot(time(1,1:21),kgy_dev_tp3(1,1:21),'LineWidth',2);
hold on
plot(time(1,1:21),kgy_dev_tp3_high(1,1:21),':r','LineWidth',3);
title('Public Capital to GDP','Fontsize',20)
xlim([0 20])
axis tight
        
subplot(3,3,7)
plot(time(1,1:21),logtaucback_dev_tp3(1,1:21),'LineWidth',2);
hold on
plot(time(1,1:21),logtaucback_dev_tp3_high(1,1:21),':r','LineWidth',3);
title('Tax rate','Fontsize',20)
xlim([0 20])
axis tight
 
subplot(3,3,8)
plot(time(1,1:21),by_dev_tp3(1,2:22),'LineWidth',2);
hold on
plot(time(1,1:21),by_dev_tp3_high(1,2:22),':r','LineWidth',3);
title('Public Debt to Annual GDP','Fontsize',20)
xlim([0 20])
ylim ([-5 15])
axis tight
        
subplot(3,3,9)
plot(time(1,1:21),TFPgrowth_tp3(1,1:21),'LineWidth',2);  
hold on
plot(time(1,1:21),TFPgrowth_tp3_high(1,1:21),':r','LineWidth',3);
title('TFP Growth Rate','Fontsize',20)
xlim([0 20])
axis tight
        
h_legend = legend([hh1,hh2],{'Average disaster','Hurricane Matthew'}, 'Location', 'Best', 'Orientation', 'vertical');
set(h_legend, 'FontSize', 16)

