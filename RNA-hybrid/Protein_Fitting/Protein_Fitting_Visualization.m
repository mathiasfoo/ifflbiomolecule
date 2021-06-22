clear all
% % close all
clc
global par p
%%% Load Results %%%
load Final_Fitting_Results.mat
p = p;

par.P_z = 1e-9;
par.P_y = 1e-9;
par.P_x = 1e-9;


atc_conv = 0.46822; 

nu_atc = [200 100 20 2].';
nu_IPTG = [1 0.1 0.02]; 
atc_conv = 0.46822; 
nol = 0; 

    
    tspan=0:1:450*60;  %%% seconds
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    
figure
hold on
for j = 1:3
   
    par. IPTG = nu_IPTG(j)*10^-3; 
    
    for k = 1:4 
        
        nol = nol + 1; 
        
        par. aTc = (nu_atc(k)/atc_conv)*10^-9;
        x0 = [0 0 0 0 0 0 0 0 par.aTc];
        [t,x]=ode23s(@(t,x)Protein_Hill_Model(t,x,p),tspan,x0, options);
        TetR = x(:,6)*10^(p(15))*1000;
        Simu_t = t./60;
        
        if j <= 2
                col = [[0 0.5 1],'r','b','g'];
                subplot(3,4,nol)
                plot(Simu_t, TetR, col(1), 'LineWidth',2)
                set(gca,'FontName','Times New Roman')
                set(gca,'FontSize',18)
                ylim([0 9000])
                hold on
        else 
                subplot(3,4,nol)
                plot(Simu_t, TetR, col(1), 'LineWidth',2)
                set(gca,'FontName','Times New Roman')
                set(gca,'FontSize',18)
                xlabel('Time (min)')
                ylim([0 9000])  
                hold on 
         end
        
    end
    
end

wid = 0.75;

load Exp_Data_V4.mat
load Exp_Valid_V1.mat
%% IPTG 1mM %%%
subplot(3, 4, 1)
plot(Time, TetR_1_200,'+','LineWidth',wid)
title('aTc = 400 nM ')
ylabel('GFP (a.u.)','fontweight','bold')


subplot(3, 4, 2)
plot(Time, TetR_1_100, '+', 'LineWidth',wid)
title('aTc = 200 nM ')

subplot(3, 4, 3)
plot(Time, TetR_1_20,'+', 'LineWidth',wid)
title('aTc = 40 nM ')


subplot(3, 4, 4)
plot(Time, TetR_1_2,'+', 'LineWidth',wid)
title('aTc = 4 nM ')



%%% IPTG 0.1 mM %%%
subplot(3, 4, 5)
plot(Time, TetR_2_200,'+', 'LineWidth',wid)
ylabel('GFP (a.u.)','fontweight','bold')


subplot(3, 4, 6)
plot(Time, TetR_2_100,'+', 'LineWidth',wid)

subplot(3, 4, 7)
plot(Time, TetR_2_20,'+', 'LineWidth',wid)


subplot(3, 4, 8)
plot(Time, TetR_2_2,'+', 'LineWidth',wid)


%%% IPTG 0.02 mM %%%
subplot(3, 4, 9)
plot(Time, TetR_02_200,'+', 'LineWidth',wid)
ylabel('GFP (a.u.)', 'fontweight','bold')
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',18)

subplot(3, 4, 10)
plot(Time, TetR_02_100,'+', 'LineWidth',wid)

subplot(3, 4, 11)
plot(Time, TetR_02_20,'+', 'LineWidth',wid)

subplot(3, 4, 12)
plot(Time, TetR_02_2,'+', 'LineWidth',wid)
legend('Model', 'Experimental Results')


%%% IPTG Label 
subplot(3,4,4)
h=text(500,2000,'I = 1.0 mM','fontweight','bold');
set(h,'Rotation',90)
set(h,'FontName','Times')
set(h,'FontSize',18)

subplot(3,4,8)
h=text(500,2000,'I = 0.1 mM','fontweight','bold');
set(h,'Rotation',90)
set(h,'FontName','Times')
set(h,'FontSize',18)

subplot(3,4,12)
h=text(500,2000,'I = 0.02 mM','fontweight','bold');
set(h,'Rotation',90)
set(h,'FontName','Times')
set(h,'FontSize',18)


hold off



