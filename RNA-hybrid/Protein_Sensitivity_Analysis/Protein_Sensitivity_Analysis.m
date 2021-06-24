clear all
% % close all
clc
global par p



par.P_z = 1e-9;
par.P_y = 1e-9;
par.P_x = 1e-9; 

alpha_z = 1;  
alpha_y = 1;
alpha_x = 0.05;

delta_z = 0.005;
delta_y = 0.05;
delta_x = 0.005;

alpha_tetR = 0.5;
delta_tetR = 0.0005;

alpha_g = 0.05;  
delta_g = 0.0005;   
  
omega = 5e3;        
gamma = 1e4;         

delta_xz = 0.001;    
delta_xy = 0.001;    

ze = 10;
nu = 0.5;
K_1 = 1e-3;
Beta = 1e4;         

 p = [alpha_z alpha_y alpha_x delta_z delta_y delta_x alpha_tetR delta_tetR ...
     alpha_g delta_g omega gamma delta_xz delta_xy ze nu K_1 Beta];

 
%%%%%%%%%%%% Nominal Value Calculation %%%%%%%%%%%%%%%%%%%%%%%% 
 
nu_atc = [200 100 20 2].';    %%% aTc Data Points
nu_IPTG = [1 0.4 0.02 0.01];  %%% IPTG Data Points
atc_conv = 0.46822;           %%% Conversion Factor for aTc (ng/ml to mM) 
nol = 0; 

    tspan=0:0.1:450*60;  %%% seconds
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    
figure
hold on
for j = 1:4
   
    par.IPTG = nu_IPTG(j)*10^-3; 
    
    for k = 1:4 
        
        nol = nol + 1; 
        
        par.aTc = (nu_atc(k)/atc_conv)*10^-9;
        x0 = [0 0 0 0 0 0 0 0 par.aTc];
        [t,x]=ode23s(@(t,x)Protein_Hill_Model(t,x,p),tspan,x0, options);
        TetR = x(:,6)*10^(p(15));
        Simu_t = t./60;
        
        if j <= 3
                col = ['k','r','b','g'];
                subplot(4,4,nol)
                plot(Simu_t, TetR, col(1), 'LineWidth',2)
                set(gca,'FontName','Times New Roman')
                set(gca,'FontSize',18)
                ylim([0 6500])
                hold on
        else 
                subplot(4,4,nol)
                plot(Simu_t, TetR, col(1), 'LineWidth',2)
                set(gca,'FontName','Times New Roman')
                set(gca,'FontSize',18)
                xlabel('Time (min)')
                ylim([0 6500])  
                hold on 
         end
        
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%% Sensativity Analysis %%%%%%%%%%%%%%%%%%%%%%%%



%%% Inital Parameter Values
 p0 = [alpha_z alpha_y alpha_x delta_z delta_y delta_x alpha_tetR delta_tetR ...
     alpha_g delta_g omega gamma delta_xz delta_xy ze nu K_1 Beta];    

lb = p0.'*0.5; %%%lower bound
ub = p0.'*1.5; %%%upper bound

dV = zeros(18,1);
Pval = zeros(18,10);

for i = 1:18
    dV(i) = (ub(i) - lb(i))/10;
    
    for j = 1:10
        Pval(i,j) = lb(i) + dV(i)*1/2 + dV(i)*(j-1);
    end
    
end


nu_atc = [200 100 20 2].'; 
nu_IPTG = [1 0.4 0.02 0.01]; 

for ij = 1:10  %%% iterating over number of parameter sets
    
    alpha_a = Pval(1,randi(10));

    alpha_b = Pval(2,randi(10));
    alpha_c = Pval(3,randi(10));
    delta_a = Pval(4,randi(10));
    delta_b = Pval(5,randi(10));
    delta_c = Pval(6,randi(10));
    alpha_tetR = Pval(7,randi(10));

    delta_tetR = Pval(8,randi(10));
    alpha_g = Pval(9,randi(10));
    delta_g = Pval(10,randi(10));
    omega = Pval(11,randi(10));
    gamma = Pval(12,randi(10));
    delta_ac = Pval(13,randi(10));
    delta_bc = Pval(14,randi(10));

    ze = 10;  %%%% Scaling Factor (Constant)
    nu = Pval(16,randi(10));
    K_1 = Pval(17,randi(10));
    Beta = Pval(18,randi(10));
  

    p = [alpha_a alpha_b alpha_c delta_a delta_b delta_c alpha_tetR delta_tetR ...
    alpha_g delta_g omega gamma delta_ac delta_bc ze nu K_1 Beta];

    
    subp = 0;
    
    for nj = 1:4
        par.IPTG = nu_IPTG(nj)*10^-3; 
            
        for mj = 1:4
            subp = subp+1;
            par.aTc = (nu_atc(mj)/atc_conv)*10^-9;
            x0 = [0 0 0 0 0 0 0 0 par.aTc];
            [t,x]=ode23s(@(t,x)Protein_Hill_Model(t,x,p),tspan,x0, options);
            TetR = x(:,6)*10^(p(15));
            Simu_t = t./60;

            subplot(4,4,subp)
            plot(Simu_t, TetR, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5)
            set(gca,'FontName','Times New Roman')
            set(gca,'FontSize',18)
            ylim([0 6500])
            hold on
        end
    end
end

%%% Title Label
subplot(4, 4, 1)
title ('aTc = 400 nM ')
subplot(4, 4, 2)
title ('aTc = 200 nM ')
subplot(4, 4, 3)
title ('aTc = 40 nM ')
subplot(4, 4, 4)
title ('aTc = 4 nM ')
h=text(500,75,'I = 1.0 mM','fontweight','bold');
set(h,'Rotation',90)
set(h,'FontName','Times')
set(h,'FontSize',18)

subplot(4,4,8)
h=text(500,75,'I = 0.4 mM','fontweight','bold');
set(h,'Rotation',90)
set(h,'FontName','Times')
set(h,'FontSize',18)

subplot(4,4,12)
h=text(500,75,'I = 0.02 mM','fontweight','bold');
set(h,'Rotation',90)
set(h,'FontName','Times')
set(h,'FontSize',18)

subplot(4,4,16)
h=text(500,75,'I = 0.01 mM','fontweight','bold');
set(h,'Rotation',90)
set(h,'FontName','Times')
set(h,'FontSize',18)


%%% GFP Label
subplot(4, 4, 1)
plot(0, 0)
ylabel('GFP (a.u.)','fontweight','bold')
subplot(4, 4, 5)
plot(0, 0)
ylabel('GFP (a.u.)','fontweight','bold')
subplot(4, 4, 9)
plot(0, 0)
ylabel('GFP (a.u.)','fontweight','bold')
subplot(4, 4, 13)
plot(0, 0)
ylabel('GFP (a.u.)','fontweight','bold')
hold off



 
 
 

