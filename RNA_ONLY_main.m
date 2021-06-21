clear all
% close all
clc

global iratio aratio p


load RNA_Exp_Data_Ara_0_1_IPTG_1_2.mat
load RNA_Exp_Data_Ara_0_1_IPTG_1_8.mat
load RNA_Exp_Data_Ara_0_1_IPTG_1_32.mat
load RNA_Exp_Data_Ara_0_1_IPTG_1_128.mat

load RNA_Exp_Data_Ara_0_1_4_IPTG_1_2.mat
load RNA_Exp_Data_Ara_0_1_4_IPTG_1_8.mat
load RNA_Exp_Data_Ara_0_1_4_IPTG_1_32.mat
load RNA_Exp_Data_Ara_0_1_4_IPTG_1_128.mat

load RNA_Exp_Data_Ara_0_1_16_IPTG_1_2.mat
load RNA_Exp_Data_Ara_0_1_16_IPTG_1_8.mat
load RNA_Exp_Data_Ara_0_1_16_IPTG_1_32.mat
load RNA_Exp_Data_Ara_0_1_16_IPTG_1_128.mat

load RNA_Exp_Data_Ara_0_1_64_IPTG_1_2.mat
load RNA_Exp_Data_Ara_0_1_64_IPTG_1_8.mat
load RNA_Exp_Data_Ara_0_1_64_IPTG_1_32.mat
load RNA_Exp_Data_Ara_0_1_64_IPTG_1_128.mat

ExpData = [outGFP_B_RJ_A2_Ara_0_1_IPTG_1_2';...
    outGFP_B_RJ_A2_Ara_0_1_IPTG_1_8';...
    outGFP_B_RJ_A2_Ara_0_1_IPTG_1_32';...
    outGFP_B_RJ_A2_Ara_0_1_IPTG_1_128';...
    outGFP_B_RJ_A2_Ara_0_1_4_IPTG_1_2';...
    outGFP_B_RJ_A2_Ara_0_1_4_IPTG_1_8';...
    outGFP_B_RJ_A2_Ara_0_1_4_IPTG_1_32';...
    outGFP_B_RJ_A2_Ara_0_1_4_IPTG_1_128';...
    outGFP_B_RJ_A2_Ara_0_1_16_IPTG_1_2';...
    outGFP_B_RJ_A2_Ara_0_1_16_IPTG_1_8';...
    outGFP_B_RJ_A2_Ara_0_1_16_IPTG_1_32';...
    outGFP_B_RJ_A2_Ara_0_1_16_IPTG_1_128';...
    outGFP_B_RJ_A2_Ara_0_1_64_IPTG_1_2';...
    outGFP_B_RJ_A2_Ara_0_1_64_IPTG_1_8';...
    outGFP_B_RJ_A2_Ara_0_1_64_IPTG_1_32';...
    outGFP_B_RJ_A2_Ara_0_1_64_IPTG_1_128'];


prCmRNA = [1.563,0.045332,0.00016796,11186.9458,0.00054103];
prMM = [0.18503,0.0015127];
prBmRNA = [0.92384,0.07247,5602.5417];
prAmRNA = [0.0032535,0.0063173];
prGFP = [0.020225,0.00019489,444619888155.746];

p0 = [prCmRNA,prMM,prBmRNA,prAmRNA,prGFP];


%%% Specifiying IPTG concentration, 1 == 500uM
iratiototal = [1/2];
iratiostring = {'500\muM'};

%%% Specifiying Arabinose concentration, 1 == 6600uM
aratiototal = [1];
aratiostring = {'6600\muM'};

for aa = 1:length(aratiototal)
    aratio = aratiototal(aa);
    for ii = 1:length(iratiototal)
        iratio = iratiototal(ii);
        p = [prCmRNA,prMM,prBmRNA,prAmRNA,prGFP];
        
        All_x = cell(1,1);
        All_t = cell(1,1);
        
        tspan = 0:(10*60):((4)*60*60); % Hour X Minutes X Seconds
        options = odeset('RelTol',1e-10,'AbsTol',1e-10);
        
        x0 = [0 0 0 0];
        [t,x] = ode23s('RNA_ONLY_ODE',tspan,x0);
        L = length(tspan);
        All_x{1,1} = x;
        All_t{1,1} = t./60;
        
        %%% Figure Plotting
        figure(1)
        plot(All_t{1,1},p(15).*All_x{1,1}(:,4),'-','LineWidth',2,'Color',[0.68,0.68,0.68])
        hold on
        
        
        %%% Experimental Data Lookup Value, ExpData(w,7:end), where w equals
        % (1) A=1, I=1/2, (2) A=1, I=1/8, (3) A=1, I=1/32, (4) A=1,I=1/128
        % (5) A=1/4, I=1/2, (6) A=1/4, I=1/8, (7) A=1/4, I=1/32, (8) A=1/4,I=1/128
        % (9) A=1/16, I=1/2, (10) A=1/16, I=1/8, (11) A=1/16, I=1/32, (12) A=1/16,I=1/128
        % (13) A=1/64, I=1/2, (14) A=1/64, I=1/8, (15) A=1/64, I=1/32, (16) A=1/64,I=1/128
        
        plot(All_t{1,1},ExpData(1,7:end),'x-','LineWidth',2)
        hold on
        title(['I = ', iratiostring{ii}, ', A = ', aratiostring{aa}],'FontSize', 9)

        xlim([0 240])
        ylim([0 200])
    end
end

