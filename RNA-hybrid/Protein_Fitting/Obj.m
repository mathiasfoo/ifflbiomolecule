function objective = Obj(p) 

global par
load Exp_Data_V4.mat;

With_TetR_1_200 = TetR_1_200./1000;
With_TetR_1_20 = TetR_1_20./1000;
With_TetR_1_2 = TetR_1_2./1000; 

With_TetR_2_200 = TetR_2_200./1000;
With_TetR_2_20 = TetR_2_20./1000;
With_TetR_2_2 = TetR_2_2./1000; 


tspan=0:0.1:430*60;  %%% seconds
options = odeset('RelTol',1e-10,'AbsTol',1e-10); 


    

atc_conv = 0.46822; 

ratio = [200 20 200 20];

for j = 1:4 
    
    Simout_TetR_j = zeros(length(tspan),1);
        
    par.aTc = (ratio(j)/atc_conv)*10^-9;
    x0 = [0 0 0 0 0 0 0 0 par.aTc];
  
    if j <= 2
        
        par. IPTG = (1*10^-3); 
        par.P_y = 1e-9;
        [t,x]=ode23s(@(t,x)Protein_Hill_Model(t,x,p),tspan,x0, options);
        Simout_TetR_j = x(1:6000:end,6)*10^(p(15));
        
    else 
        par. IPTG = 0.1*10^-3; 
        par.P_y = 1e-9;
        [t,x]=ode23s(@(t,x)Protein_Hill_Model(t,x,p),tspan,x0, options);
        Simout_TetR_j = x(1:6000:end,6)*10^(p(15));
    end

 
    if j == 1 
    Simout_1 = Simout_TetR_j;
    end
    if j == 2
    Simout_2 = Simout_TetR_j;
    end
    if j == 3 
    Simout_3 = Simout_TetR_j;
    end
    if j == 4
    Simout_4 = Simout_TetR_j;
    end


    
    
    
end


objective = sum((Simout_1 - With_TetR_1_200).^2)/max(With_TetR_1_200)+ sum((Simout_2 - With_TetR_1_20).^2)/max(With_TetR_1_20)...
  + sum((Simout_3 - With_TetR_2_200).^2)/max(With_TetR_2_200) + sum((Simout_4 - With_TetR_2_20).^2)/max(With_TetR_2_20);





disp(['SSE Objective: ' num2str(objective)])
disp(['Simu1_End: ' num2str(Simout_1(end)),  'Exp1: ' num2str(With_TetR_1_200(end))])
disp(['Simu2_End: ' num2str(Simout_2(end)),  'Exp2: ' num2str(With_TetR_1_20(end))])
disp(['Simu4_End: ' num2str(Simout_3(end)),  'Exp4: ' num2str(With_TetR_2_200(end))])
disp(['Simu5_End: ' num2str(Simout_4(end)),  'Exp5: ' num2str(With_TetR_2_20(end))])




end
