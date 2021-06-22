function dx = RNA_ONLY_ODE(t,x)

global iratio aratio p

dx = zeros(4,1);

%%%%molecule species
X = x(1); %%% X mRNA 
Y = x(2); %%% Y mRNA
Z = x(3); %%% Z mRNA
GFP = x(4); %%% GPF expression, simplified

%% For Parameter Estimation

parAra = (6.6)*1e-3; %Arabinose 0.1% -> 6.6 mM 
parIPTG = (1e-3); %%% Molar
parP_x = 1e-9;
parP_y = 1e-9; 
parP_z = 1e-9; 

Py_star = parP_y*((iratio*parIPTG)^p(6)/((p(7))^p(6)+(iratio*parIPTG)^p(6))); 
Pz_star = parP_z*((iratio*parIPTG)^p(6)/((p(7))^p(6)+(iratio*parIPTG)^p(6))); 


%%% X mRNA 
dx(1) = real(p(1)*parP_x*((aratio*parAra)^p(2)/(abs(p(3))^p(2)+(aratio*parAra)^p(2))) - p(4)*X*Pz_star - p(4)*X*Py_star - p(5)*X);

%%% Y mRNA
dx(2) = real(p(8)*p(4)*X*Py_star - p(9)*Y - (p(10))*Y*Z);

%%% Z mRNA 
dx(3) = real(p(11)*p(4)*X*Pz_star - p(12)*Z - (p(10))*Y*Z);

%%% Simplified degradable GFP
dx(4) = real(p(13)*Z - p(14)*GFP);
