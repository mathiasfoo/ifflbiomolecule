function dx= Protein_Hill_Model(t,x,p)
global par 
dx = zeros(9,1);


%%%%molecule species

X = x(1);  
Y = x(2);       
TetR = x(3);
Pz_rep = x(4);
Z = x(5);     
GFP = x(6);
XZ = x (7);
XY = x (8);
aTc = x(9);

Pz = par.P_z - Pz_rep;

dx(1) = (par.IPTG^p(16)/(p(17)^p(16)+par.IPTG^p(16)))*par.P_x*p(3) - p(12)*X*Z - p(12)*X*Y - p(6)*X;

dx(2) = par.P_y*p(2) - p(12)*X*Y - p(5)*Y;

dx(3) = p(7)*XY- p(8)*TetR - p(11)*TetR*Pz - p(18)*aTc*TetR;

dx(4) = p(11)*TetR*Pz;

dx(5) = p(1)*Pz*(par.IPTG^p(16)/(p(17)^p(16)+ par.IPTG^p(16))) - p(12)*X*Z - p(4)*Z;

dx(6) = p(9)*XZ - p(10)*GFP;

dx (7) = p(12)*X*Z - p(13)*XZ; 

dx (8) = p(12)*Y*X - p(14)*XY;

dx (9) = - p(18)*aTc*TetR; 
