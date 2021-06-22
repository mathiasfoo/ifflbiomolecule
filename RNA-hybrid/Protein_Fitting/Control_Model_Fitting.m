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

delta_xy = 0.001; 
delta_xz = 0.001; 

ze = 10;
nu = 0.5;
K_1 = 1e-3; 
Beta = 1e4;

 p0 = [alpha_z alpha_y alpha_x delta_z delta_y delta_x alpha_tetR delta_tetR ...
     alpha_g delta_g omega gamma delta_xz delta_xy ze nu K_1 Beta];


A = [];
b = [];
Aeq = [];
beq = [];
nlcon = [];
lb = [0.01 0.01 0.01 0.0001 0.0001 0.0001 0.0005 0.00005 0.0005 0.00005 1e2 1e2 0.0001 0.0001 1 0.1 1e-6 1e2]; %%%lower bound
ub = [10 10 10 0.1 0.1 0.1 0.5 0.05 0.5 0.05 1e6 1e6 0.1 0.1 100 4 1e-3 1e6]; %%%upper bound

Evalall = [];

p = fmincon(@Obj,p0,A,b,Aeq,beq,lb,ub,nlcon);


% show final objective
disp(['Final SSE Objective: ' num2str(Obj(p))])

col = ['k','r','b','g'];






