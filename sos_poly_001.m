clear all,close all,clc;
%% sostools400
dir1=pwd
cd('C:\Users\Artun\Documents\MATLAB\SOSTOOLS400')
run('addsostools.m')
cd(dir1)
clear all,close all,clc;
%% this is for "SDPT3"
current_dir=pwd;
cd('C:\Users\Artun\Documents\MATLAB\sdpt3');
run('install_sdpt3.m')
cd(current_dir);
clear all,close all,clc;
%% mosek 2015a
addpath(genpath('C:\Program Files\Mosek\9.3\toolbox\R2015a'))
%% yalmip
addpath(genpath('C:\Users\Artun\Documents\MATLAB\YALMIP-master'))
% "clear classes"
run('yalmiptest.m')
%%






%% install MOSEK or SDPT3
%% install YALMIP

%% YALMIP sos study
%% problem 1: is this poly sos ?
clear all,close all,clc;yalmip('clear');
x = sdpvar(1,1);
y = sdpvar(1,1);
p = (1+x)^4 + (1-y)^2;
% generate the monomials
Z = monolist([x y],degree(p)/2); % sdisplay(Z)
Q = sdpvar(length(Z));
p_sos = Z'*Q*Z;
% optimize
F = [coefficients(p-p_sos,[x y]) == 0, Q >= 0];
optimize(F)
value(Q)
eig(value(Q))
















%% YALMIP sos study
%% problem 2: minimize a poly [UNCONSTRAINED]
clear all,close all,clc;yalmip('clear');
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
p = (4-2.1*x1^2+(x1^4)/3)*x1^2+x1*x2+(-4+4*x2^2)*x2^2;
%%
Z = monolist([x1 x2],degree(p)/2); % sdisplay(Z)
Q = sdpvar(length(Z));
p_sos = Z'*Q*Z;
%%
sdpvar t
F = [coefficients((p-t)-p_sos,[x1 x2]) == 0, Q >= 0];
optimize(F,-t)
sdisplay(p_sos)
value(t)








%% YALMIP sos study
%% problem 3: minimize a poly [CONSTRAINED] polynomial optimization
% sostools-400-manual-pdf-page-40
clear all,close all,clc;yalmip('clear');
x1 = sdpvar(1,1);       x2 = sdpvar(1,1);
f = x1+x2;          % term to be minimized
g1=x1;              % in-eq-constraint-1
g2=x2-0.5;          % in-eq-constraint-2
h1=x1^2+x2^2-1;     % eq-constraint-2
h2=x2-x1^2-0.5;     % eq-constraint-2
deg1=6;
% sigma_0
m_sigma_0 = monolist([x1 x2],deg1);       % sdisplay(m_sigma_0)
Q_sigma_0 = sdpvar(length(m_sigma_0));
sigma_0 = m_sigma_0'*Q_sigma_0*m_sigma_0; % sdisplay(sigma_0)
% sigma_1
m_sigma_1 = monolist([x1 x2],deg1);
Q_sigma_1 = sdpvar(length(m_sigma_1)); 
sigma_1 = m_sigma_1'*Q_sigma_1*m_sigma_1;
% sigma_2
m_sigma_2 = monolist([x1 x2],deg1);
Q_sigma_2 = sdpvar(length(m_sigma_2));
sigma_2 = m_sigma_2'*Q_sigma_2*m_sigma_2;
% sigma_1_2
m_sigma_1_2 = monolist([x1 x2],deg1);
Q_sigma_1_2 = sdpvar(length(m_sigma_1_2));
sigma_1_2 = m_sigma_1_2'*Q_sigma_1_2*m_sigma_1_2;
% lambda_1
[lambda_1,coeff_1,mono_list_1] = polynomial([x1;x2],deg1); % sdisplay(lambda_1)
% lambda_2
[lambda_2,coeff_2,mono_list_2] = polynomial([x1;x2],deg1);
%
sdpvar t
F = [coefficients((f-t)-1*[sigma_0+lambda_1*h1+lambda_2*h2+sigma_1*g1+sigma_2*g2+sigma_1_2*g1*g2],[x1 x2]) == 0,...
    Q_sigma_0 >= 0, Q_sigma_1 >= 0, Q_sigma_2 >= 0, Q_sigma_1_2 >= 0];
options = sdpsettings('solver','bisection','bisection.solver','mosek'); % 'sdpt3'
% diagnostics = bisection(Constraint,Objective,options)
diagnostics = bisection([F],[-t],options)
value(t)

















%