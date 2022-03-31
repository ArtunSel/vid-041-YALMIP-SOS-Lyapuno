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



%% add MOSEK AND YALMIP
%% problem : find a Lyap [UNCONSTRAINED] with integer-coeffs 2D-example
clear all,close all,clc;yalmip('clear');
x1 = sdpvar(1,1);   x2 = sdpvar(1,1);  % state-variables "indeterminates"
f1 = -x1-10*x2^2;   f2 = -2*x2;        % the dynamics
f1 = x2;   f2 = -2*x1-3*x2;        % the dynamics --------------------------------------------------
deg1=4;             eps1=1e-2;
eps1_poly=eps1*(x1^2+x2^2);            % small-sos-poly
% V
[V,coeff_1,mono_list_1] = polynomial([x1;x2],deg1); sdisplay(mono_list_1)
% Vdot
mVdot=-jacobian(V,[x1;x2])*[f1;f2];
options = sdpsettings('solver','mosek'); % 'sdpt3' 'mosek'
% constraints "V>0","-Vdot>0","V(0)=0","Vdot(0)=0"
F=[];
F=[F;replace(V,[x1,x2],[0,0])==0];
F=[F;replace(mVdot,[x1,x2],[0,0])==0];
x_val=[];
x_val=[x_val;randn(1e2,2)*1e-3];  
x_val=[x_val;randn(1e2,2)*1e-2];
x_val=[x_val;randn(1e2,2)*1e-1];
x_val=[x_val;randn(1e2,2)*1];
x_val=[x_val;randn(1e2,2)*1e1];
x_val=[x_val;randn(1e2,2)*1e2]; % plot(x_val(:,1),x_val(:,2),'r*')
for ii=1:1:size(x_val,1)
    x1_val=x_val(ii,1);x2_val=x_val(ii,2);
    F=[F;replace(V,[x1,x2],[x1_val,x2_val])>=eps1*(x1_val.^2+x2_val.^2)];
    F=[F;replace(mVdot,[x1,x2],[x1_val,x2_val])>=eps1*(x1_val.^2+x2_val.^2)]; % '0'
end
F=[F;replace(V,[x1,x2],[1,0])>=1];
F=[F;replace(V,[x1,x2],[1,0])<=10];
F=[F; 0<=coeff_1<=20];
optimize([F;integer(coeff_1)],[],options) 
%% get "V"
coeff_1_value=value(coeff_1)
V_char=cell2mat(sdisplay(V));
if contains(V_char,'internal(1)')
    V_char=replace(V_char,'internal(1)','x1');
end
if contains(V_char,'internal(2)')
    V_char=replace(V_char,'internal(2)','x2');
end
if contains(V_char,'f1')
    V_char=replace(V_char,'f1','x2');
end
for ii=1:1:length(coeff_1)
    if abs(coeff_1_value(ii))<1e-3
        coeff_1_value(ii)=0;
    end
    V_char=replace(V_char,['coeff_1(',num2str(ii),')'],[num2str(coeff_1_value(ii))]);
end
V_char
%% get "V" as sym
clearvars -except V_char
syms x1 x2 real
evalin('base',['V=',V_char,';'])
V=simplify(V)
vpa(V,3)
%% is V pos-def ?
V_char=char(vpa(V,3))
clearvars -except V_char
yalmip('clear');
x1 = sdpvar(1,1);   x2 = sdpvar(1,1);
evalin('base',['V=',V_char,';']); sdisplay(V)

Z = monolist([x1 x2],degree(V)/2);
Q = sdpvar(length(Z));
V_sos = Z'*Q*Z;

F = [coefficients(V-V_sos,[x1 x2]) == 0, Q >= 0];
optimize(F)
%% is mVdot pos-def ?
V_char
clearvars -except V_char
yalmip('clear');
x1 = sdpvar(1,1);   x2 = sdpvar(1,1);
evalin('base',['V=',V_char,';']); sdisplay(V)

f1 = x2;   f2 = -2*x1-3*x2;        % the dynamics --------------------------------------------------
mVdot=-jacobian(V,[x1;x2])*[f1;f2];
Z = monolist([x1 x2],degree(mVdot)/2);
Q = sdpvar(length(Z));
V_sos = Z'*Q*Z;

F = [coefficients(V-V_sos,[x1 x2]) == 0, Q >= 0];
optimize(F)
%% finally, output the V and mVdot expressions
V_char
clearvars -except V_char
syms x1 x2 real
evalin('base',['V=',V_char,';'])
f1 = x2;   f2 = -2*x1-3*x2;        % the dynamics --------------------------------------------------
mVdot=-jacobian(V,[x1;x2])*[f1;f2];

V=expand(V);
mVdot=expand(mVdot);
vpa(V,3)
vpa(mVdot,3)










%