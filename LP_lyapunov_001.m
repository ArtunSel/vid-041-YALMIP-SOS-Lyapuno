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

%% BEFORE RUNNING THE CODE ADD "YALMIP" AND "SDPT3" libraries
%     set(findall(gcf,'type','line'),'linewidth',[1]);
% yalmip('clear');
%%

%% add MOSEK
%% add YALMIP

%% YALMIP sos study
%% problem : find a Lyap [UNCONSTRAINED] with integer-coeffs
clear all,close all,clc;yalmip('clear');
tic
x1 = sdpvar(1,1);   x2 = sdpvar(1,1);  % state-variables "indeterminates"
f1 = -x1-10*x2^2;   f2 = -2*x2;        % the dynamics
deg1=4;             eps1=1e-3;
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
x_val=[x_val;[rand(1e2,2)-0.5]*1e-3];  
x_val=[x_val;[rand(1e2,2)-0.5]*1e-2];
x_val=[x_val;[rand(1e2,2)-0.5]*1e-1];
x_val=[x_val;[rand(1e2,2)-0.5]*1];
x_val=[x_val;[rand(1e2,2)-0.5]*1e1];
x_val=[x_val;[rand(1e2,2)-0.5]*1e2];
x_val=[x_val;[rand(1e2,2)-0.5]*1e3];% plot(x_val(:,1),x_val(:,2),'r*')
for ii=1:1:size(x_val,1)
    x1_val=x_val(ii,1);x2_val=x_val(ii,2);
    F=[F;replace(V,[x1,x2],[x1_val,x2_val])>=eps1*(x1_val.^2+x2_val.^2)];
    F=[F;replace(mVdot,[x1,x2],[x1_val,x2_val])>=0]; % '0'
end
F=[F;replace(V,[x1,x2],[1,0])==1];
% F=[F;replace(V,[x1,x2],[0,1])>=1];F=[F;replace(V,[x1,x2],[0,1])<=100];
% F=[F;replace(V,[x1,x2],[1,1])>=2];F=[F;replace(V,[x1,x2],[1,1])<=100];
% F=[F;replace(mVdot,[x1,x2],[1,0])>=1];F=[F;replace(mVdot,[x1,x2],[1,0])<=100];
% F=[F;replace(mVdot,[x1,x2],[0,1])>=1];F=[F;replace(mVdot,[x1,x2],[0,1])<=100];
% F=[F;replace(mVdot,[x1,x2],[1,1])>=2];F=[F;replace(mVdot,[x1,x2],[1,1])<=100];
F=[F; 0<=coeff_1<=20];
optimize([F],[sum(coeff_1)],options) % '[norm(coeff_1,1)]'
toc
%% get "V"
coeff_1_value=value(coeff_1)
V_char=cell2mat(sdisplay(V));
if contains(V_char,'internal(1)')
    V_char=replace(V_char,'internal(1)','x1');
end
if contains(V_char,'internal(2)')
    V_char=replace(V_char,'internal(2)','x2');
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

f1 = -x1-10*x2^2;   f2 = -2*x2;
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
f1 = -x1-10*x2^2;   f2 = -2*x2;
mVdot=-jacobian(V,[x1;x2])*[f1;f2];

V=expand(V);
mVdot=expand(mVdot);
vpa(V,3)
vpa(mVdot,3)

%% if we are able to find a lyap-fcn, we are done at this stage!
%% now we have the "monomial-terms", and maybe we can find a Lyap-fcn with int-coeffs
% '5.22*x1*x2 + 1.55*x1*x2^2 + 0.0246*x1^2*x2 + 0.00209*x1*x2^3 + 0.999*x1^2 + 20.0*x2^2 + 11.5*x2^3 + 20.0*x2^4'
clear all,close all,clc;yalmip('clear');
x1 = sdpvar(1,1);   x2 = sdpvar(1,1);  % state-variables "indeterminates"
f1 = -x1-10*x2^2;   f2 = -2*x2;        % the dynamics

while(true)
coeff_1=randi([0,10],5,1);
% V=[x1*x2,x1*x2^2,x1^2*x2,x1*x2^3,x1^2,x2^2,x2^3,x2^4]*coeff_1;
V=[x1*x2,x1*x2^3,x1^2,x2^2,x2^4]*coeff_1;
if coeff_1(1)==0 && coeff_1(2)==0 && coeff_1(3)==0
    continue;
end
if coeff_1(1)==0 && coeff_1(2)==0 && coeff_1(4)==0 && coeff_1(5)==0
    continue;
end
mVdot=-jacobian(V,[x1;x2])*[f1;f2];
% "V>0" check
Z1 = monolist([x1 x2],degree(V)/2);
Q1 = sdpvar(length(Z1));
V_sos = Z1'*Q1*Z1;
F1 = [coefficients(V-V_sos,[x1 x2]) == 0, Q1 >= 0];
diagnostics=optimize(F1) % diagnostics = optimize(Constraints,Objective,options)
if diagnostics.problem==0
    coeff_1

    Z2 = monolist([x1 x2],degree(mVdot)/2);
    Q2 = sdpvar(length(Z2));
    mVdot_sos = Z2'*Q2*Z2;
    F2 = [coefficients(mVdot-mVdot_sos,[x1 x2]) == 0, Q2 >= 0];
    diagnostics=optimize(F2) 
    if diagnostics.problem==0
        disp('solved!');
        coeff_1
        break;
    else
        continue;
    end
else
    continue;
end
end




%