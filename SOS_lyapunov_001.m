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
%% problem : find a Lyap [UNCONSTRAINED] with SOS "WORKING" 3D-problem
clear all,close all,clc;yalmip('clear');
tic
x1 = sdpvar(1,1);   x2 = sdpvar(1,1); x3 = sdpvar(1,1);  % state-variables "indeterminates"
f1 = -x1-10*x2^2;   f2 = -2*x2;     f3 = -3*x3;    % the dynamics
deg1=4;             eps1=1e-2;
eps1_poly=eps1*(x1^2+x2^2+x3^2);            % small-sos-poly
% V
[V,coeff_1,mono_list_1] = polynomial([x1;x2;x3],deg1); % sdisplay(mono_list_1)
% Vdot
mVdot=-jacobian(V,[x1;x2;x3])*[f1;f2;f3];
options = sdpsettings('solver','mosek'); % 'sdpt3' 'mosek'
% constraints "V>0","-Vdot>0","V(0)=0","Vdot(0)=0"
F=[];
F=[F;replace(V,[x1,x2,x3],[0,0,0])==0];
F=[F;replace(mVdot,[x1,x2,x3],[0,0,0])==0];
F=[F;replace(V,[x1,x2,x3],[1,0,0])==1];
F=[F; 0<=coeff_1<=100];
F=[F;sos(V-eps1_poly)];
F=[F;sos(mVdot-eps1_poly)];
% [sol,u,Q] = solvesos(Constraints,Objective,options,decisionvariables)
[sol,u,Q] = solvesos(F,[],options,[coeff_1])
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
if contains(V_char,'internal(3)')
    V_char=replace(V_char,'internal(3)','x3');
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
syms x1 x2 x3 real
evalin('base',['V=',V_char,';'])
V=simplify(V)
vpa(V,3)
%% is V pos-def ?
V_char=char(vpa(V,3))
clearvars -except V_char
yalmip('clear');
x1 = sdpvar(1,1);   x2 = sdpvar(1,1); x3 = sdpvar(1,1);
evalin('base',['V=',V_char,';']); sdisplay(V)

Z = monolist([x1 x2 x3],degree(V)/2);
Q = sdpvar(length(Z));
V_sos = Z'*Q*Z;

F = [coefficients(V-V_sos,[x1 x2 x3]) == 0, Q >= 0];
optimize(F)
%% is mVdot pos-def ?
V_char
clearvars -except V_char
yalmip('clear');
x1 = sdpvar(1,1);   x2 = sdpvar(1,1); x3 = sdpvar(1,1);
evalin('base',['V=',V_char,';']); sdisplay(V)

f1 = -x1-10*x2^2;   f2 = -2*x2; f3 = -3*x3;    % the dynamics
mVdot=-jacobian(V,[x1;x2;x3])*[f1;f2;f3];
Z = monolist([x1 x2 x3],degree(mVdot)/2);
Q = sdpvar(length(Z));
V_sos = Z'*Q*Z;

F = [coefficients(mVdot-V_sos,[x1 x2 x3]) == 0, Q >= 0];
optimize(F)
%% finally, output the V and mVdot expressions
V_char
clearvars -except V_char
syms x1 x2 x3 real
evalin('base',['V=',V_char,';'])
f1 = -x1-10*x2^2;   f2 = -2*x2; f3 = -3*x3;    % the dynamics
mVdot=-jacobian(V,[x1;x2;x3])*[f1;f2;f3];

V=expand(V);
mVdot=expand(mVdot);
vpa(V,3)
vpa(mVdot,3)



%%



%% let us look at some random trajectories
clear all,close all,clc;
fig1=figure(1);fig1.Color=[1,1,1];
ax1=axes('Parent',fig1);
    set(0,'CurrentFigure',fig1);
    set(fig1,'currentaxes',ax1);
for ii=1:1:10
    tspan=[0:0.01:6]; x0=[randi([-5,5],3,1)];
    wt=tspan;
    f=randi([1,10],1,1); w=sin(2*pi*f*tspan);% w=square(tspan);
    [t,x]=ode45(@(t,x) odefcn(t,x,wt,w),tspan,x0);
    plot(t,x(:,1),'r-','LineWidth',[1],"Parent",ax1); hold on;
    plot(t,x(:,2),'g-','LineWidth',[1],"Parent",ax1); hold on;
    plot(t,x(:,3),'b-','LineWidth',[1],"Parent",ax1); hold on;
end
yline(1);hold on;yline(-1);hold on;
axis square

function xdot=odefcn(t,x,wt,w)
w=interp1(wt,w,t);
xdot=zeros(3,1);
x1=x(1);x2=x(2);x3=x(3);
f1 = -x1-10*x2^2;   f2 = -2*x2; f3 = -3*x3;    % the dynamics

xdot(1)=f1;
xdot(2)=f2;
xdot(3)=f3;
end