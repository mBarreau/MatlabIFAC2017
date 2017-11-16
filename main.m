%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Simulation on a time-delay system with an average on the output     %%
%                                           Version 1.0 / November 2017   %
%                                                                         %
% By M. Barreau                                                           %
%                                                                         %
%   If you are using or modifying this code, please cite the following    %
%   reference:                                                            %
%   M. Barreau, A. Seuret, F. Gouaisbaut,                                 %
%   Wirtinger-based Exponential Stability for Time-Delay Systems,         %
%   IFAC World Congress, Toulouse, Volume 50, Issue 1, 2017               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reset
clear
close all
clc

%% System definition

% Definition of the system
s = tf('s');
A = [4,11,-30;1,0,0;0,1,0];
B = [1 1; 1 0; 0 1];
C = eye(size(A,1));
ts = 0.1; % Sampling time for the noise

%% Systems

% LTI system
h = 0.2;
K = place(A, B,[-2, -2.1, -2.2]*1.4);
Fr = -pinv(C/(A-B*K)*B);

% Observer
ho = 0.2;
Ko = [5.65637568067253,9.81891791571489,-23.7635778046581;
    0.930479228115671,2.51720684281594,1.78630954592076;
    0.00980709730736698,1.11438020419894,2.37954156181310];
L = place(A, B,[-2, -2.1, -2.2]*1.5);
Fro = -pinv(C/(A-B*L)*B);

% Controller
hc = 0.2; % alpha = 0.5
Kc = [-2.61211513880737,-4.46536632682474,10.4643618039906;
    0.211800076405891,-0.519360003926725,-3.36528795968091];
AD = B*Kc/hc;
Ad = zeros(size(A,1));
Frc = -pinv(C/(A+B*Kc)*B);

%% Simulation

% Parameters
init = 10;
ref = 0;
time = 5;
initialState = init*ones(size(A,1),1)*0.95; % For the observer

% Noise delay
tn = max([h hc ho]);

% Simulation without noise
P = 0; % Power of noise set to 0
sim('simu');
t = xNorm.time-max([h hc ho]);
select = t>(time*0.75);

f = figure;
plot(t, yNorm.signals.values, t, yoNorm.signals.values, '--', ...
    t, ycNorm.signals.values, '-.');
xlim([min(t) max(t)]);
legend('Non delayed system', 'Observer + Controller', 'Controller',...
    'Location', 'Best');
title('Convergence in norm of x without noise');
xlabel('Time [s]');
ylabel('||x||_2');

mean([xNorm.signals.values(select),...
    xoNorm.signals.values(select),...
    xcNorm.signals.values(select)])

% Simulation with noise
P = 0.01; % Power of noise non null
sim('simu');
t = xNorm.time-max([h hc ho]);
select = t>(time*0.75);

yNorm = yNorm.signals.values;
yoNorm = yoNorm.signals.values;
ycNorm = ycNorm.signals.values;

figure
plot(t, yNorm, t, yoNorm, '--', t, ycNorm, '-.');
xlim([min(t) max(t)]);
legend('Non delayed system', 'Observer + Controller', 'Controller', ...
    'Location', 'Best');
title('Convergence in norm of y with noise');
xlabel('Time [s]');
ylabel('||y||_2');

mean([yNorm(select) yoNorm(select) ycNorm(select)])
