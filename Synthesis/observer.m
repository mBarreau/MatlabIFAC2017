function [pres,K,P,R,S] = observer( A, C, h, e1, e2, e3, e4, alpha )
%OBSERVER Generates a static feedback gain for an average system
%   This function generates a static feedback gain for a system with an
%   average on the state as the output to make an observer.
%   The error e between the state and the observed one is as follow:
%
%   de/dt = Ae - 1/h*K*C*int_{t-h}^t e(s)ds
%
%   A and C are matrices of the observed system:
%
%   dx/dt = Ax +Bu
%   y = 1/h*C*int_{t-h}^t x(s)ds
%
%   h >= 0 is the window for average
%   e1, e2, e3 and e4 are non null real for the algorithm
%   alpha is a non negative real related to the exponential stability
%
%   PRES is the primal residual of the optimization process. All the values
%   need to be negative for the solution to be correct.
%   K is the feedback gain
%   P, R and S are the matrices of the Lyapunov functionnal. 
%
%   [PRES, K, P, R, S] = OBSERVER(A,C,h) computes the feedback gain for
%   e1=e2=e3=e4=1 and alpha=0.
%   [PRES, K, P, R, S] = OBSERVER(A,C,h,alpha) computes the feedback gain for
%   e1=e2=e3=e4=1.
%   [PRES, K, P, R, S] = OBSERVER(A,C,h,e1,e2,e3,e4) computes the feedback gain for
%   alpha=0.
%   [PRES, K, P, R, S] = OBSERVER(A,C,h,e1,e2,e3,e4,alpha) computes the feedback gain.
%
%   Version 1.0 / November 2017
% 
%   If you are using or modifying this code, please cite the following
%   reference:
%   M. Barreau, A. Seuret, F. Gouaisbaut,
%   Wirtinger-based Exponential Stability for Time-Delay Systems,
%   IFAC World Congress, Toulouse, Volume 50, Issue 1, 2017, Pages 11984-11989
%
%   See also CONTROLLER

%% Discussion on the argument
switch nargin
    case 3
        e1 = 1;
        e2 = 1;
        e3 = 1;
        e4 = 1;
        alpha = 0;
    case 4
        alpha = e1;
        e1 = 1;
        e2 = 1;
        e3 = 1;
        e4 = 1;
    case 5
        e1 = 1;
        e2 = 1;
        e3 = 1;
        e4 = 1;
    case 7
        alpha = 0;
    case 8
    otherwise
        warning('OBSERVER does not have a correct number of inputs');
end

%% Variables assignment
n = size(A,1); % size of the system

% SDP variables
S = sdpvar(n);
P = sdpvar(2*n);
R = sdpvar(n);
Z = sdpvar(n, n, 'full');
Kbar = sdpvar(n,size(C, 1), 'full');
beta1 = sdpvar(1);

% General matrices
Sbar = blkdiag(S, -S*exp(-2*alpha*h), zeros(2*n));
Rtilde = blkdiag(R, 3*R);
W = [e1*eye(n) e2*eye(n) e3*eye(n) e4*eye(n)];

% Constraints matrices
F0 = [zeros(n) zeros(n) eye(n) zeros(n);
    eye(n) -eye(n) zeros(n) zeros(n)];
F1 = [eye(n) zeros(n) zeros(n) zeros(n);
    zeros(n) zeros(n) zeros(n) h*eye(n)];
F2 = [eye(n) -eye(n) zeros(n) zeros(n);
    eye(n) eye(n) zeros(n) -2*eye(n)];
F3 = [zeros(n) zeros(n) eye(n) zeros(n)];
N = [A zeros(n) -eye(n) zeros(n)];
N2 = [zeros(n, 3*n) -Kbar*C];

% M Matrix
He = F1'*P*F0;
M = He + He' + 2*alpha*F1'*P*F1 + Sbar + h*F3'*R*F3;
if h > 0
    M = M - (exp(-2*alpha*h)/h)*F2'*Rtilde*F2;
end

% T Matrix
T = N'*Z*W + N2'*W;
T = T + T';

% Optimization matrices
Q1 = P+2*exp(-2*h*alpha)*[h*R -R; -R R/h]+...
    blkdiag(zeros(n), S*exp(-2*alpha*h)/h)...
    -beta1*blkdiag(eye(n),zeros(n));

% Solve Lyap
problem = [M+T < -1e-5, S > 1e-5, R > 1e-5, P > 1e-5, Q1>0, beta1 > 1e-5];
objective = {};

%% Optimization
option = sdpsettings('verbose', 0, 'solver', 'sdpt3');
optimize(problem, objective, option);
[pres,~] = checkset(problem);

%% Values
K = value(Z')\value(Kbar);
P = value(P);
R = value(R);
S = value(S);

end

