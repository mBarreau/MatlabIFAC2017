function [pres,K,P,R,S] = controller(A, B, h, e1, e2, e3, e4, alpha)
%CONTROLLER Generates a static feedback gain for an average system
%   This function generates a static feedback gain for a system with an
%   average on the state as the output.
%   The controlled system is as follow:
%
%   dx/dt = Ax - 1/h*B*K*int_{t-h}^t x(s)ds
%
%   A and B are matrices as stated above.
%   h >= 0 is the window for average
%   e1, e2, e3 and e4 are non null real for the algorithm
%   alpha is a non negative real related to the exponential stability
%
%   PRES is the primal residual of the optimization process. All the values
%   need to be negative for the solution to be correct.
%   K is the feedback gain
%   P, R and S are the matrices of the Lyapunov functionnal. 
%
%   [PRES, K, P, R, S] = CONTROLLER(A,B,h) computes the feedback gain for
%   e1=e2=e3=e4=1 and alpha=0.
%   [PRES, K, P, R, S] = CONTROLLER(A,B,h,alpha) computes the feedback gain for
%   e1=e2=e3=e4=1.
%   [PRES, K, P, R, S] = CONTROLLER(A,B,h,e1,e2,e3,e4) computes the feedback gain for
%   alpha=0.
%   [PRES, K, P, R, S] = CONTROLLER(A,B,h,e1,e2,e3,e4,alpha) computes the feedback gain.
%
%   Version 1.0 / November 2017
% 
%   If you are using or modifying this code, please cite the following
%   reference:
%   M. Barreau, A. Seuret, F. Gouaisbaut,
%   Wirtinger-based Exponential Stability for Time-Delay Systems,
%   IFAC World Congress, Toulouse, Volume 50, Issue 1, 2017, Pages 11984-11989
%
%   See also OBSERVER

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
        warning('CONTROLLER does not have a correct number of inputs');
end

%% Variables assignment
n = size(A,1); % size of the system

% SDP variables
S = sdpvar(n);
R = sdpvar(n);
P = sdpvar(2*n);
X = sdpvar(n, n, 'full');
Kbar = sdpvar(size(B, 2), n, 'full');
beta1 = sdpvar(1);

% General variables
Sbar = blkdiag(S, -S*exp(-2*alpha*h), zeros(2*n));
Rtilde = blkdiag(R, 3*R);
Xtilde = blkdiag(X,X,X,X);
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
N2 = [zeros(n, 3*n) B*Kbar];

% M Matrix
He = F1'*P*F0;
M = He + He' + 2*alpha*F1'*P*F1 + Sbar + h*F3'*R*F3;
if h > 0
    M = M - (exp(-2*alpha*h)/h)*F2'*Rtilde*F2;
end

% T matrix
T = ((N*Xtilde+N2)'*W);
T = T + T';

% Optimization matrices
Q1 = P+2*exp(-2*h*alpha)*[h*R -R; -R R/h]+...
    blkdiag(zeros(n), S*exp(-2*alpha*h)/h)...
    -beta1*blkdiag(eye(n),zeros(n));

% Solve Lyap
problem = [M+T < -1e-5, S > 1e-5, R > 1e-5, Q1 > 1e-5, beta1 > 1e-5];
objective = {};

%% Optimization
option = sdpsettings('verbose', 0, 'solver', 'sdpt3');
optimize(problem, objective, option);
[pres,~] = checkset(problem);

%% Values
X = value(X);
Kbar = value(Kbar);
K = Kbar/X;
P = value(P);
R = value(R);
S = value(S);

end