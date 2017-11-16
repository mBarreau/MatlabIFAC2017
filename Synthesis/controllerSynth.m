%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Synthesys of a feedback for a system with an average on the state   %%
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
warning('off','YALMIP:strict')

% System def
A = [0.2 0; 0.2 0.1];
B = [-1 0; -1 -1];

% Parameters
alphaMax = 2;
alphaStart = 0;
alphaStep = 0.5;
nbDecimal = 4; % Precision of the answer
hStop = 10; % Maximal admissible delay h

% Rules:
e1 = 1;
e2 = 0;
e3 = 1;
e4 = 1;

% Storage matrix H
if(alphaStep > 0)
    H = zeros(2,floor(alphaMax/alphaStep)+1);
else
    H = zeros(2,1);
end

% Start of the algorithm
j = 0;
for alpha = alphaStart:alphaStep:alphaMax
    h = 10^(-nbDecimal);
    epsilon = 1;
    
    hmin = -1;
    hmax = -1;
    i = 0;
    while (1)
        [pres,K] = controller(A, B, h, e1, e2, e3, e4, alpha);
        if(sum(pres > 0) == length(pres)) % There exists a controler for this (h,alpha)
            h
            if hmin < 0 % hmin has never been defined previously
                hmin = h;
                display('hmin found');
            end
        else % There is no controller
            if hmin >= 0 % it was feasible for a smaller h
                h = h-epsilon; % We reduce the h to a feasible value
                if i == nbDecimal % If we get the maximum number of decimal, we stop
                    display('hmax found');
                    hmax = h;
                    [pres,K,gamma] = controller(A, B, hmax, e1, e2, e3, e4, alpha);
                    break
                end
                i = i + 1;
                epsilon = epsilon*0.1; % We increase the precision
            else
                h
            end
        end
        h = h + epsilon; % We increase the h for the new loop
        if h > hStop
            display('Stop before the end. hStop reached.');
            alpha = alphaMax;
            break;
        end
    end
    
    if h > hStop
        break;
    end
    j=j+1
    H(:, j) = [alpha; hmax];
    display('----------------------');
    hmax
    
end

display('----------------------');
hmin
hmax
