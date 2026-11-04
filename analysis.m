clc, clearvars, close all

%% Osbservability Analysis of KF
u_d = 7;
[T,K] = nomoto_first_order(u_d);
A = [   
    0 1 0;
    -K/T 0 -K/T;
    0 0 0;
    ];
B = [0 K/T 0 ]';
E = [
    0 0 0;
    0 1 0;
    0 0 1;
    ];
C = [1 0 0];
K = obsv(A,C); %kalman obs matrix
observabilityRank = rank(K);
isObservable = observabilityRank == size(A, 1);
disp(observabilityRank)
disp(isObservable)

%% Discretization
h = 0.1;                % Sampling time (s)
[Ad,Bd] = c2d(A,B,h);
[Ad,Ed] = c2d(A,E,h);
disp(Ad);
disp(Bd);
disp(Ed);
