% CADiS configuration coarrays with different compressed factor p &
% displaced factor l
% reference:
% 1. Qin S, Zhang Y D, Amin M G. Generalized 
% coprime array configurations for direction-of-arrival estimation[J]. 
% Signal Processing, IEEE Transactions on, 2015, 63(6): 1377-1390.
close all; clear; clc;

%% initialization
f = 3400;
c = 340;
lamda = c/f;
d = lamda/2;
phi = [-90:0.1:90];
M = 6;
N = 7;
p = 2;               % compressed fator
L = M/p+N;           % displaced factor
%% coordinate
p_N = zeros(2,N);    % storing the 2-d coordinates of subarray sensors
p_M = zeros(2,M-1);

% the positions of subarray sensors
p_N(2,:) = 0:M/p*d:M*(N-1)*d/p;
p_M(2,:) = M*(N-1)/p*d+L*d:N*d:(M-2)*N*d+M*(N-1)/p*d+L*d;

% the positions of self-difference
l_N = zeros(1,N);
l_M = zeros(1,M);
l_N = 0:M/p:M*(N-1)/p;
l_M = 0:N:N*(M-1);
L_s = union(l_N,l_M);

% the positions of cross-difference
k=1;
for i=1:M
    for j=1:N                      
       L_c(1,k) = N*(i-1)-M/p*(j-1);
       L_c1(1,k) = N*(i-1)-M*(j-1);
       k=k+1;
    end
end
L_P = union(union(L_s, -1*L_s), union(L_c, -1*L_c));

[~,n] = size(L_P);
p = zeros(2,n);
p(2,:) = L_P*d;

%% the positions of virtul sensors
figure ;
plot(L_P,0,'r*');
hold on;
axis([-37 37 -1 1]);
grid minor;


%% the positions of array sensors
figure;
plot(p_N(2,:),0,'ro');
hold on;
plot(p_M(2,:),0,'b*');
title('CADiS(M=6,N=7,p=2)');
axis([0 3 -0.5 0.5]);
grid on;
