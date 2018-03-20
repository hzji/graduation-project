% CACIS configuration coarrays, for different compression factor p
% M = 6, N = 7; p = 2,3,6 

close all; clear; clc;
%% 参数
phi0 = 0;
phi = -90:0.1:90;
M = 6;
N = 7;
p = 6;               %压缩系数
f = 2000;
c = 1500;
lamda = c/f;
d = lamda/2;


%% 坐标
p_N = zeros(2,N);    %存储子阵列二维坐标，这里x维均为0
p_M = zeros(2,M);

p_N1 = zeros(2,N);   %存储子阵列二维坐标，这里x维均为0
p_M1 = zeros(2,M);

% 子阵列阵元位置
p_N(2,:) = 0:M/p*d:M*(N-1)*d/p;
p_M(2,:) = 0:N*d:(M-1)*N*d;

p_N1(2,:) = 0:M*d:M*(N-1)*d;
p_M1(2,:) = 0:N*d:(M-1)*N*d;

% p = zeros(2,M*N);
% p(2,:) = 0:d:(M*N-1)*d;
%% 自时延正轴
%CACIS
l_N = zeros(1,N);
l_M = zeros(1,M);
l_N = 0:M/p:M*(N-1)/p;
l_M = 0:N:N*(M-1);
L_s = union(l_N,l_M);
% 原型互质阵列
l_N1 = zeros(1,N);
l_M1 = zeros(1,M);
l_N1 = 0:M:M*(N-1);
l_M1 = 0:N:N*(M-1);
L_s1 = union(l_N1,l_M1);
%% 交叉时延
k=1;
for i=1:M
    for j=1:N                      
       L_c(1,k) = N*(i-1)-M/p*(j-1);
       L_c1(1,k) = N*(i-1)-M*(j-1);
       k=k+1;
    end
end
L_P = union(union(L_s, -1*L_s), union(L_c, -1*L_c));
L_P1 = union(union(L_s1, -1*L_s1), union(L_c1, -1*L_c1));
[~,n] = size(L_P);
[~,n1] = size(L_P1);
p = zeros(2,n);
p1 = zeros(2,n1);
p(2,:) = L_P*d;
p1(2,:) = L_P1*d;

%% 虚拟阵元
figure ;
plot(L_P,0.1,'r*');
hold on;
plot(L_P1,-0.1,'b*');
hold on;
axis([-37 37 -1 1]);
grid minor;


%% 物理阵元
figure;
plot(p_N(2,:),0.1,'ro');
hold on;
plot(p_M(2,:),0.1,'b*');
% hold on;
% plot(p_N1(2,:),-0.1,'bo');
% hold on;
% plot(p_M1(2,:),-0.1,'b*');
axis([0 15 -0.5 0.5]);
grid on;
