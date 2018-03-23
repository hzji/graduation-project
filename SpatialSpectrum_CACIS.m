% DOA estimation for CACIS config. using  MUSIC algorithm
% reference:
% 1. Qin S, Zhang Y D, Amin M G. Generalized 
% coprime array configurations for direction-of-arrival estimation[J]. 
% Signal Processing, IEEE Transactions on, 2015, 63(6): 1377-1390.
% 2. Pal P, Vaidyanathan P P. Coprime sampling and the music algorithm[C] 
% Digital Signal Processing Workshop and IEEE Signal Processing Education 
% Workshop. IEEE, 2011:289-294.
clc;clear all;close all;

%% ----------------------initialization---------------------------
c = 340;
f = 3400;
lambda = c/f;
d = lambda/2;                          % 阵元间距离，取为入射波长的一半
T = 2000;                              % snapshots
N_Sig = 30;                            % 入射信号数量
M = 6;                                 % number of sensor for subarray A 
N = 7;                                 % number of sensor for subarray B
DOA = linspace(-60, 60, N_Sig);        % 入射信号角度
p = 6;                                 % compressed factor p = 2,3,6
SNR = 15;                              % 信噪比

Fc=[2*10^3:2*10^3/(N_Sig-1):5*10^3];                 % 入射信号频率

thetatest=(-90*pi/180:1*pi/180:90*pi/180);               % theta角度搜索范围  
thetanum=length(thetatest);


%% ------------------------信号生成--------------------------------
T_Vector=(1:T)/f;  
p_N = [0:M/p:M*(N-1)/p];   
p_M = [0:N:(M-1)*N];
P = union(p_N,p_M)                                     % 物理阵元位置

A = zeros(length(P),N_Sig);          
SignalVector = zeros(N_Sig,T);

%% 构造非均匀线阵信号
for Q=1:N_Sig        
        A(:,Q)=exp(-j*P'*2*pi*d*sin(DOA(Q)*pi/180)/lambda); 
        SignalVector(Q,:)=exp(1j*2*pi*Fc(Q).*T_Vector);  
end
X0=A*SignalVector;
X=awgn(X0,SNR,'measured');                              % 信号加噪

%% ---------------------S_MUSIC Algorithm-------------------------
                                                                       
R=X*X'/T;                                                            
z=R(:);                                         % 矢量化 matrix--->vector

%% 提取连续虚拟阵列对应的z1   
% 虚拟阵列在-M*N+M*(N-1)/factor+1,M*N-M*(N-1)/factor-1上连续
% 抽出相应的行向量，先确定向量号
%flag = 1;           %不同互质结构的标识符,1表示CACIS结构,0表示CADiS
row_num = extract_row1(P,M,N,p);
z1 = z(row_num);

MM = M*N-M*(N-1)/p;
Ri = zeros(MM,MM,MM);
for i = 1:MM
    zi = z1(i:i+MM-1);
    Ri(:,:,i) = zi*zi';
end
Rzz = sum(Ri,3)/MM;

%% ---------------------MUSIC algorithm---------------------------
[U,S,V] = svd(Rzz);
lamda = diag(S);
ratio = lamda(1:end-1)./lamda(2:end);
if MM >= N_Sig+1
    ix = N_Sig+1;
else
    ix = find(ratio>10) + 1;  % find signal subspace & noise subspace
    ix = ix(1);
end
Un = U(:,ix:end);             % noise subspace

% Spatial spectrum
N_scan = 1801;
az_scan = linspace(-pi/2, pi/2, N_scan);
manifold_scan = exp( -2j*pi/lambda*d* [0:MM-1].'*sin(az_scan) );
BeamPattern = 1./diag( manifold_scan'*(Un*Un')*manifold_scan );
BeamPattern_norm = abs(BeamPattern)/max(abs(BeamPattern));
BeamPattern_db = 10*log10(BeamPattern_norm);

figure; 
plot(az_scan/pi*180, BeamPattern_norm);
xlabel('\theta(deg)');
ylabel('Spectrum(Normalized)');
