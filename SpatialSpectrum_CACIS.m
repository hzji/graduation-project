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
d = c/f/2;                          % unit spacing = half-wavelength
T = 2000;                              % snapshots
N_Sig = 33;                            % the num of input signals
M = 6;                                 % number of sensor for subarray A 
N = 7;                                 % number of sensor for subarray B
DOA = linspace(-60, 60, N_Sig);        % direction of input signals
p = 6;                                 % compresssion factor p = 2,3,6
SNR = 10;                               % 信噪比

Fc=[2*10^3:2*10^3/(N_Sig-1):5*10^3];   % frequencies of input sig

%% ----------------------Generate signals-------------------------
T_Vector=(1:T)/f;  
p_N = [0:M/p:M*(N-1)/p];   
p_M = [0:N:(M-1)*N];
P = union(p_N,p_M);                    % the positions of physical sensors

A = zeros(length(P),N_Sig);            % manifold    
SigVec = zeros(N_Sig,T);

%% Generate the imping Sources
for Q = 1:N_Sig        
        A(:,Q) = exp(-j*P'*2*pi*d*sin(DOA(Q)*pi/180)*f/c); 
        SigVec(Q,:) = exp(1j*2*pi*Fc(Q).*T_Vector);  
end
x0 = A*SigVec; 
% Generate noises  on Signals
x = awgn(x0,SNR,'measured');          % Received data

%% construct covariance matrix                                                                  
R = x*x'/T;                                                            
z = R(:);       % vectorization: matrix--->vector
z1 = CACIS_Sort(z,P,M,N,p);

MM = M*N-M*(N-1)/p;
Ri = zeros(MM,MM,MM);
for i = 1:MM
    zi = z1(i:i+MM-1);
    Ri(:,:,i) = zi*zi';
end
Rz = sum(Ri,3)/MM;

%% ---------------------MUSIC algorithm---------------------------
[U,S,V] = svd(Rz);
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
manifold_scan = exp( -2j*pi*f/c*d* [0:MM-1].'*sin(az_scan) );
BeamPattern = 1./diag( manifold_scan'*(Un*Un')*manifold_scan );
BeamPattern_norm = abs(BeamPattern)/max(abs(BeamPattern));
BeamPattern_db = 10*log10(BeamPattern_norm);

figure; 
plot(az_scan/pi*180, BeamPattern_norm);
xlabel('\theta(deg)');
ylabel('Spectrum(Normalized)');
