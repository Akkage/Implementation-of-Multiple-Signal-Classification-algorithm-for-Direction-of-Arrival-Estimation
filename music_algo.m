clc; clear; close all;
%% Project : Implementation of MUSIC (Multiple Signal Classification) algorithm for Direction of Arrival Estimation(DOA).
% Name     : Akshay Pramod Age
% Email Id : akshayage27@gmail.com
% Type     : Self Project
% Concept  : Based on Linear Algebra

%%
L     = 20;   % # antennas
P     = 5;    % # sources
M     = 100;  % # samples
fc    = 10^6; % centre frequency
fs    = 10^7; % sampling frequency
c     = 3*10^8; % speed of light
x_var = 1; % transmitted symbol variances
n_var = 1; % noise variance
doa   = [20 50 85 110 145]; %randi([0 180], 1,P); % direction of arrival
d     = 150; % antenna spacing
theta = 0:1:180;
%%
x = zeros(P,M); % narrowband complex exponential source signals
for p1 = 1:P
    x(p1,:) = sqrt(x_var)*randn(1,M).*exp(1i*2*pi*fc*[1:M]/fs);
end
A = zeros(L,P); % matrix of steering vectors
for p1 = 1:P
    A(:,p1) = exp(-1i*2*pi*d*(fc/c)*[0:L-1]'*cosd(doa(p1)));
end

% signal received at the antenna array
y = A*x + sqrt(n_var)*randn(L,M);

% Covariance matrix of antenna array
R = (y*y')/M;

% Eigen value decomposition
[V,D] = eig(R);

% Subspace corresponding to small eigen values
U = V(:,1:(L-P));

% 
fn  = zeros(length(theta),1);
A_d = zeros(L,length(theta));

for th = 1 : length(theta)
    A_d(:,th) = exp(-1i*2*pi*d*(fc/c)*[0:L-1]'*cosd(theta(th)));
    fn(th,1)    = 1/(norm(A_d(:,th)'*U)^2);
end

% The peaks in the graph are the estimated DOAs
% Compare the peaks in graph with the values from the vector doa
plot(theta,fn)
xlabel('Angle (deg)');
ylabel('1/Norm^2');
title ('Estimation of DOAs ')
% findpeaks(fn)
