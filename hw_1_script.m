clear all;
close all;
clc;

%% PARAMETERS OF THE HOMEWORK
N = 32; % number of antenna

c=3*10^8; % velocity of light
f = 77*10^9; % frequency
Band = 2*10^9; % Bandwidth
rho_R = Band*2/c; 

lambda = c/f; % wavelength

dx = lambda / 4; % antenna spacing
max_res = 3.6; %degrees

x_s = (0:N-1)*dx; %x_s of the sensor
x_s = x_s - mean(x_s); % trick to obtain vector of sensors
y_s = zeros(size(x_s)); 

% target position
x_t = 30;
y_t = 50;

wid = .1; %rect amplitude
r_ax = 0:wid/10:300;

%% Received compressed signal

%data matrix for received signal
a = zeros(length(r_ax), N);
for ii = 1:N
    R = sqrt((x_s(ii)-x_t).^2 + (y_s(ii)-y_t).^2); %pitagora theorem
    
    a(:,ii) = sinc((r_ax-R)/rho_R)*exp(-1j.*4*pi*R/lambda);
end

figure; imagesc(1:N, r_ax, abs(a));
xlabel("Sensors"); ylabel("Distance [m]")
title('Received Compressed Signal')
%% DoA of the target estimation 

a=a';
A=my_dft(a,r_ax,1024);
y = -90:1:90;
figure; imagesc(r_ax,y,abs(A));
colorbar
title('DOA estimation of the target')
xlabel('range [m]')
ylabel('angle [deg]')
%% Placing the image target
%placing the target distant from an obstacle
x_i = -x_t -20;
y_i = y_t;

%% DoA of the image target

%range-compress data matrix for image target
b = zeros(length(r_ax), N);
for ii = 1:N

    R_i = sqrt((x_s(ii)-x_i).^2 + (y_s(ii)-y_i).^2); %pitagora theorem
    
    b(:,ii) = sinc((r_ax-R_i)/rho_R)*exp(-1j.*4*pi*R_i/lambda);

end

b=b';
B=my_dft(b,r_ax,1024);

figure; imagesc(r_ax,y,abs(B)); % from this we retrive f ---> substitute with 2/lambda * sin(teta)
title('DOA of the image target')
ylabel('angle [deg]')
xlabel('range [m]')
%% Ghost

% ghost target obtained by multipath
g = zeros(length(r_ax), N);
for ii = 1:N

    R = sqrt((x_s(ii)-x_t).^2 + (y_s(ii)-y_t).^2); %pitagora theorem
    R_i = sqrt((x_s(ii)-x_i).^2 + (y_s(ii)-y_i).^2); %pitagora theorem
    R_g = R + R_i; 
    
    
    g(:,ii) = sinc((r_ax-R_g)/rho_R)*exp((-1j.*2*pi*R_g)/lambda);

end

g = g';
G = my_dft(g,r_ax,1024);
figure; imagesc(r_ax,y,abs(G));
title('Ghost Target');
ylabel('angle [deg]')
xlabel('range [m]')
%% Bistatic
clear all;

%% PARAMETERS OF THE HOMEWORK
N = 32; % number of antenna
y = -90:1:90;
c=3*10^8;
f = 77*10^9;
Band = 2*10^9;
rho_R = Band*2/c;

lambda = c/f; % wavelength

dx = lambda / 2; % antenna spacing
max_res = 3.6; %degrees

x_s = (0:N-1)*dx; %x_s of the sensor
x_s = x_s - mean(x_s);
y_s = zeros(size(x_s));


x_t = 30;
y_t = 50;

%campionamento
wid = .1; %rect amplitude
r_ax = 0:wid/10:300;

%% Real Target

a = zeros(length(r_ax), N);
for ii = 1:N
    R = sqrt((x_s(ii)-x_t).^2 + (y_s(ii)-y_t).^2); %pitagora theorem
    
    a(:,ii) = sinc((r_ax-R)/rho_R)*exp(-1j.*2*pi*R/lambda);
end

figure; imagesc(1:N, r_ax, abs(a));
xlabel("Sensor"); ylabel("Distance [m]")
title('Received compressed signal (Bistatic)')

%% DoA of the target

a=a';
%
A=my_dft(a,r_ax,1024);

figure; imagesc(r_ax,y,abs(A)); 

title ('DOA of target Bistatic ');
xlabel('range [m]')
ylabel('angle [deg]')
%% Placing the image target
% it is important to underline that between that the vehicle should be
% distant from the guardrail or an obstacle
x_i = -x_t -20;
y_i = y_t;

%% DoA of the image target

%range-compress data matrix for image target
b = zeros(length(r_ax), N);
for ii = 1:N

    R_i = sqrt((x_s(ii)-x_i).^2 + (y_s(ii)-y_i).^2); %pitagora theorem
   
    b(:,ii) = sinc((r_ax-R_i)/rho_R)*exp(-1j.*2*pi*R_i/lambda);

end

b=b';
B=my_dft(b,r_ax,1024);

figure; imagesc(r_ax,y,abs(B)); % from this we retrive f ---> substitute with 2/lambda * sin(teta)
title('DOA of image target Bistatic');
xlabel('range [m]')
ylabel('angle [deg]')

%% Ghost

g = zeros(length(r_ax), N);
for ii = 1:N

    R = sqrt((x_s(ii)-x_t).^2 + (y_s(ii)-y_t).^2); %pitagora theorem
    R_i = sqrt((x_s(ii)-x_i).^2 + (y_s(ii)-y_i).^2); %pitagora theorem
    R_g = R + R_i; 
    g(:,ii) = sinc((r_ax-R_g)/rho_R)*exp((-1j.*2*pi*R_g)/lambda); 
end

g = g';
G = my_dft(g,r_ax,1024);
figure; imagesc(r_ax,y,abs(G));
title('DOA of Ghost Target Bistatic');
ylabel('angle [deg]')
xlabel('range [m]')


