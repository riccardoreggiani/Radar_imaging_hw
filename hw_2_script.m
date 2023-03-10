clc;
clear all;
close all;

%PROBLEM DATA
h = 693e3;                               %hight of the satellite
b = 200;                                 %baseline
f = 5.4e9;                               %frequency
SRR = 5;                                 %slant range resolution
teta_i = ((35*pi) / 180);                %incident angle in rad
teta_d = pi/2 - teta_i;                  %depression angle (calculation from the model)
c = physconst('lightspeed');             %speed of light
lambda = c/f;                            %wavelength
rho_rg = SRR/4;                          %sampling in ground range direction

%Master Satellite Position
Y_m = 0;                                 % y axis is the horizontal 
Z_m = h;  
%Slave Satellite Position 
Y_s = 0+b; 
Z_s = h;
%Generate ground range
p = h*tan(teta_i);                       %middle point of mountain (Pitagorean theorem)
y_i = (-500:rho_rg:500)+p;               %spanning 1km along mid point p;

%Generate complex reflectivity
a_t = rand(1,801);                       %abs
f_t = 2*pi.*rand(1,801);                 %phase
t_i = a_t.*exp(1i*f_t);
                                     

%Generate elevation profile
z_i = 6000*normpdf(y_i,p-20/10,100);
figure
plot(y_i,z_i)
title('Elevation profile')
ylabel('Heigth [m]')
xlabel('Slant range [m]')
%Generate slant range axis
sr = h/cos(teta_i);
sr_axis = sr-2000:SRR/10:sr+2000;
%%
%Distance from n_th sensor to target pi=[yi,zi,ti]
R_m = sqrt((Y_m-y_i).^2 + (Z_m-z_i).^2);                             %for MASTER
R_s = sqrt((Y_s-y_i).^2 + (Z_s-z_i).^2);                             %for SLAVE

% inizialization of a matrix with length = length(sr_axis)
I_m = zeros(length(sr_axis),1);
I_s = zeros(length(sr_axis),1);

for ii = 1:length(sr_axis)
    % formula (2) 
    I_m(ii,1) = sum(t_i.*sinc((sr_axis(ii)-R_m)/SRR).*exp(-1i*4*pi*R_m/lambda));
    I_s(ii,1) = sum(t_i.*sinc((sr_axis(ii)-R_s)/SRR).*exp(-1i*4*pi*R_s/lambda));
end

figure;
subplot(2,1,1);
plot(sr_axis,abs(I_m));
ylabel('Height [m]');
xlabel('Slant Range [m]');grid on;grid minor;
title('Master')
subplot(2,1,2)
plot(sr_axis,abs(I_s));
ylabel('Height [m]');
xlabel('Slant Range [m]');grid on;grid minor;
title('Slave')

%% COREGISTRATION

%Coord. of reference surface
z_iref = z_i.*0;
y_iref = y_i;

%formula (4)
R_mref = sqrt((Y_m-y_iref).^2 + (Z_m-z_iref).^2);                    %distance from the master satellite to ref surface

z_r = interp1(R_mref,z_iref,sr_axis,'linear');
y_r = interp1(R_mref,y_iref,sr_axis,'linear');

r_mref = sqrt((Y_m-y_r).^2+(Z_m-z_r).^2);                            %in radar coordinates (slant range)
r_sref = sqrt((Y_s-y_r).^2+(Z_s-z_r).^2);

%formula (6)
delta_r1 = r_mref-r_mref;
delta_r2 = r_sref-r_mref;

I_n_c_r1 = interp1(sr_axis,I_m,sr_axis+delta_r1,'linear');
I_n_c_r2 = interp1(sr_axis,I_s,sr_axis+delta_r2,'linear');

figure;
subplot(2,1,1);
plot(sr_axis,abs(I_n_c_r1));
ylabel('Height [m]');
xlabel('Slant Range [m]');grid on;grid minor;
axis([8.45e5 8.47e5 0 2]);
title('Coregistration in slant range')
                                                                     
subplot(2,1,2);
plot(sr_axis,abs(I_n_c_r2));
ylabel('Height [m]');
xlabel('Slant Range [m]');grid on;grid minor;
axis([8.45e5 8.47e5 0 2]);

%% INTERFEROMETRIC
%Point 4.1

I = I_n_c_r1.* conj(I_n_c_r2);                       %interferogram
figure
plot(sr_axis,angle(I));
title('Interferometry')
xlabel('Slant range [m]');
ylabel('Phase [rad]')
