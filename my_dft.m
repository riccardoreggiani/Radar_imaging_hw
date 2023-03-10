function [X,f] = my_dft(g,t, Nf)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


G = fft(g,Nf,1);
G = fftshift(G,1);

dt = t(2)-t(1);

% ATTENTION: there is a mistake in this point in the code we wrote in
% class! It was df = 1/(dt/Nf); which is clearly wrong
df = 1/(dt*Nf); % fs/Nf

if mod(Nf,2) % Odd number of samples
    f = (-(Nf-1)/2 : (Nf-1)/2)*df;
else
    f = (-Nf/2:Nf/2-1)*df;
end

phase = -2*pi*f(:)*t(1);
w = exp(1j.*phase);

X = G.*w*dt;


end