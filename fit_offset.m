function y = fit_offset(x,offset)
%Shifts time-zero of Michelson with (offset) on the grid (x)
%Uses globals to define Michelson shape and scaling
%not good at finding global minimum, need to start within one cycle of the Michelson time-zero fringe

offset = 1e-15.*offset; %scale offset back to fs

global mu_guess;
global sig_guess;
global scale_guess;
global shift_guess;

%scale mu and sig back up to rad/s
mu = 1e15.*mu_guess;
sig = 1e14.*sig_guess;

x=x-offset; %shift x-axis to account for offset in time-zero

c = 299792458;
w_p = 2.*pi.*c./406e-9;
w_mesh_u = mu+4.*sig;
w_mesh_l = w_p - w_mesh_u;
w_mesh = linspace(w_mesh_l,w_mesh_u,1001); %mesh of angular freqs defining SPDC spectral density, centered around 2*pump frequency
w_mesh = permute(w_mesh, [3 1 2]); %need to do a 3d array for matlab's fitting algorithm to work for whatever reason
S = (1./(sig.*sqrt(2.*pi))).*exp(-(w_mesh-mu).^2./(2.*sig.^2)); %standard gaussian
S = S + flip(S); %mirror SPDC about 2*pump frequency
g_ff = (1+exp(1i.*w_p.*x)+exp(1i.*w_mesh.*x)+exp(1i.*(w_p-w_mesh).*x)); %coincidence transfer function
integrand = S.*g_ff.*conj(g_ff);
P_ff = real(sum(integrand,3)); %Perform integration as riemann sum (don't worry about absolute value since it's all relative)
y = scale_guess.*P_ff./max(P_ff);
y = y + shift_guess;
end