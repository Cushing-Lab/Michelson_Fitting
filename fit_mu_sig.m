function y = fit_mu_sig(x,mu,sig)
%Generates Michelson given by nondegenerate SPDC spectrum with mean (mu) and standard deviation (sig) on a grid (x)
%Uses globals to define Michelson scaling and offset

%scale mu and sig back up to rad/s
mu = 1e15.*mu;
sig = 1e14.*sig;

global offset_guess;
global shift_guess;
global scale_guess;

%scale offset back to fs
offset = offset_guess.*1e-15;

x = x-offset;

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
y = shift_guess + y;
end