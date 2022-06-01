clc
clear
close all

%Testing multiple iterations of fitting by trying to fit to a noisy,
%shifted michelson

%Make a sample michelson with sim_full
c_const = 299792458;
w_p = 2.*pi.*c_const./406e-9;
d_tau = (1./c_const).*50e-9; %step size of stage
num_steps = ceil(40e-15./d_tau);
x = (-num_steps.*d_tau:d_tau:num_steps.*d_tau)';
xerr = d_tau./16.*randn(length(x),1); %adds noise to x-grid (simulating error in stage)
x = x+xerr;

%Fine grid for plotting fitting michelson
d_tau_fine = (1./c_const).*5e-9;
num_steps_fine = ceil(40e-15./d_tau_fine);
x_fine = (-num_steps_fine.*d_tau_fine:d_tau_fine:num_steps_fine.*d_tau_fine)';

w_spdc_mu = 2.5455e15; %gives michelson similar to paper
w_spdc_sig = 1.0247e14;
fs_offset = -0.1e-15; %add small offset to time-zero to throw off fitting

scaled_mu = w_spdc_mu./1e15; %scale down mu, sig, offset
scaled_sig = w_spdc_sig./1e14;
scaled_offset = fs_offset./1e-15;

y = 1.1.*sim_full(x, scaled_mu, scaled_sig, scaled_offset,1); %scale by 1.1 to throw off fitting
y = y + 0.05.*randn(length(y),1); %add noise to y data
y = y + 0.1; %shift data in y direction to throw off fitting

%Note that all fitted variables are scaled to be close to one - this makes
%it so that the fitting algorithm puts them all on the same footing, and
%that the minimum/maximum step sizes are tuned correctly as well

%Passing a variable to a function will cause it to get fitted. I found that
%I can't fit everything at once, need to break it up into a few iterations.
%Use globals to pass variables to a function without being treated as
%inputs - if you have a better way please let me know, using globals seems
%pretty clunky and easy to cause a mistake
global mu_guess;
mu_guess = scaled_mu.*0.95;
global sig_guess;
sig_guess = scaled_sig.*0.8;
global offset_guess;
offset_guess = 0;
global scale_guess;
scale_guess = 1;
global shift_guess;
shift_guess = 0;
%Made all of the initial guesses be off from the simulation 

%see matlab docs on fitting to a function file:
%https://www.mathworks.com/help/curvefit/fittype.html#btpaend-1-expression
ft_offset = fittype("fit_offset(x,offset)");
fo_offset = fitoptions(ft_offset);
fo_offset.StartPoint = offset_guess;
fo_offset.Upper = 3; %time-zero in data needs to be within 3fs of zero in x-axis
fo_offset.Lower = -3;
fo_offset.DiffMaxChange = 0.01;
fo_offset.DiffMinChange = 1e-8;

ft_mu = fittype("fit_mu_sig( x, mu, sig)");
fo_mu = fitoptions(ft_mu);
fo_mu.StartPoint = [mu_guess,sig_guess];
fo_mu.Upper = [w_p./1e15, 100]; %pump freq is upper bound
fo_mu.Lower = [w_p./2./1e15, 0.01]; %pump freq/2 
fo_mu.DiffMaxChange = 0.01;
fo_mu.DiffMinChange = 1e-8;

ft_ss = fittype("fit_scale_shift( x, scale, shift)");
fo_ss = fitoptions(ft_ss);
fo_ss.StartPoint = [scale_guess, shift_guess];
fo_ss.Upper = [3,1];
fo_ss.Lower = [0.3,-1];
fo_ss.DiffMaxChange = 0.01;
fo_ss.DiffMinChange = 1e-8;

sse = [];

for i = 1:5 %5 iterations of fitting - seems to converge after 2 or 3
    fo_offset.StartPoint = offset_guess; %update start point
    [f_offset, gof_offset, output_offset] = fit(x,y,ft_offset,fo_offset) %fit offset
    offset_guess = f_offset.offset; %update offset guess
    [y_out,norm_coarse] = sim_full(x,mu_guess,sig_guess,offset_guess,scale_guess); %get coarse scaling
    [y_out,norm_fine] = sim_full(x_fine,mu_guess,sig_guess,offset_guess,scale_guess); %get fitted michelson
    y_out = y_out.*norm_fine./norm_coarse; %scale fine michelson to match coarse michelson
    y_out = y_out + shift_guess; %shift fine michelson
    sse = [sse,gof_offset.sse];

    figure
    title("offset "+num2str(i))
    plot(x.*1e15, y, '.', x_fine.*1e15, y_out);
    legend('raw','fit')

    fo_mu.StartPoint = [mu_guess, sig_guess];
    [f_mu, gof_mu, output_mu] = fit(x,y,ft_mu,fo_mu)
    mu_guess = f_mu.mu;
    sig_guess = f_mu.sig;
    [y_out,norm_coarse] = sim_full(x,mu_guess,sig_guess,offset_guess,scale_guess);
    [y_out,norm_fine] = sim_full(x_fine,mu_guess,sig_guess,offset_guess,scale_guess);
    y_out = y_out.*norm_fine./norm_coarse;
    y_out = y_out + shift_guess;
    sse = [sse, gof_mu.sse];

    figure
    title("mu "+num2str(i))
    plot(x.*1e15, y, '.', x_fine.*1e15, y_out);
    legend('raw','fit')

    fo_ss.StartPoint = [scale_guess, shift_guess];
    [f_ss, gof_ss, output_ss] = fit(x,y,ft_ss,fo_ss)
    shift_guess = f_ss.shift;
    scale_guess = f_ss.scale;
    [y_out,norm_coarse] = sim_full(x,mu_guess,sig_guess,offset_guess,scale_guess);
    [y_out,norm_fine] = sim_full(x_fine,mu_guess,sig_guess,offset_guess,scale_guess);
    y_out = y_out.*norm_fine./norm_coarse;
    y_out = y_out + shift_guess;
    sse = [sse, gof_ss.sse];

    figure
    title("ss "+num2str(i))
    plot(x.*1e15, y, '.', x_fine.*1e15, y_out);
    legend('raw','fit')    
end

figure
semilogy((1:length(sse))./3,sse)
xlabel('iteration number')
ylabel('SSE')

%take fourier transform of raw data and compare to FT of fitted data
f_y = 1./d_tau;
T_y = d_tau;
L_y = length(x);
ft_y = fft(y);
P2_y = abs(ft_y./L_y);
P1_y = P2_y(1:floor(L_y/2)+1);
P1_y(2:end-1) = 2.*P1_y(2:end-1);
f_axis = f_y.*(0:(L_y/2))./L_y;
f_axis = f_axis.*2.*pi; %convert to rad/s
figure
plot(f_axis, P1_y)
xlim([0,w_p.*1.3])

f_y_fine = 1./d_tau_fine;
T_y_fine = d_tau_fine;
L_y_fine = length(x_fine);
ft_y_fine = fft(y_out);
P2_y_fine = abs(ft_y_fine./L_y_fine);
P1_y_fine = P2_y_fine(1:floor(L_y_fine/2)+1);
P1_y_fine(2:end-1) = 2.*P1_y_fine(2:end-1);
f_axis_fine = f_y_fine.*(0:(L_y_fine/2))./L_y_fine;
f_axis_fine = f_axis_fine.*2.*pi;
hold on
plot(f_axis_fine, P1_y_fine);
xlim([0,w_p.*1.3]);
legend('raw','fit')
