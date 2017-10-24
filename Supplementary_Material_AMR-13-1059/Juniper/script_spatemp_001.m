function output = script_spatemp_001()

% Perform a spatio-temporal stability anaylsis
% (c) Matthew Juniper 2013

param.scheme = 'single';

% Input the parameters
param.L = 0.9;
param.S = 1;
% param.Sig = 1;
% param.h = 1;

% Set the sampling values of the wavenumber, k
[kr,ki] = meshgrid(0:0.02:2,-1:0.02:1);
k = complex(kr,ki);

% Calculate omega1 and omega2 at these values of k
[w1,w2] = fun_eval_w(param,k);

% Contour plot theses side by side
figure(1);clf
% Plot the real and imaginary parts of omega1
subplot(2,2,1)
cla; hold all; set(gca,'FontSize',14,'FontName','Times')
colormap(gray)
xlabel('$k_r$','Interpreter','Latex')
ylabel('$k_i$','Interpreter','Latex')
title('$\omega_1(k)$','Interpreter','Latex')
contourf(kr,ki,imag(w1),20)
contour(kr,ki,imag(w1),[0 0],'LineColor','k','LineWidth',2)
contour(kr,ki,real(w1),20,'LineColor',[0.6 0.6 0.6])
axis image
colorbar

% Plot the real and imaginary parts of omega2
subplot(2,2,2)
cla; hold all; set(gca,'FontSize',14,'FontName','Times')
xlabel('$k_r$','Interpreter','Latex')
ylabel('$k_i$','Interpreter','Latex')
title('$\omega_2(k)$','Interpreter','Latex')
contourf(kr,ki,imag(w2),20)
contour(kr,ki,imag(w2),[0 0],'LineColor','k','LineWidth',2)
contour(kr,ki,real(w2),20,'LineColor',[0.6 0.6 0.6])
axis image
colorbar

subplot(2,2,4)
surf(kr,ki,imag(w2))

subplot(2,2,3)
surf(kr,ki,imag(w1))
print -depsc2 'script_spatemp_001.eps'

output = {};

end
