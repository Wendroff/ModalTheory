function output = script_spatial_001()

% Perform a spatial stability anaylsis on a single shear layer
% (c) Matthew Juniper 2013

clear
clc

% Input the parameters
param.scheme = 'single';
param.L = 0.5;
param.S = 1;

% Set the sampling values of the angular frequency, w
w = [0.01:0.01:1]';
% Input an initial guess for k
k_init = complex(0.1,-0.02);
% Input the tolerance for convergence to the final value of k
tol = 10^-6;
% Initialize k vector
k = 0*w;

% Loop through all the values of omega, finding the relevant value of k
for nn = 1:length(w)
  k(nn) = fun_eval_k(param,k_init,w(nn),tol);
end

% Plot the results on the same chart
figure(1);clf
subplot(2,2,1)
title(['Spatial; '...
  '   $\Lambda = ',num2str(param.L), ...
  '$, $S = ',num2str(param.S), ...
  '$'],'Interpreter','Latex')
%
cla; hold on; set(gca,'FontSize',14,'FontName','Times')
plot(real(w),imag(k),'k')
ylabel('Imag($k$)','Interpreter','Latex')
axis([0 1 -0.15 0.05])
%
subplot(2,2,3)
cla; hold on; set(gca,'FontSize',14,'FontName','Times')
plot(real(w),real(k),'k')
xlabel('angular frequency, $\omega$','Interpreter','Latex')
ylabel('Real($k$)','Interpreter','Latex')
axis([0 1 0 1])

% print -deps 'script_spatial_001.eps'

output = {};

end
