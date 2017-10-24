function output = script_spatial_002()

% Perform a spatial stability anaylsis on a single shear layer
% both directly and via Gaster's transformation.
% (c) Matthew Juniper 2013

clear
clc

% Input the parameters
param.scheme = 'single';
param.L = 0.5;
param.S = 1;

%%% TEMPORAL ANALYSIS
% Set the sampling values of the wavenumber, k
k = (0:0.01:2)';
% Evaluate the phase velocity
[w1,w2] = fun_eval_w(param,k);

%%% GASTER'S TRANSFORMATION
% Evaluate the group velocity
[dwdk1,dwdk2] = fun_eval_dwdk(param,k);
% Find kr and ki (approximately) through Gaster's transformation
krG = k;
kiG = -imag(w1)./(real(dwdk1));

%%% SPATIAL ANALYSIS (taken from script_spatial_001.m)
% Set the sampling values of the angular frequency, w
w = [0.01:0.01:1]';
% Input the tolerance for convergence to the final value of k
tol = 10^-6;
% Initialize k vector
k = 0*w;
% Loop through all the values of omega, finding the relevant value of k
for nn = 1:length(w)
  % Find the index of w1 that is closest to w
  ind = dsearchn(real(w1),w(nn));
  % Find the value of kiG at that index
  k_init = complex(krG(ind),kiG(ind));
  % Find the exact value of k from the spatial analysis
  k(nn) = fun_eval_k(param,k_init,w(nn),tol);
end

% Plot the results on the same chart
figure(1);clf
subplot(2,2,1)
title(['Spatial'...
  ' ; $\Lambda = ',num2str(param.L), ...
  '$, $S = ',num2str(param.S), ...
  '$'],'Interpreter','Latex')
%
cla; hold on; set(gca,'FontSize',14,'FontName','Times')
plot(real(w),imag(k),'k--')
plot(real(w1),kiG,'k')
ylabel('Imag($k$)','Interpreter','Latex')
legend('Exact','Gaster',4)
axis([0 1 -0.15 0.05])
%
subplot(2,2,3)
cla; hold on; set(gca,'FontSize',14,'FontName','Times')
plot(real(w),real(k),'k--')
plot(real(w1),krG,'k')
xlabel('angular frequency, $\omega$','Interpreter','Latex')
ylabel('Real($k$)','Interpreter','Latex')
legend('Exact','Gaster',2)
axis([0 1 0 1])

% print -deps 'script_spatial_002.eps'

output = {};

end
