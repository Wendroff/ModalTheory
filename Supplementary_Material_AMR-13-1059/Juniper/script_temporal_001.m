function output = script_temporal_001()

% Perform a temporal stability anaylsis
% (c) Matthew Juniper 2013

% Set the dispersion relation
param.scheme = 'single';
% Set the parameters
param.L = 0.9;
param.S = 1;
% Set the sampling values of the wavenumber, k
k = (0:0.01:4)';
% Evaluate the phase velocity
[w1,w2] = fun_eval_w(param,k);

% Plot the results on the same chart
figure(1);clf
subplot(2,1,1)
cla; hold on; set(gca,'FontSize',20,'FontName','Times')
plot(k,imag(w1),'k-')
plot(k,imag(w2),'k--')
xlabel('wavenumber, $k$','Interpreter','Latex')
ylabel('growth rate, $\omega_i$','Interpreter','Latex')
title(['Single shear layer; ',...
  '$\Lambda = ',num2str(param.L), ...
  '$, $S = ',num2str(param.S),'$'],'Interpreter','Latex')
legend('First eigenvalue','Second eigenvalue')
%
subplot(2,1,2)
cla; hold on; set(gca,'FontSize',20,'FontName','Times')
plot(k,real(w1),'k-')
plot(k,real(w2),'k--')
xlabel('wavenumber, $k$','Interpreter','Latex')
ylabel('frequency, $\omega_r$','Interpreter','Latex')
legend('First eigenvalue','Second eigenvalue',2)
%
% print -deps 'script_temporal_001.eps'

output = {};

end
