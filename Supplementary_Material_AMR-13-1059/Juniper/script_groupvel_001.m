function output = script_groupvel_001()

% Perform a temporal stability anaylsis
% and plot the phase and group velocities of the permitted waves.
% (c) Matthew Juniper 2013

% Set the dispersion relation
param.scheme = 'single';
% Set the parameters
param.L = 1;
param.S = 1;
param.Sig = 0.001;
param.h = 1;
% Set the sampling values of the wavenumber, k
k = [0:0.01:4]';
% Evaluate the phase velocity
[c1,c2] = fun_eval_c(param,k);
% Evaluate the group velocity
[dwdk1,dwdk2] = fun_eval_dwdk(param,k);

% Plot the results on the same chart
figure(1);clf
subplot(2,2,1)
title(['Phase and group velocities '...
  ' , $\Lambda = ',num2str(param.L), ...
  '$, $S = ',num2str(param.S), ...
  '$'],'Interpreter','Latex')
%
cla; hold on; set(gca,'FontSize',14,'FontName','Times')
plot(k,imag(c1),'r-')
plot(k,imag(c2),'b-')
plot(k,imag(dwdk1),'r--')
plot(k,imag(dwdk2),'b--')
xlabel('wavenumber, $k$','Interpreter','Latex')
ylabel('Imag($c_g$)','Interpreter','Latex')
set(gca,'YLim',[-4 4])
han = legend('$c_1$','$c_2$','$c_{g1}$','$c_{g2}$');
set(han,'Interpreter','Latex')
%
subplot(2,2,2)
cla; hold on; set(gca,'FontSize',14,'FontName','Times')
plot(k,real(c1),'r-')
plot(k,real(c2),'b-')
plot(k,real(dwdk1),'r--')
plot(k,real(dwdk2),'b--')
xlabel('wavenumber, $k$','Interpreter','Latex')
ylabel('Real($c_g$)','Interpreter','Latex')
set(gca,'YLim',[-4 4])
han = legend('$c_1$','$c_2$','$c_{g1}$','$c_{g2}$',3);
set(han,'Interpreter','Latex')
%
% print -deps 'script_groupvel_001.eps'

output = {};

end
