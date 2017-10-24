function [dwdk1,dwdk2] = fun_eval_dwdk(param,k)

% Evaluate dw1/dk and dw2/dk by finite difference, given k
% (c) Matthew Juniper 2013
%
% (n.b. this version will break if k is on a branch cut)

% Input incremental k
delta = 1e-6;

% Calculate w1 and w2 at these values of k
[w11,w21] = fun_eval_w(param,k);
[w12,w22] = fun_eval_w(param,k+delta);

% Calculate the gradients
dwdk1 = (w12-w11)/delta;
dwdk2 = (w22-w21)/delta;

end % end function
