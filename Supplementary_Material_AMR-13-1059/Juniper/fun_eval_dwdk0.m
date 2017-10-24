function [k,w] = fun_eval_dwdk0(param,k_init,tol)

% Find the value of k that gives dw/dk = 0, given k_init
% (c) Matthew Juniper 2013

% Set the options for fsolve
options = optimset('Display','Off','TolFun',tol);

% Converge to the exact value of k using Matlab's fsolve
k = fsolve(@myfun,k_init,options);

% Evaluate omega at this point
[dwdk01,dwdk02] = fun_eval_dwdk(param,k);
[w01,w02] = fun_eval_w(param,k);
[~,i] = min(abs([dwdk01,dwdk02]));
w = (i==1)*w01 + (i==2)*w02;

% Embedded function to evaluate dwdk
  function [dwdk] = myfun(k)
    [dwdk1,dwdk2] = fun_eval_dwdk(param,k);
    dwdk = min(dwdk1,dwdk2);
  end
% End of embedded function

end % end function
