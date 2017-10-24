function output = check_fun_eval_c()

% Evaluate c1 and c2, given k, and check it against the dispersion relation
% (c) Matthew Juniper 2013

clear
clc

% Input the parameters
param.L = -0.9;
param.S = 1;
param.Sig = 1.3;
param.h = 2.32;

% Input the value of k
k = complex(2*rand(1),0);

% Create a cell array containing the names of the disp. rels to check
scheme{1} = 'single';
scheme{2} = 'unconf_jetwake_sin';
scheme{3} = 'unconf_jetwake_var';
scheme{4} = 'conf_jetwake_sin';
scheme{5} = 'conf_jetwake_var';

% Cycle through the dispersion relations to check them
for nn = 1:length(scheme)
  param.scheme = scheme{nn};
  % Calculate omega1 and omega2 (w1 and w2 for short)
  [c1,c2] = fun_eval_c(param,k);
  w1 = c1*k;
  w2 = c2*k;
  % Calculate D1 and D2
  D1 = fun_eval_D(param,k,w1);
  D2 = fun_eval_D(param,k,w2);
  % Check that both are zero to machine precision
  if (abs(D1)<10*eps)&&(abs(D2)<10*eps)
    disp(['Correct dispersion relation for ',param.scheme])
  else
    disp(['Incorrect dispersion relation for ',param.scheme])
    beep
  end
end

output = param;

end
