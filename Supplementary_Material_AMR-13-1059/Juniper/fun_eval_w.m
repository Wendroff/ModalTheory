function [w1,w2] = fun_eval_w(param,k)

% Evaluate w1 and w2, given k
% (c) Matthew Juniper 2013
%
% param.scheme contains the name of the dispersion relation

switch param.scheme
  case 'single'
    [c1,c2] = fun_eval_c_single(param,k);
  case 'unconf_jetwake_sin'
    [c1,c2] = fun_eval_c_unconf_jetwake_sin(param,k);
  case 'unconf_jetwake_var'
    [c1,c2] = fun_eval_c_unconf_jetwake_var(param,k);
  case 'conf_jetwake_sin'
    [c1,c2] = fun_eval_c_conf_jetwake_sin(param,k);
  case 'conf_jetwake_var'
    [c1,c2] = fun_eval_c_conf_jetwake_var(param,k);
  otherwise
    disp('param.scheme not recognized in fun_eval_c')
end % end switch

% Evaluate w1 and w2 from c1 and c2
w1 = k.*c1; w2 = k.*c2;

end % end function