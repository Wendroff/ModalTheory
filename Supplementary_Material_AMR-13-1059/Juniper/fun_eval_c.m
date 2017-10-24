function [c1,c2] = fun_eval_c(param,k)

% Evaluate c1 and c2, given k
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

end % end function
