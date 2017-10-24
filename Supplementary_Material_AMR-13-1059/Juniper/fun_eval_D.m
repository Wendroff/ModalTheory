function [D] = fun_eval_D(param,k,w)

% Evaluate D given w and k
% (c) Matthew Juniper 2013
%
% param.scheme contains the name of the dispersion relation

switch param.scheme
  case 'single'
    [D] = fun_eval_D_single(param,k,w);
  case 'unconf_jetwake_sin'
    [D] = fun_eval_D_unconf_jetwake_sin(param,k,w);
  case 'unconf_jetwake_var'
    [D] = fun_eval_D_unconf_jetwake_var(param,k,w);
  case 'conf_jetwake_sin'
    [D] = fun_eval_D_conf_jetwake_sin(param,k,w);
  case 'conf_jetwake_var'
    [D] = fun_eval_D_conf_jetwake_var(param,k,w);
  otherwise
    disp('param.scheme not recognized in fun_eval_D')
end % end switch

end % end function
