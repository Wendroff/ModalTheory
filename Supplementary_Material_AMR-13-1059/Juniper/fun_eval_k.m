function [k_out] = fun_eval_k(param,k_init,w,tol)

% Iterate towards k, given omega
% (c) Matthew Juniper 2013
%
% param.scheme is the name of the dispersion relation
% k_init       is the initial guess for k (complex)
% tol          is the tolerance required for the iteration

% Set the options for fsolve
options = optimset('Display','Off','TolFun',tol);

% Converge to the exact value of k using Matlab's fsolve
[k_out] = fsolve(@fun_eval_D_loc,k_init,options);

% Embedded function to call the correct dispersion relation
  function [D] = fun_eval_D_loc(k)
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
  end
% End of embedded function

end % end function
