function output = script_click(param)

% With param already in memory,
% click on a point in the figure and find the nearest saddle point.
% Determine whether it is Convectively or Absolutely unstable.
% (c) Matthew Juniper 2013

% Click on a point on the figure
click = ginput(1);
k_init = complex(click(1),click(2));
tol = 10^-6;

% Converge to the nearest point, where dwdk = 0
[k,w] = fun_eval_dwdk0(param,k_init,tol);

% Plot the point accurately

% Determine whether absolutely or convectively unstable
if imag(w) > 0
  disp('Absolutely unstable')
  plot(k,'ro','MarkerFaceColor','r')
elseif imag(w) < 0
  disp('Convectively unstable')
  plot(k,'go','MarkerFaceColor','g')
else
  disp('Marginally unstable')
  plot(k,'yo','MarkerFaceColor','y')
end

output = {};

end
