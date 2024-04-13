# GREEN
clc;
k_1 = 0.2562;
k_0 = 0.1984;

function [G] = Green1(x,xp,k_0,k_1)

%% INPUTS
% x = observaion point
% xp = source point
% k_0 = wavevector in x<0 domain [1/length]
% k_1 = wavevector in x>0 domain

%% OUTPUT
% G = Greens function value
% if k_0 = k_1, G = formula for hoogeneous medium
% if k_0 ~= k_1, G = formula for junction between two homogeneous media

tol = 1e-12 ;

junction = true;

k_diff = (k_1 - k_0);

if(abs(k_diff)<tol)
   junction = false;
end

% Homogeneous medium

if ~junction
  G = exp(1i*k_1*abs(x-xp))/(2*1i*k_1);
  return;
end

% Junction between two homogeneous media

k_sum = k_1 + k_0 ;

reflex = k_diff/k_sum ;

if (x>0 && xp>0)
  G = exp(1i*k_1*abs(x-xp))/(2*1i*k_1) + reflex * exp(1i*k_1*(x+xp))/(2*1i*k_1);
  return;
end

if (x<0 && xp>0)
  G = exp(-1i*k_0*x + 1i*k_1*xp) / (1i*k_sum);
  return;
end

if ( x>0 && xp<0)
  error(['Green1: x>0 and xp<0 not implemented because ' ...
           'should not be needed.'])

end

end % ebd of function Green1









