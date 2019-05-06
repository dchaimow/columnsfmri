function p=meanpower(s)
% MEANPOWER calculates the mean power.
%   p = MEANPOWER(s)

p = mean(abs(s(:)).^2);
end