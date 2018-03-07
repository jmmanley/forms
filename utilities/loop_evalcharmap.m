function [ xy deriv sd ] = loop_evalcharmap(valency, u, v)

% LOOP_EVALCHARMAP  Evaluate the characteristic map of the Loop subdivision
%   scheme.

nc = loop_natconf(valency);

if valency ~= 6
    [ xy deriv sd ] = limit_extraordinary(u, v, valency);
else
    [ xy deriv sd ] = limit_ordinary(u, v);
end

xy    = xy    * nc;
deriv = deriv * nc;
sd    = sd    * nc;

end

