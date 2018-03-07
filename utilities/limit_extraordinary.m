function [ limitpoint deriv sd ] = limit_extraordinary(u, v, valency)

% LIMIT_EXTRAORDINARY  Find limit point on a subdivision surface
%   where the required point lies in an extraordinary surface patch (i.e.
%   next to an extraordinary vertex)
%
%   See also LIMIT_ORDINARY

limitpoint = zeros(1, 6 + valency);

if u == 0 && v == 0
    [c e] = loop_limitstencil(valency);
    
    limitpoint(9) = c;
    limitpoint([5 6 8]) = e;
    
    limitpoint(10:valency + 6) = e;
    
    deriv    = zeros(3, valency + 6);
    sd       = zeros(3, valency + 6);
    
    sd(1, 9) = -valency;
    sd(1, [5 6 8]) = 1;
    sd(1, 10:valency + 6) = 1;
    
    if valency > 3
        deriv(1, [5 8 11:valency + 6 10 6]) = ...
            -cos(2 * pi * (0:valency - 1) / valency);
        deriv(2, [6 5 8 11:valency + 6 10]) = ...
            cos(2 * pi * (0:valency - 1) / valency);
        deriv(3, [6 5 8 11:valency + 6 10]) = ...
            (2 * sin(pi / valency)) * ...
            sin(2 * pi * (0:valency - 1) / valency - pi / valency);
        
        sd(2, [6 5 8 11:valency + 6 10]) = (valency / 3) * ...
            cos(4 * pi * (0:valency - 1) / valency - 2 * pi / valency);
        sd(3, [6 5 8 11:valency + 6 10]) = (valency / 3) * ...
            cos(4 * pi * (0:valency - 1) / valency - 2 * pi / valency ...
                                                   - pi / 2);
    else
        deriv(1, [5 8 6]) = -cos(2 * pi * (0:2) / 3);
        deriv(2, [6 5 8]) = cos(2 * pi * (0:2) / 3);
        deriv(3, [6 5 8]) = (2 * sin(pi / 3)) * ...
                            sin(2 * pi * (0:2) / 3 - pi / 3);
                        
        % No saddle component for valency 3: leave second derivative
        % components at 0.
    end
    
    % This is just so that the extraordinary case matches the regular one
    % (otherwise, the scaling is arbitrary)
    deriv = deriv / 3;
    
    % Convert from cup/saddle basis to basis of second derivatives in
    % parametric directions
    sd    = [ 2 -2 2 ; -1 4 -1 ; -sqrt(3) 0 sqrt(3) ] \ sd;
else
    [subdiv pick123 pick4] = loop_subdivmatrix(valency);
    extended    = subdiv;
    subdiv      = pick4 * subdiv;
    raiseto     = floor(log(u + v) / log(0.5));
    subdivsteps = (subdiv ^ raiseto);
    u           = u * 2^raiseto;
    v           = v * 2^raiseto;
    derivdenom  = (2 ^ (raiseto + 1));
    
    if v < 0.5 && u < 0.5 && u + v >= 0.5
        [lt dr sd] = limit_ordinary(1 - 2 * u, 1 - 2 * v);
        limitpoint = lt * pick123(:, :, 2) * extended * subdivsteps;
        dr         = dr * pick123(:, :, 2) * extended * subdivsteps;
        sd         = sd * pick123(:, :, 2) * extended * subdivsteps;
        dr         = -dr;
    elseif u >= 0.5
        [lt dr sd] = limit_ordinary(2 * u - 1, v * 2);
        limitpoint = lt * pick123(:, :, 3) * extended * subdivsteps;
        dr         = dr * pick123(:, :, 3) * extended * subdivsteps;
        sd         = sd * pick123(:, :, 3) * extended * subdivsteps;
    elseif v >= 0.5
        [lt dr sd] = limit_ordinary(u * 2, 2 * v - 1);
        limitpoint = lt * pick123(:, :, 1) * extended * subdivsteps;
        dr         = dr * pick123(:, :, 1) * extended * subdivsteps;
        sd         = sd * pick123(:, :, 1) * extended * subdivsteps;
    end
    
    deriv = dr * derivdenom;
    sd    = sd * (derivdenom ^ 2);
end

end

