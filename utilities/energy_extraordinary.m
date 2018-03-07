function energy = energy_extraordinary(valency)

% ENERGY_EXTRAORDINARY  Thin-plate energy for an extraordinary Loop surface
%   patch.

[subdiv pick123 pick4] = loop_subdivmatrix(valency);
extended = subdiv;
subdiv   = pick4 * subdiv;
    
subdiv_energy1 = extended' * pick123(:, :, 1)' * ...
                 energy_ordinary() * pick123(:, :, 1)  * extended;
subdiv_energy2 = extended' * pick123(:, :, 2)' * ...
                 energy_ordinary() * pick123(:, :, 2)  * extended;
subdiv_energy3 = extended' * pick123(:, :, 3)' * ...
                 energy_ordinary() * pick123(:, :, 3)  * extended;

subdiv_energy = 4 * (subdiv_energy1 + subdiv_energy2 + subdiv_energy3);

if valency == 3
    J = zeros(9,9);
    J(1,1) = 1.0/8.0; J(2,2) = 1.0/4.0; J(3,3) = 1.0;
    J(4,4) = 1.0/16.0; J(4,5) = 1.0; J(5,5) = 1.0/16.0;
    J(6,6) = 1.0/16.0; J(7,7) = 1.0/8.0; J(8,8) = 1.0/4.0;
    J(9,9) = 1.0/8.0;

    P = zeros(9,9);
    P(1,1) = -1.0e1; P(1,2) = 8.0/3.0; P(1,3) = 1.0/5.0;
    P(1,4) = -1.0/1.6e1; P(1,5) = 1.07e2/1.5e1;
    P(1,6) = 8.405e3/2.562e3; P(1,7) = 7.9145e4/5.49e2;
    P(1,8) = -1.20308e5/9.15e2; P(1,9) = -1.78324e5/1.647e3;
    P(2,1) = -1.0e1; P(2,2) = 1.0; P(2,3) = 1.0/5.0;
    P(2,5) = -6.0/5.0; P(2,7) = 7.7498e4/5.49e2;
    P(2,8) = -1.4581e4/3.05e2; P(2,9) = -1.76677e5/1.647e3;
    P(3,1) = -7.0; P(3,2) = -4.0/3.0; P(3,3) = 1.0/5.0;
    P(3,4) = -1.0/1.6e1; P(3,5) = 1.22e2/1.5e1;
    P(3,6) = 3.281e3/2.562e3; P(3,7) = 5.40839e5/5.49e3;
    P(3,8) = 6.1984e4/9.15e2; P(3,9) = -2.66453e5/3.294e3;
    P(4,2) = 1.0; P(4,3) = 1.0/5.0; P(4,5) = -6.0/5.0; P(4,7) = 3.0;
    P(4,8) = -1.5496e4/3.05e2; P(4,9) = -1.0; P(5,2) = 2.0/3.0;
    P(5,3) = 1.0/5.0; P(5,5) = 2.0/1.5e1; P(5,8) = -3.0077e4/9.15e2;
    P(6,2) = -1.0/3.0; P(6,3) = 1.0/5.0; P(6,5) = 2.0/1.5e1;
    P(6,8) = 1.5496e4/9.15e2; P(7,1) = 3.0; P(7,2) = -2.0;
    P(7,3) = 1.0/5.0; P(7,5) = -6.0/5.0; P(7,7) = -7.8047e4/1.83e3;
    P(7,8) = 3.0077e4/3.05e2; P(7,9) = 2.8967e4/1.098e3;
    P(8,2) = -1.0/3.0; P(8,3) = 1.0/5.0; P(8,5) = 2.0/1.5e1;
    P(8,8) = 1.4581e4/9.15e2; P(9,3) = 1.0/5.0; P(9,5) = -1.0/5.0;

    infsumcoeffs = 4 .* repmat(diag(J), 1, 9) .* ...
                        repmat(diag(J)', 9, 1);
    infsumcoeffs(infsumcoeffs >= 1 - 1e-15) = 0;
    infsumcoeffs = infsumcoeffs ./ (1 - infsumcoeffs);

    Z = (P' * subdiv_energy * P);
    infsum = Z .* infsumcoeffs;

    % 1024 / 3969 ==
    %   \sum_{i = 1}^\infty 4^i \frac{1}{16}^i
    %                         i \frac{1}{16}^{i - 1}
    infsum(:, 5) = infsum(:, 5) + 1024 * Z(:, 4) / 3969;
    infsum(5, :) = infsum(5, :) + 1024 * Z(4, :) / 3969;

    % 1064960 / 250047 ==
    %   \sum_{i = 1}^\infty 4^i \frac{1}{16}^{i - 1}
    %                       i^2 \frac{1}{16}^{i - 1}
    infsum(5, 5) = infsum(5, 5) + 1064960 * Z(4, 4) / 250047;

    energy = infsum / P;
    energy = (energy' / P)';
    energy = subdiv_energy + energy;
else
    [outofeigen evals] = eig(subdiv);
    outofeigen = real(outofeigen);
    evals = real(evals);

    infsumcoeffs = 4 .* repmat(diag(evals), 1, 6 + valency) .* ...
                        repmat(diag(evals)', 6 + valency, 1);

    % See "Efficient, Fair Interpolation using Catmull-Clark Surfaces",
    % Halstead et al.
    % To calculate this infinite sum we just set the divergent terms to
    % zero in this matrix.
    infsumcoeffs(infsumcoeffs >= 1 - 1e-15) = 0;
    infsumcoeffs = infsumcoeffs ./ (1 - infsumcoeffs);

    infsum = (outofeigen' * subdiv_energy * outofeigen) .* infsumcoeffs;

    energy = infsum / outofeigen;
    energy = (energy' / outofeigen)';
    energy = subdiv_energy + energy;
end

end