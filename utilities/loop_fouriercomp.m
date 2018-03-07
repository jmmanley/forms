function [ zf nzf1 nzf2 ] = loop_fouriercomp(valency, extra_rings)

% LOOP_FOURIERCOMP  Returns the Loop subdivision matrices for the zero,
%   unit and double frequency components of a Discrete Fourier Transform

if nargin < 2, extra_rings = 0;
else assert(extra_rings >= 0);
end

nzf_size = (extra_rings + 2) * (extra_rings + 3) / 2;
zf_size = nzf_size + 1;

zf   = zeros(zf_size);
nzf1 = zeros(nzf_size);
nzf2 = zeros(nzf_size);

omega = exp(2 * pi * sqrt(-1) / valency);
omesq = omega ^ 2;

alpha = (0.375 + cos(2 * pi / valency) / 4) ^ 2 + 0.375;

zf(1:4,1:4) = [ alpha 1 - alpha 0 0 ;
                6     10        0 0 ;
                1     12        1 2 ;
                2     12        0 2 ];

nzf1(1:3,1:3) = [  6 + 2 * omega + 2 * conj(omega)  0   0               ;
                  10 +     omega +     conj(omega)  1   1 + conj(omega) ;
                   6 + 6 * omega                    0   2               ];

nzf2(1:3,1:3) = [  6 + 2 * omesq + 2 * conj(omesq)  0   0               ;
                  10 +     omesq +     conj(omesq)  1   1 + conj(omesq) ;
                   6 + 6 * omesq                    0   2               ];
               
for pt = 4:size(nzf1, 1)
    i = floor(0.5 + sqrt(2 * pt - (1/4)));
    j = pt - (1 + i * (i - 1) / 2);
    
    if mod(i, 2) == 0 && mod(j, 2) == 0
        % vertex point
        update(pt, i / 2, j / 2, 10);
        update(pt, (i / 2) - 1, j / 2, 1);
        update(pt, (i / 2) - 1, (j / 2) - 1, 1);
        update(pt, (i / 2), (j / 2) - 1, 1);
        update(pt, (i / 2), (j / 2) + 1, 1);
        update(pt, (i / 2) + 1, j / 2, 1);
        update(pt, (i / 2) + 1, (j / 2) + 1, 1);
    else
        % edge point
        if mod(i, 2) == 1 && mod(j, 2) == 1
            update(pt, floor(i / 2), floor(j / 2), 6);
            update(pt,  ceil(i / 2),  ceil(j / 2), 6);
            update(pt, floor(i / 2),  ceil(j / 2), 2);
            update(pt,  ceil(i / 2), floor(j / 2), 2);
        elseif mod(i, 2) == 1
            update(pt, floor(i / 2), j / 2, 6);
            update(pt,  ceil(i / 2), j / 2, 6);
            update(pt, floor(i / 2), (j / 2) - 1, 2);
            update(pt,  ceil(i / 2), (j / 2) + 1, 2);
        else % mod(j, 2) == 1
            update(pt, i / 2, floor(j / 2), 6);
            update(pt, i / 2, ceil(j / 2), 6);
            update(pt, (i / 2) - 1, floor(j / 2), 2);
            update(pt, (i / 2) + 1,  ceil(j / 2), 2);
        end
    end
end
               
% Normalize
denoms = sum(zf, 2);
zf = zf ./ repmat(denoms, 1, size(zf, 2));
denoms = denoms(2:end);
nzf1 = nzf1 ./ repmat(denoms, 1, size(nzf1, 2));
nzf2 = nzf2 ./ repmat(denoms, 1, size(nzf2, 2));

    function update(to, fromi, fromj, w)
        if fromj >= fromi
            from = 1 + fromi * (fromi - 1) / 2;
            nzf1(to, from) = nzf1(to, from) + omega * w;
            nzf2(to, from) = nzf2(to, from) + omesq * w;
        elseif fromj < 0
            fromi = fromi + 1;
            from = fromi + fromi * (fromi - 1) / 2;
            nzf1(to, from) = nzf1(to, from) + conj(omega) * w;
            nzf2(to, from) = nzf2(to, from) + conj(omesq) * w;
        else
            from = 1 + fromj + fromi * (fromi - 1) / 2;
            nzf1(to, from) = nzf1(to, from) + w;
            nzf2(to, from) = nzf2(to, from) + w;
        end
        zf(1 + to, 1 + from) = zf(1 + to, 1 + from) + w;
    end
end