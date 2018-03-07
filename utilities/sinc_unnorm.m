function ys = sinc_unnorm(xs)

% SINC_UNNORM  The unnormalized sinc function sin(x) / x

    ys          = sin(xs) ./ xs;
    ys(xs == 0) = 1;
end