function x = period2whole_vector(x_p, n)
% x_p is p-dim vector (col vector), p is period. n is the dimension of x.
% Output: x = [x_p; x_p; ...; x_pr] in n-dim.

    p = length(x_p);
    x = [];
    m = floor(n / p);
    r = mod(n, p);
    for i = 1:m
        x = [x; x_p];
    end
    x = [x; x_p(1:r)];
end

