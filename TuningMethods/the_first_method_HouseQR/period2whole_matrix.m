function A = period2whole_matrix(A_p, N, n)
% A_p is p*p matrix, p is period. A has size N*n.
% Output: extend A_p to A.
    
    % Extend columns.
    p = size(A_p, 1);
    A = [];
    m = floor(n / p);
    r = mod(n, p);
    for i=1:m
        A = [A, A_p];
    end
    A = [A, A_p(:,1:r)];
    
    % Extend rows.
    m = floor(N / p);
    r = mod(N, p);
    A_p = A;
    for i=1:(m-1)
        A = [A; A_p];
    end
    A = [A; A_p(1:r,:)];
    
end

