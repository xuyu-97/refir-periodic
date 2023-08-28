function x = semi_vector_multiplication(Ut,Vt,x)
% K = S(U,V) = tril(U*V') + triu(V*U',1)
% Output: overwrite x by K*x
    [k, n] = size(Ut);
    v_bar = zeros(k, 1);
    u_bar = Ut * x;
    
    for i = 1:n
        v_bar = v_bar + Vt(:, i) * x(i);
        u_bar = u_bar - Ut(:, i) * x(i);
        x(i) = Ut(:, i)' * v_bar + Vt(:, i)' * u_bar;
    end
end

