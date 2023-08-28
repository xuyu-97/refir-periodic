function X = semi_matrix_multiplicaiton(Ut,Vt,X)
% K = S(U,V) = tril(U*V') + triu(V*U',1)
% Output: overwrite x by K*X
    k = size(X,2);
    for i=1:k
        X(:,i) = semi_vector_multiplication(Ut,Vt, X(:,i));
    end

end

