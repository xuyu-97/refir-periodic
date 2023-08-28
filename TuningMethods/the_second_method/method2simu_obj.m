function res = method2simu_obj(N, R, y, le, Ut, Vt, alpha)
% Purpose: H = Phi * K * Phi' + alpha * I
% We want to compute x'*inv(H)*x + logdet(H)
% Change part computation.

    p = size(R, 1);
    
   % Compute (R*K*R'+alpha*I_p)
   Rnew = R * (semi_matrix_multiplicaiton(Ut,Vt,R'));
   Rnew = Rnew + alpha * eye(p);
   
   % x'*inv(H)*x
   z1 = semi_vector_multiplication(Ut, Vt, y);
   z2 = y' * z1;
   z1 = R * z1;
   z1 = z1' * (Rnew \ z1);
   z = (z2-z1)/alpha;
   res_inv = (le - z)/alpha;
   
   % logdet(H)
   res_logdet = (N-p)*log(alpha) + logdet(Rnew);
   
   res = res_inv + res_logdet;
   
end

