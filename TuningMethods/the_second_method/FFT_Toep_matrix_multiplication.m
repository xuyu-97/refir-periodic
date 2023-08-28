function R = FFT_Toep_matrix_multiplication(T,X)
% T: toeplitz matrix with m*n dimension
% X: matrix with n * k dimension where n, k << m
% Calculate using fast fourier transformation, time is k * (m+(n-1)) *
% log(m+(n-1))
    [m, n] = size(T);
    k = size(X, 2);
    
    if n ~= size(X, 1)
        error("The dimension is wrong!") 
    end
    
    v = fft([T(:,1); T(1,end:-1:2)']);
    Y = fft([X ; zeros(m-1, k)]);
    R = ifft(v.*Y);
    R = R(1:m, 1:k);
    return
end

