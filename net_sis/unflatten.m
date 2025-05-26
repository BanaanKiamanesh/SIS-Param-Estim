function A = unflatten(A_flat)
    A_flat = A_flat(:);  % Ensure column vector
    n = sqrt(numel(A_flat));
    if mod(n,1) ~= 0
        error('Input length must be a perfect square.');
    end

    A = reshape(A_flat, n, n);
end