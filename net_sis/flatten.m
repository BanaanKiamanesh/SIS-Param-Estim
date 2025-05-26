function A_flat = flatten(A)
    [r, c] = size(A);
    if r ~= c
        error('Input matrix must be square.');
    end
    A_flat = A(:);
end