%% Flatten a Matrix Row-wise
% Ex: A = [1 3; 4 2] -> [1 3 4 2]'
function A_flat = flatten(A)
    [r, c] = size(A);
    if r ~= c
        error('Input matrix must be square.');
    end
    A_flat = reshape(A.', [], 1);
end