%% Calculate the Effective Transmission Rate (R)
% 
function R = calculate_R(A, beta_, gamma_)
    lambda_max = max(eigs(A));
    R = (beta_ * lambda_max) / gamma_;
end