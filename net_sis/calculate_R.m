%% Calculate the Effective Transmission Rate (R)
% 
function R = calculate_R(A, gamma_)
    lambda_max = max(eigs(A));
    R = lambda_max / gamma_;
end