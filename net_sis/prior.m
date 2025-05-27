%% Prior Function
%
% Currently assuming that the prior for `x` are uniformly distributed between
% 0 and 2.
function p = prior(x)
    if all(x >= 0 & x <= 2)
        %p = (1/2)^length(x);
        p = -length(x) * log(2);
    else
        p = -inf;
    end
end