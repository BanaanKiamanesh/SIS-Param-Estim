%% Prior Function
% Calculate the Log Probability for a particular set of components, 
% assuming that they are all uniformily distributed between 0 and 3.
function p = prior(x)
    if all(x >= 0 & x <= 3)
        p = -length(x) * log(3);
    else
        p = -inf;
    end
end