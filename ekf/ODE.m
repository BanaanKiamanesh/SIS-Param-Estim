% Function to Simulate the System ODE
function dX = ODE(X)
    % State Unpack
    % X = [x1 x2 a11 a12 a21 a22 Gamma]'
    x1  = X(1);  
    x2  = X(2);
    
    a11 = X(3);  
    a12 = X(4);
    a21 = X(5);  
    a22 = X(6);
    
    Gamma  = X(7);

    % State Derivative Vector
    dx1 = (1 - x1)*(a11*x1 + a12*x2) - Gamma*x1;
    dx2 = (1 - x2)*(a21*x1 + a22*x2) - Gamma*x2;

    dX   = [dx1; dx2; zeros(5,1)];
end