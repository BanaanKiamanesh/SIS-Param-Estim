% Function to Calculate the Jacobian of ODE
function J = Jacobian(X)
    % State Unpack
    % X = [x1 x2 a11 a12 a21 a22 Gamma]'
    x1 = X(1);  
    x2 = X(2);
    
    a11 = X(3);  
    a12 = X(4);
    a21 = X(5);  
    a22 = X(6);
    
    Gamma  = X(7);

    % Partials for Node 1
    df1dx1   =  a11 - 2*a11*x1 - a12*x2 - Gamma;
    df1dx2   = (1 - x1)*a12;
    df1da11  = (1 - x1)*x1;
    df1da12  = (1 - x1)*x2;
    df1dg    = -x1;

    % Partials for Node 2
    df2dx1   = (1 - x2)*a21;
    df2dx2   =  a22 - a21*x1 - 2*a22*x2 - Gamma;
    df2da21  = (1 - x2)*x1;
    df2da22  = (1 - x2)*x2;
    df2dg    = -x2;

    % Assemble 7×7 Jacobian
    J        = zeros(7);
    J(1,1)   = df1dx1;
    J(1,2)   = df1dx2;
    J(1,3)   = df1da11;
    J(1,4)   = df1da12;
    J(1,7)   = df1dg;

    J(2,1)   = df2dx1;
    J(2,2)   = df2dx2;
    J(2,5)   = df2da21;
    J(2,6)   = df2da22;
    J(2,7)   = df2dg;
end