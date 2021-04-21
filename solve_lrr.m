function [Z,E] = solve_lrr(X,S,lambda,lambda_2)
Q = orth(X');
A = X*Q;
S = Q'*S;
% A=X;
[Z,E] = lrra(X,A,S,lambda,lambda_2);
Z = Q*Z;
