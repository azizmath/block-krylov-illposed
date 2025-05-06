function [A,rhs,sol] = fulltomo(n)

rand('seed',0)
[A rhs sol] = tomo(round(sqrt(n)));
A = full(A);

