function [ J_opt, u_opt_ind ] = LinearProgramming( P, G )
%LINEARPROGRAMMING Value iteration
%   Solve a stochastic shortest path problem by Linear Programming.
%
%   [J_opt, u_opt_ind] = LinearProgramming(P, G) computes the optimal cost
%   and the optimal control input for each state of the state space.
%
%   Input arguments:
%
%       P:
%           A (K x K x L)-matrix containing the transition probabilities
%           between all states in the state space for all control inputs.
%           The entry P(i, j, l) represents the transition probability
%           from state i to state j if control input l is applied.
%
%       G:
%           A (K x L)-matrix containing the stage costs of all states in
%           the state space for all control inputs. The entry G(i, l)
%           represents the cost if we are in state i and apply control
%           input l.
%
%   Output arguments:
%
%       J_opt:
%       	A (K x 1)-matrix containing the optimal cost-to-go for each
%       	element of the state space.
%
%       u_opt_ind:
%       	A (K x 1)-matrix containing the index of the optimal control
%       	input for each element of the state space.

% put your code here
K = size(G, 1); % Number of states
L = size(G, 2); % Number of control inputs
u_opt_ind=zeros(K, 1);

for l=1:L
    if l==1
        b= G(:, l);
        A = eye(K, K)-P(:, :, l);
        continue;
    end
    b= [b; G(:, l)];
    A = [A; eye(K, K)-P(:, :, l)];
end

A_filtered = A((b~=inf), :);
b_filtered = b(b~=inf);
f=-1*ones(1, K);

J_opt = linprog(f, A_filtered , b_filtered);

for i=1:K
    J_min = inf;
    for l=1:L
        J_temp = G(i, l)+P(i, :, l)*J_opt;
        if J_temp < J_min
            u_opt_ind(i) = l;
            J_min = J_temp;
        end
    end
end
end