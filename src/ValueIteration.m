function [ J_opt, u_opt_ind ] = ValueIteration( P, G )
%VALUEITERATION Value iteration
%   Solve a stochastic shortest path problem by Value Iteration.
%
%   [J_opt, u_opt_ind] = ValueIteration(P, G) computes the optimal cost and
%   the optimal control input for each state of the state space.
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

%% Initialize variables

assert(size(P,1) == size(P,2), ...
    'Please check dimension of Probability function');
assert(size(P,1) == size(G,1), ...
    'States between Probability and Cost function are different');
assert(size(P,3) == size(G,2), ...
    'Inputs between Probability and Cost function are different');

% variables
num_states = size(P,1);             % number of states
num_inputs = size(P,3);             % number of inputs
J_opt = zeros(num_states, 1);       % the optimal cost-to-go
u_opt_ind = zeros(num_states, 1);   % index of the optimal control input
err = 1e-5;                         % iteration error bound
%% Perform


% temporarily used variables
J_prev = J_opt;                           % old optimal cost-to-go
J_candidates = zeros(num_states, num_inputs);  % optimal cost-to-go include inputs

while(true)
    for i = 1:num_states
        for u = 1:num_inputs
            J_candidates(i,u) =  G(i, u) + P(i,:,u) * J_prev;
        end
        [J_opt(i), u_opt_ind(i)] = min(J_candidates(i,:));
    end
    if max(abs(J_opt - J_prev)) < err 
        break;
    else
        J_prev = J_opt;
    end
end
end