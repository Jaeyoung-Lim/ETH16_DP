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
iter = 0;                           % number of interations
err = 0.5;                            % iteration error bound


%% Perform


% temporarily used variables
J_old = J_opt;                           % old optimal cost-to-go
J_temp = zeros(num_states, num_inputs);  % optimal cost-to-go include inputs
cost = zeros(num_states, num_inputs);    % total cost i to j

while(1)
    iter = iter + 1;
    
    for i = 1:num_states
        for u = 1:num_inputs            
            J_temp(i,u) =  G(i, u) + P(i,:,u) * J_old;
        end
        [J_opt(i), u_opt_ind(i)] = min(J_temp(i,:));
        
    end
    
    %disp(['That Value: ', num2str(sum((abs(J_opt - J_old))./(abs(J_old))))]);
    
% When you run this value interation, as you can see
% The problem is J_opt and J_old are not changing anymore
% But they are not same. So, it means J are not converging....
% I don't know why... help me XD
    
    if sum((abs(J_opt - J_old))./(abs(J_old))) < err
     %   disp('Value Iteration is done');
        break;
    elseif iter > 10000 
        disp('Value Iteration failed');
        break;
    else
        J_old = J_opt;
    end       
end
end