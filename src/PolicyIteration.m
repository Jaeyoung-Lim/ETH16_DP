function [ J_opt, u_opt_ind ] = PolicyIteration( P, G )
%POLICYITERATION Value iteration
%   Solve a stochastic shortest path problem by Policy Iteration.
%
%   [J_opt, u_opt_ind] = PolicyIteration(P, G) computes the optimal cost and
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

% Initialize
K = size(G,1); %Number of states
L = size(G, 2);
terminate = false;
A = zeros(K, K);
b = zeros(K,1);




% update cost based on Initial policy

u_opt_ind = ones(K,1) * L; % Stationary Policy
for i=1:K % Policy evaluation step
    b(i) = G(i, L); %Stationary policy
    A(i, :)= P(i, :, L);
end

J_opt = (eye(K) - A)\b; %Calculate inital optimal cost to go
J_temp = J_opt;

while(~terminate)
    % find new policy policy using the above cost
    
     for m=1:K
         cost_min = G(m, 1) + P(m,:,1)*J_opt;
         J_temp(m) = cost_min;
         u_opt_ind(m) = 1;
         
         for u=2:5
             c_trail = G(m, u) + P(m,:,u)*J_opt;
             
             if( c_trail <= cost_min)
                 cost_min = c_trail; %Update minimum cost to go
                 J_temp(m) = cost_min;  %For updating J_opt
                 u_opt_ind(m) = u; %Find the optimal u
             end
         end
     end
    % check condition for termination
    if max(J_opt - J_temp) < 1e-5
         terminate = true;
         u_opt_ind = u_opt_ind_temp;
    else
        J_opt = J_temp; %Upate
        u_opt_ind_temp = u_opt_ind;
    end
    
end
end