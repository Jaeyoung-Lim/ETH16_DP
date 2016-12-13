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

K = size(G,1);

% initial guess on policy: taking photo everywhere
policy = ones(K,1) * 5;
cost = zeros(K,1);

terminate = false;
itr_cnt = 0;

while(!terminate)

    printf("iteration: %d\n", itr_cnt)

    % update cost based on current policy
    A = eye(K);
    b = zeros(K,1);
    for m=1:K
        u = policy(m);
        if(G(m,u)==inf)
            b(m) = 1e10;
        else
            b(m) = G(m, u);
        end
        for n=1:K
            if( (P(m,n,u)>0) )
                A(m,n) = A(m,n) - P(m,n,u);
            end
        end
    end
    cost_tmp = inv(A)*b;

    % find new policy policy using the above cost
    policy_tmp = zeros(K,1);
    for m=1:K
        c = G(m, 1) + P(m,:,1)*cost_tmp;
        policy_tmp(m) = 1;
        for u=2:5
            c_trail = G(m, u) + P(m,:,u)*cost_tmp;
            if( c_trail < c)
                c = c_trail;
                policy_tmp(m) = u;
            end
        end
    end

    % check condition for termination
    itr_cnt = itr_cnt + 1;
    if( isequal(cost, cost_tmp) || (itr_cnt>1000) )
        terminate = true;
    end

    % update policy
    policy = policy_tmp;
    cost = cost_tmp;
end

J_opt = cost;
u_opt_ind = policy;

end

