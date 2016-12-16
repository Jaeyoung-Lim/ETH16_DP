function Check(  stateSpace, controlSpace, ...
        map, gate, mansion, cameras, tol_P, tol_G, tol_VI, tol_PI, tol_LP )
    
    % P
    P_DH = ComputeTransitionProbabilitiesDH(...
        stateSpace, ...
        controlSpace, ...
        map, ...
        gate, ...
        mansion, ...
        cameras );
    P = ComputeTransitionProbabilities(...
        stateSpace, ...
        controlSpace, ...
        map, ...
        gate, ...
        mansion, ...
        cameras );
    
    assert(isequal(size(P_DH), size(P)), 'P matrix dim error');
    Pdiff = abs(P_DH - P);
    disp(sprintf('\tP: %s', DispResult(Pdiff, tol_P)))

    % G
    G_DH = ComputeStageCostsDH(...
        stateSpace, ...
        controlSpace, ...
        map, ...
        gate, ...
        mansion, ...
        cameras );
    G = ComputeStageCosts(...
        stateSpace, ...
        controlSpace, ...
        map, ...
        gate, ...
        mansion, ...
        cameras );
    
    assert(isequal(size(G_DH), size(G)), 'G matrix dim error');
    Gdiff = abs(G_DH - G);
    Gdiff(isnan(Gdiff)) = 0;
    disp(sprintf('\tG: %s', DispResult(Gdiff, tol_G)));
    
    % VI
    [ J_opt_vi_DH, u_opt_ind_vi_DH ] = ValueIterationDH( P_DH, G_DH );
    [ J_opt_vi, u_opt_ind_vi ] = ValueIteration( P, G );

    assert(isequal(size(J_opt_vi_DH), size(J_opt_vi)), 'J_opt_vi matrix dim error');
    assert(isequal(size(u_opt_ind_vi_DH), size(u_opt_ind_vi)), 'u_opt_ind_vi matrix dim error');

    J_opt_vi_diff = abs(J_opt_vi_DH - J_opt_vi);
    u_opt_ind_vi_diff = abs(u_opt_ind_vi_DH - u_opt_ind_vi);
    disp(sprintf('\tVI:'))
    disp(sprintf('\t\tJ_opt : %s\n\t\tu_opt : %s', ...
        DispResultJ(J_opt_vi_diff, tol_VI), ...
        DispResultU(u_opt_ind_vi_diff)));
    
    % PI
    [ J_opt_pi_DH, u_opt_ind_pi_DH ] = PolicyIterationDH( P_DH, G_DH );
    [ J_opt_pi, u_opt_ind_pi ] = PolicyIteration( P, G );
    
    assert(isequal(size(J_opt_pi_DH), size(J_opt_pi)), 'J_opt_pi matrix dim error');
    assert(isequal(size(u_opt_ind_pi_DH), size(u_opt_ind_pi)), 'u_opt_ind_pi matrix dim error');
    
    J_opt_pi_diff = abs(J_opt_pi_DH - J_opt_pi);
    u_opt_ind_pi_diff = abs(u_opt_ind_pi_DH - u_opt_ind_pi);
    disp(sprintf('\tPI:'))
    disp(sprintf('\t\tJ_opt : %s\n\t\tu_opt : %s', ...
        DispResultJ(J_opt_pi_diff, tol_PI), ...
        DispResultU(u_opt_ind_pi_diff)));
    
    % LP
    [ J_opt_lp_DH, u_opt_ind_lp_DH ] = LinearProgrammingDH( P_DH, G_DH );
    [ J_opt_lp, u_opt_ind_lp ] = LinearProgramming( P, G );
    
    assert(isequal(size(J_opt_lp_DH), size(J_opt_lp)), 'J_opt_lp matrix dim error');
    assert(isequal(size(u_opt_ind_lp_DH), size(u_opt_ind_lp)), 'u_opt_ind_lp matrix dim error');
    
    J_opt_lp_diff = abs(J_opt_lp_DH - J_opt_lp);
    u_opt_ind_lp_diff = abs(u_opt_ind_lp_DH - u_opt_ind_lp);
    disp(sprintf('\tLP:'))
    disp(sprintf('\t\tJ_opt : %s\n\t\tu_opt : %s', ...
        DispResultJ(J_opt_lp_diff, tol_LP), ...
        DispResultU(u_opt_ind_lp_diff)));
end

% subfunction for disp
function str = DispResult(diff, tol)
    yn = all(diff(:) < tol);
    score = sum(diff(:) < tol);
    
    if yn == true
        str = sprintf('score = %8d/%8d\t\t(O)', score, numel(diff));
    else
        str = sprintf('score = %8d/%8d\t\t(X)', score, numel(diff));
    end 
end

function str = DispResultJ(diff, tol)
    yn = all(diff(:) < tol);
    score = sum(diff(:) < tol);
    
    if yn == true
        str = sprintf('score = %6d/%6d\t(O)', score, numel(diff));
    else
        str = sprintf('score = %6d/%6d\t(X)', score, numel(diff));
    end 
end

function str = DispResultU(diff)
    yn = all(diff(:) == 0);
    score = sum(diff(:) == 0);
    
    if yn == true
        str = sprintf('score = %6d/%6d\t(O)', score, numel(diff));
    else
        str = sprintf('score = %6d/%6d\t(X)', score, numel(diff));
    end 
end


% to make my codes secured....
% P and G
function P = ComputeTransitionProbabilitiesDH( stateSpace, controlSpace, map, gate, mansion, cameras )
    %COMPUTETRANSITIONPROBABILITIES Compute transition probabilities.
    % 	Compute the transition probabilities between all states in the state
    %   space for all control inputs.
    %
    %   P = ComputeTransitionProbabilities(stateSpace, controlSpace,
    %   map, gate, mansion, cameras) computes the transition probabilities
    %   between all states in the state space for all control inputs.
    %
    %   Input arguments:
    %
    %       stateSpace:
    %           A (K x 2)-matrix, where the i-th row represents the i-th
    %           element of the state space.
    %
    %       controlSpace:
    %           A (L x 1)-matrix, where the l-th row represents the l-th
    %           element of the control space.
    %
    %       map:
    %           A (M x N)-matrix describing the terrain of the estate map.
    %           Positive values indicate cells that are inaccessible (e.g.
    %           trees, bushes or the mansion) and negative values indicate
    %           ponds or pools.
    %
    %   	gate:
    %          	A (1 x 2)-matrix describing the position of the gate.
    %
    %    	mansion:
    %          	A (F x 2)-matrix indicating the position of the cells of the
    %           mansion.
    %
    %    	cameras:
    %          	A (H x 3)-matrix indicating the positions and quality of the
    %           cameras.
    %
    %   Output arguments:
    %
    %       P:
    %           A (K x K x L)-matrix containing the transition probabilities
    %           between all states in the state space for all control inputs.
    %           The entry P(i, j, l) represents the transition probability
    %           from state i to state j if control input l is applied.
    %
    %   written by Dongho
    
    K = size(stateSpace, 1);    % num of states
    L = size(controlSpace, 1);  % num of controls
    
    % P matrix init
    P = zeros(K, K, L);
    
    camMat = SecurityCameraMatrix(map, cameras);
    papMat = PaparazziCameraMatrix(map, mansion);
    
    % global vars
    global p_c pool_num_time_steps
    
    for i=1:L 
        % iteration for control input
        controlInput = controlSpace(i);
        
        for j=1:K
            % iteration for states
            
            % currentState
            currentState = stateSpace(j,:);
            
            % next state by control input
            nextState = NextStateByControlInput(map, currentState, controlInput);
            x = nextState(1);
            y = nextState(2);
            
            if controlInput == 'p'
                % taking picture case
                
                % caculte prob to take successful picture
                probToSucceed = max(1 - prod(1 - papMat(y, x, :)), p_c);
                probToFailed = 1 - probToSucceed;
                
                % cacluate prob to be caught by security camera
                probToBeCaught = 1 - prod(1 - camMat(y, x, :));
                probNotToBeCaught = 1 - probToBeCaught;
                
                % update P matrix
                % caught! (and failed to take picture) go to gate
                P(j, FindStateIndex(gate, stateSpace), i) = ...
                    P(j, FindStateIndex(gate, stateSpace), i) + probToFailed * probToBeCaught;            
                
                % not caught! (and failed to take picture) stay here
                P(j, FindStateIndex(nextState, stateSpace), i) = ...
                    P(j, FindStateIndex(nextState, stateSpace), i) + probToFailed * probNotToBeCaught;
            else
                 % otherwise
                if map(y, x) < 0
                    % pond/pool case 
                    % calculate prob to be caught by security camera
                    
                    % else, time step will be pool_num_time_steps
                    probNotToBeCaught = prod(1 - camMat(y, x, :)) ^ pool_num_time_steps;
                    probToBeCaught = 1 - probNotToBeCaught;
                    
                else
                    % calculate prob to be caught by security camera
                    probNotToBeCaught = prod(1 - camMat(y, x, :));
                    probToBeCaught = 1 - probNotToBeCaught;
                end
                                                
                % update P matrix
                % caught! go to gate
                P(j, FindStateIndex(gate, stateSpace), i) = ...
                    P(j, FindStateIndex(gate, stateSpace), i) + probToBeCaught;             
                
                % not caught! go to nextState
                P(j, FindStateIndex(nextState, stateSpace), i) = ...
                    P(j, FindStateIndex(nextState, stateSpace), i) + probNotToBeCaught;  
            end
        end
    end
end

function G = ComputeStageCostsDH( stateSpace, controlSpace, map, gate, mansion, cameras )
    %COMPUTESTAGECOSTS Compute stage costs.
    % 	Compute the stage costs for all states in the state space for all
    %   control inputs.
    %
    %   G = ComputeStageCosts(stateSpace, controlSpace, map, gate, mansion,
    %   cameras) computes the stage costs for all states in the state space
    %   for all control inputs.
    %
    %   Input arguments:
    %
    %       stateSpace:
    %           A (K x 2)-matrix, where the i-th row represents the i-th
    %           element of the state space.
    %
    %       controlSpace:
    %           A (L x 1)-matrix, where the l-th row represents the l-th
    %           element of the control space.
    %
    %       map:
    %           A (M x N)-matrix describing the terrain of the estate map.
    %           Positive values indicate cells that are inaccessible (e.g.
    %           trees, bushes or the mansion) and negative values indicate
    %           ponds or pools.
    %
    %   	gate:
    %          	A (1 x 2)-matrix describing the position of the gate.
    %
    %    	mansion:
    %          	A (F x 2)-matrix indicating the position of the cells of the
    %           mansion.
    %
    %    	cameras:
    %          	A (H x 3)-matrix indicating the positions and quality of the
    %           cameras.
    %
    %   Output arguments:
    %
    %       G:
    %           A (K x L)-matrix containing the stage costs of all states in
    %           the state space for all control inputs. The entry G(i, l)
    %           represents the cost if we are in state i and apply control
    %           input l.
    %
    %   written by Dongho
    
    K = size(stateSpace, 1);    % num of states
    L = size(controlSpace, 1);  % num of controls
    
    % G matrix init
    G = Inf(K, L);

    camMat = SecurityCameraMatrix(map, cameras);
    papMat = PaparazziCameraMatrix(map, mansion);
    
    % global vars
    global p_c pool_num_time_steps detected_additional_time_steps
    
    for i=1:L
        % iteration for control input
        controlInput = controlSpace(i);
        
        for j=1:K
            % iteration for states
            
            % currentState
            currentState = stateSpace(j,:);
            
            % next state by control input
            nextState = NextStateByControlInput(map, currentState, controlInput);
            x = nextState(1);
            y = nextState(2);
            
            if controlInput == 'p'
                % taking picture case
                
                % caculte prob to take successful picture
                probToSucceed = max(1 - prod(1 - papMat(y, x, :)), p_c);
                probToFailed = 1 - probToSucceed;
                
                % cacluate prob to be caught by security camera
                probToBeCaught = 1 - prod(1 - camMat(y, x, :));
                probNotToBeCaught = 1 - probToBeCaught;
          
                % costs
                costToSucceed = 1;
                costToBeCaught = 1 + detected_additional_time_steps;
                costNotToBeCaught = 1;
                
                % update G matrix
                G(j, i) = probToSucceed * costToSucceed ...
                    + (probToFailed * probToBeCaught) * costToBeCaught ...
                    + (probToFailed * probNotToBeCaught) * costNotToBeCaught;
            else
                % otherwise
                if map(y, x) < 0
                    % pond/pool case 
                    % calculate prob to be caught by security camera
                    
                    % else, time step will be pool_num_time_steps
                    probNotToBeCaught = prod(1 - camMat(y, x, :)) ^ pool_num_time_steps;
                    probToBeCaught = 1 - probNotToBeCaught;
                    
                    costToMove = pool_num_time_steps;
                    costToBeCaught = pool_num_time_steps + detected_additional_time_steps;
                    
                else
                    % calculate prob to be caught by security camera
                    probToBeCaught = 1 - prod(1 - camMat(y, x, :));
                    probNotToBeCaught = 1 - probToBeCaught;
                                        
                    costToMove = 1;
                    costToBeCaught = 1 + detected_additional_time_steps;
                end
                
                if nextState == currentState
                    % inf when next state is inaccessible
                    costToMove = inf;
                    costToBeCaught = 1 + detected_additional_time_steps;
                end
                
                % update G matrix
                G(j, i) = probToBeCaught * costToBeCaught ...
                    + probNotToBeCaught * costToMove;
            end
        end
    end
 end

function S = SecurityCameraMatrix ( map, cameras )
    %   compute S matrix (size: M x N x H) which contains probability to catch
    %   by security camera for each cells.
    %
    %   written by Dongho 
    
    M = size(map, 1);
    N = size(map, 2);
    H = size(cameras, 1);
    
    % S matrix init
    S = zeros(M, N, H);
    
    % security camera obstacle map 
    % in this case, same with map
    map_obs = map;
    
    for i=1:H
        % iterate for each camera.
        cam = cameras(i,:);
        
        x = cam(1); % x coord
        y = cam(2); % y coord
        r = cam(3); % resolution
        
        for j=1:x-1
            % iterate for left-side of camera
            col = x-j;
            row = y;
            
            % obstacle! break iteration
            if map_obs(row, col) > 0; break; end 
                         
            S(row,col,i) = abs(r/j);
        end
        
        for j=1:N-x
            % iterate for right-side of camera
            col = x+j;
            row = y;
            
            % obstacle! break iteration
            if map_obs(row, col) > 0; break; end 
            
            S(row,col,i) = abs(r/j);
        end
        
        for j=1:y-1
            % iterate for down-side of camera
            col = x;
            row = y-j;
            
            % obstacle! break iteration
            if map_obs(row, col) > 0; break; end 
            
            S(row,col,i) = abs(r/j);
        end
        
        for j=1:M-y 
            % iterate for up-side of camera
            col = x;
            row = y+j;
            
            % obstacle! break iteration
            if map_obs(row, col) > 0; break; end 
            
            S(row,col,i) = abs(r/j);
        end 
    end
end

function C = PaparazziCameraMatrix ( map, mansion )
    %   compute C matrix (size: M x N) which contains probability to take 
    %   successful photo by paparazzi
    %
    %   written by Dongho 
    
    M = size(map, 1);
    N = size(map, 2);
    H = size(mansion, 1);
    
    global gamma_p;
    
    % C matrix init
    C = zeros(M, N, H);
    
    % paparazzi camera obstacle map
    % in this case, same with map
    map_obs = map;
    
    for i=1:H
        % iterate for each mansion block.
        man = mansion(i,:);
        
        x = man(1); % x coord
        y = man(2); % y coord
                
        for j=1:x-1
            % iterate for left-side of mansion block
            col = x-j;
            row = y;
            
            % obstacle! break iteration
            if map_obs(row, col) > 0; break; end 
                         
            C(row,col,i) = abs(gamma_p/j);
        end
        
        for j=1:N-x
            % iterate for right-side of mansion block
            col = x+j;
            row = y;
            
            % obstacle! break iteration
            if map_obs(row, col) > 0; break; end 
            
            C(row,col,i) = abs(gamma_p/j);
        end
        
        for j=1:y-1
            % iterate for down-side of mansion block
            col = x;
            row = y-j;
            
            % obstacle! break iteration
            if map_obs(row, col) > 0; break; end 
            
            C(row,col,i) = abs(gamma_p/j);
        end
        
        for j=1:M-y 
            % iterate for up-side of mansion block
            col = x;
            row = y+j;
            
            % obstacle! break iteration
            if map_obs(row, col) > 0; break; end 
            
            C(row,col,i) = abs(gamma_p/j);
        end 
    end    
end

function nextState = NextStateByControlInput( map, currentState, controlInput )
    %   return nextstate considering map boundary, obstacles (inaccessible)
    %   but not considering security camera
    %
    %   written by Dongho
    
    x = currentState(1);
    y = currentState(2);
    
    M = size(map, 1);
    N = size(map, 2);
    
    % init
    nextState = currentState;  
    
    % transition by control input
    switch controlInput
        case 'n'
            nextState = currentState + [0 1];
        case 'w'
            nextState = currentState + [-1 0];
        case 's'
            nextState = currentState + [0 -1];
        case 'e'
            nextState = currentState + [1 0];
        case 'p'
            nextState = currentState;
        otherwise
            nextState = currentState;
    end
    
    if ~(N >= nextState(1) && nextState(1) > 0 ...
            && M >= nextState(2) && nextState(2) > 0)
        % boundary check!
        nextState = currentState;
    elseif map(nextState(2), nextState(1)) > 0
        % obstacle check!
        nextState = currentState;
    end
end

function index = FindStateIndex( state, stateSpace )
    % written by Dongho 
    index = find(ismember(stateSpace, state, 'rows') == true);
end

% VI
function [ J_opt, u_opt_ind ] = ValueIterationDH( P, G )
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
    
    %   written by Dongho
    
    K = size(P, 1); % num of states
        
    % J_opt and u_opt_ind init
    J_opt = zeros(K, 1);        % if there's better IC, change it!
    u_opt_ind = ones(K, 1);
    
    % do first iteration
    [J_opt_next, u_opt_ind_next] = VIteration( J_opt, P, G );
    
    while ~(max(abs(J_opt - J_opt_next)) < 1e-5) 
        % iterate until converge
        J_opt = J_opt_next;
        u_opt_ind = u_opt_ind_next;
        
        [J_opt_next, u_opt_ind_next] = VIteration( J_opt, P, G );
    end
end

function [ J_opt_next, u_opt_ind_next] = VIteration( J_opt, P, G )
    % written by Dongho
    
    K = size(P, 1); % num of states
    L = size(P, 3); % num of controls
        
    % save calculated cost
    cost = zeros(K, L);
    
    % calculate cost
    for l=1:L
        cost(:, l) = G(:, l) + P(:, :, l) * J_opt;
    end
    
    % find minimum and update J_opt
    [J_opt_next, u_opt_ind_next] = min(cost, [], 2);
end

% PI
function [ J_opt, u_opt_ind ] = PolicyIterationDH( P, G )
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
    
    %   written by Dongho
    
    K = size(P, 1); % num of states
    L = size(P, 3); % num of controls
    
    % J_opt and u_opt_ind init
    J_opt = zeros(K, 1);       
    u_opt_ind = L * ones(K, 1);     % in PI, IC is important!
    
    % solve linear system equation first
    [G_u, P_u] = FindForU(P, G, u_opt_ind);
    J_opt = (eye(K) - P_u)\G_u;            
    
    % start iteration
    [J_opt_next, u_opt_ind_next] = PIteration( J_opt, P, G );
    
    while ~(max(abs(J_opt - J_opt_next)) < 1e-5)
        % iterate until converge
        J_opt = J_opt_next;
        u_opt_ind = u_opt_ind_next;
        
        [J_opt_next, u_opt_ind_next] = PIteration( J_opt, P, G );
    end
end

function [ J_opt_next, u_opt_ind_next] = PIteration( J_opt, P, G )
    % written by Dongho
    
    K = size(P, 1); % num of states
    L = size(P, 3); % num of controls
        
    % save calculated cost
    cost = zeros(K, L);
    
    % calculate cost
    for l=1:L
        cost(:, l) = G(:, l) + P(:, :, l) * J_opt;
    end
    
    % find minimum and update J_opt
    [J_opt_next, u_opt_ind_next] = min(cost, [], 2);
end

function  [G_u, P_u] = FindForU( P, G, u_opt_ind )
    % find G(:, u) and P(:, :, u)
    %
    % written by Dongho
    
    K = size(P, 1); % num of states
    
    G_u = zeros(K, 1);
    P_u = zeros(K, K);
    
    for i=1:K
        G_u(i) = G(i, u_opt_ind(i));
        P_u(i, :) = P(i, :, u_opt_ind(i));
    end
end

% LP
function [ J_opt, u_opt_ind ] = LinearProgrammingDH( P, G )
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
    
    % written by Dongho
    
    K = size(P, 1); % num of states
    L = size(P, 3); % num of controls
        
    A = zeros(K, K, L);
    B = reshape(G, K*L, 1); 
    
    for l=1:L
        % A(:,:,l) is K x K matrix
        A(:,:,l) = eye(K) - P(:,:,l);
    end
    
    % reshape A 
    A_reshape = zeros(K*L, K);
    
    for l=1:L
        % vertically concat
        A_reshape((l-1)*K+1:l*K,:) = A(:,:,l);
    end
    
    % remove inf terms
    A_reshape = A_reshape((B~=inf),:);
    B = B(B~=inf);
    
    % to make max, f should be negative
    f = -ones(K, 1);
    
    % find optimal cost-to-go
    options = optimset('Display','none');
    J_opt = linprog(f,A_reshape,B, [], [], [], [], [], options);
        
    % find optimal action
    cost = zeros(K, L);
    
    for l=1:L
        cost(:, l) = G(:, l) + P(:, :, l) * J_opt;
    end
    
    [~, u_opt_ind] = min(cost, [], 2);
end