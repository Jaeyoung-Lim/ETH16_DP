function P = ComputeTransitionProbabilities( stateSpace, controlSpace, map, gate, mansion, cameras )
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

l = size(controlSpace,1);
[K, ~] = size(stateSpace);
[M, N] = size(map);
P = zeros(K, K, l); %Initialize P
u = [1, 0;
    0, -1;
    -1, 0;
    0, 1;
    0, 0];
    
    

state_gate = find(ismember(stateSpace, gate), 1);

for i=1:N
   for j=1:M % State K+1 is the termination state
       % find index of current state
       if map(j, i)<=0
        cur_pos = [j, i];
        cur_state = find(ismember(stateSpace, cur_pos), 1);
       else
        continue;
       end
        
       %Probability of getting good picture
       % Check Probability to get a good picture of the mansion
       for m=1: size(mansion, 1)
            err_mansion = [j, i]-mansion(m);
            if err_mansion(1)==0 || err_mansion(2)==0 %Camera is aligned with position
                if err_mansion(1)==0
                    %Values along the camera-position path
                    path_mansion = map(j, min(i, i-sign(err_mansion(2))*norm(err_mansion(2))):max(i, i-sign(err_mansion(2))*norm(err_mansion(2))));
                elseif err_mansion(2)==0
                    path_mansion = map(min(j, j-sign(err_mansion(1))*norm(err_mansion(1))):max(j, j-sign(err_mansion(1))*norm(err_mansion(1))), i);
                else
                    
                end
                if any(path_mansion(:)<=0) %There is no occlusion
                    P(cur_state, cur_state, 5) = 1-max(0.001, 0.5/norm(err_mansion));
                    %terminal state is not defined
                end
            end
       end
       
       for l=4:1 %Probability for moving inputs
           future_pos = cur_pos+u(l);
            if map(future_pos(1), future_pos(2)) < 0 % Check if there are obstacle
               % Check Probability to be caught on Each camera
               for c=1: size(cameras, 1)
                err_cam = future_pos-cameras(c);
                    if err_cam(1)==0 || err_cam(2)==0%Camera is aligned with position
                        if err_cam(1)==0
                            %Values along the camera-position path
                            path_cam = map(j, min(i,i-sign(err_cam(2))*norm(err_cam(2))):max(i,i-sign(err_cam(2))*norm(err_cam(2))));

                        elseif err_cam(2)==0
                            path_cam = map(min(j,j-sign(err_cam(1))*norm(err_cam(1))):max(j,j-sign(err_cam(1))*norm(err_cam(1))), i);
                        end
                        if any(path_cam(:)<=0) %There is no occlusion
                            P_notcaught=P_notcaught*(1-cameras(c, 3)*norm(err_cam));
                        end
                    end
               end              
               P(cur_state, find(ismember(stateSpace, future_pos), 1), l)=P_notcaught;
            end
        end
    end
end