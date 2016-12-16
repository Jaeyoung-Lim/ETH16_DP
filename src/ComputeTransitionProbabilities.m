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
%% Initialize
L = size(controlSpace,1);
K = size(stateSpace, 1);

M = size(map, 1);
N = size(map, 2);

P = zeros(K, K, L); %Initialize P

u = [0, 1;
    -1, 0;
    0, -1;
    1, 0;
    0, 0]; % Control space to postion inputs
    
gg = find(ismember(stateSpace, gate, 'rows')); % Find state of gate

%% Compute for each position
for k=1:K %Iteration for steps
    
    cur_pos = stateSpace(k, :);
    cur_state = k;
   %% Compute Probability of getting good picture
   P_picture =0.001;
   for m=1: size(mansion, 1)
        err_mansion = cur_pos-mansion(m, :);
        if err_mansion(1)==0 || err_mansion(2)==0      
            if err_mansion(1) == 0  %Camera is aligned with current position
                % Values along the camera-position path
                path_mansion = map(min(cur_pos(2), cur_pos(2)-err_mansion(2)+sign(err_mansion(2))): max(cur_pos(2), cur_pos(2)-err_mansion(2)+sign(err_mansion(2))), cur_pos(1));
            elseif err_mansion(2) == 0
                path_mansion = map(cur_pos(2), min(cur_pos(1), cur_pos(1)-err_mansion(1)+sign(err_mansion(1))): max(cur_pos(1), cur_pos(1)-err_mansion(1)+sign(err_mansion(1))));
            else
                path_mansion = [1, 1];
            end
            
            if all(path_mansion(:)<=0) %There is no occlusion
                P_picture = max(0.001, 0.5/norm(err_mansion));
                %terminal state is not defined
            end      
        end
   end

   %% Compute probabilty for moving to the next cell
   for l=1:L % Iterate for each inputs
    future_pos = cur_pos+u(l, :);
    
    if future_pos(1)<=0 || future_pos(1)>N || future_pos(2)<=0 || future_pos(2)>M
        future_pos = cur_pos;
    end
    
    future_state = find(ismember(stateSpace, future_pos, 'rows'));
    
    if map(future_pos(2), future_pos(1))>0 % Check if there is an obstacle  
        P(k, future_state, l)=0;% Cannot move to a new cell which is occupied
        future_state = cur_state;
        future_pos = cur_pos;
    end
    %% Check Probability to be caught on Each camera
    
    P_notcaught =1;
    
    for c=1: size(cameras, 1) %Iterate for each camera
        err_cam = future_pos-cameras(c, 1:2);
        if err_cam(1)==0 || err_cam(2)==0 %Camera is aligned with position
            if err_cam(1)==0
                %Patch cam has values along the camera-position path
                path_cam = map(min(future_pos(2), future_pos(2)-err_cam(2)+sign(err_cam(2))):max(future_pos(2), future_pos(2)-err_cam(2)+sign(err_cam(2))), future_pos(1));
            elseif err_cam(2)==0
                path_cam = map(future_pos(2), min(future_pos(1), future_pos(1)-err_cam(1)+sign(err_cam(1))):max(future_pos(1), future_pos(1)-err_cam(1)+sign(err_cam(1))))';
            end   
            if all(path_cam(:)<=0) %There is no occlusion
                if map(future_pos(2), future_pos(1))<0 && l~=5
                    P_notcaught=P_notcaught*((1-cameras(c, 3)/norm(err_cam))^4);
                else
                    P_notcaught=P_notcaught*(1-cameras(c, 3)/norm(err_cam));
                end
            end
        end
    end

     if l==5 %When taking a picture
         if future_state == gg
             P(k, future_state, l)=1-P_picture;         
         else
            P(k, future_state, l)=P_notcaught*(1-P_picture);
            P(k, gg, l)=(1-P_notcaught)*(1-P_picture);% Failed to take picture and got caught by the camera
         end
    else %Wwhen moving to a new state
        if future_state == gg
            P(k, future_state, l)=P_notcaught+(1-P_notcaught);
        else
            P(k, future_state, l)=P_notcaught;
            P(k, gg, l)=1-P_notcaught;
        end
    end
   end
end
