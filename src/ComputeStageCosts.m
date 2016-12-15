
function G = ComputeStageCosts( stateSpace, controlSpace, map, gate, mansion, cameras)
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


% if stay: cost = 1
% if move: cost = 1
% if moving into a pool: cost = 1
% if moving out of pool: cost = 4


% put your code here

K = size(stateSpace, 1);
L = size(controlSpace, 1);
G = zeros(K, L);

% position of gate
gg = find(ismember(stateSpace, gate, 'rows'), 1); 

P = ComputeTransitionProbabilities( stateSpace, controlSpace, map, gate, mansion, cameras );

% size of the map
M = size(map, 1);
N = size(map, 2);
Ncam = size(cameras, 1);

% construct field-of-view of cameras
FOV = findFOV(cameras, map);

% states that can see the mansion
MV = findMansionView(mansion, map);

% control inputs
% 1-> n
% 2-> w
% 3-> s
% 4-> e
% 5-> p

for k = 1:K
     m = stateSpace(k,2);
     n = stateSpace(k,1);
     ii = find(ismember(stateSpace, [n, m], 'rows'), 1);

     % cost of moving north
     G(k, 1) = 0;
     if(m>=M) 
         G(k, 1) = inf;
     elseif( map(m+1, n)>0 )
         G(k, 1) = inf;
     else
         jj = find(ismember(stateSpace, [n, m+1], 'rows'), 1);
         Pcam = findPCamera(cameras, FOV, n, m+1);
         if( map(m+1, n)<0 )
             G(k,1) = G(k,1) + 4*P(ii,jj,1); % if in pool
             if( jj ~= gg )
                 G(k,1) = G(k,1) + 7*Pcam*(1-Pcam) + 8*Pcam*(1-Pcam)^2 + ...
                           9*Pcam*(1-Pcam)^3 + 10*Pcam*(1-Pcam)^4;
             else
                 G(k,1) = G(k,1) + 1*(1-Pcam) + 7*Pcam;
             end
         else    % if not in pool
             if( jj ~= gg )
                 G(k,1) = G(k,1) + 1*P(ii,jj,1) + 7*P(ii,gg,1);
             else
                 G(k,1) = G(k,1) + 1*(1-Pcam) + 7*Pcam;
             end
         end
     end

     % cost of moving west
     G(k, 2) = 0;
     if(n<=1) 
         G(k, 2) = inf;
     elseif( map(m, n-1)>0 )
         G(k, 2) = inf;
     else
         jj = find(ismember(stateSpace, [n-1, m], 'rows'), 1);
         Pcam = findPCamera(cameras, FOV, n-1, m);
         if( map(m, n-1)<0 )
             G(k,2) = G(k, 2) + 4*P(ii,jj,2); % if in pool
             if( jj ~= gg )
                 G(k,2) = G(k, 2) + 7*Pcam*(1-Pcam) + 8*Pcam*(1-Pcam)^2 + ...
                           9*Pcam*(1-Pcam)^3 + 10*Pcam*(1-Pcam)^4;
             end
         else    % if not in pool
             if( jj ~= gg )
                 G(k,2) = G(k, 2) + 1*P(ii,jj,2) + 7*P(ii,gg,2);
             else
                 G(k,2) = G(k, 2) + 1*(1-Pcam) + 7*Pcam;

             end
         end
     end

     % cost of moving south
     G(k, 3) = 0;
     if(m<=1) 
         G(k, 3) = inf;
     elseif( map(m-1, n)>0 )
         G(k, 3) = inf;
     else
         jj = find(ismember(stateSpace, [n, m-1], 'rows'), 1);
         Pcam = findPCamera(cameras, FOV, n, m-1);
         if( map(m-1, n)<0 )
             G(k,3) = G(k,3) + 4*P(ii,jj,3); % if in pool
             if( jj ~= gg )
                 G(k,3) = G(k,3) + 7*Pcam*(1-Pcam) + 8*Pcam*(1-Pcam)^2 + ...
                           9*Pcam*(1-Pcam)^3 + 10*Pcam*(1-Pcam)^4;
             end
         else    % if not in pool
             if( jj ~= gg )
                 G(k,3) = G(k,3) + 1*P(ii,jj,3) + 7*P(ii,gg,3);
             else
                 G(k,3) = G(k,3) + 1*(1-Pcam) + 7*Pcam;

             end
         end
     end

     % cost of moving east
     G(k, 4) = 0;
     if(n>=N) 
         G(k, 4) = inf;
     elseif( map(m, n+1)>0 )
         G(k, 4) = inf;
     else
         jj = find(ismember(stateSpace, [n+1, m], 'rows'), 1);
         Pcam = findPCamera(cameras, FOV, n+1, m);
         if( map(m, n+1)<0 )
             G(k,4) = G(k,4) + 4*P(ii,jj,4); % if in pool
             if( jj ~= gg )
                 G(k,4) = G(k,4) + 7*Pcam*(1-Pcam) + 8*Pcam*(1-Pcam)^2 + ...
                           9*Pcam*(1-Pcam)^3 + 10*Pcam*(1-Pcam)^4;
             end 
         else    % if not in pool
             if( jj ~= gg )
                 G(k,4) = G(k,4) + 1*P(ii,jj,4) + 7*P(ii,gg,4);
             else
                 G(k,4) = G(k,4) + 1*(1-Pcam) + 7*Pcam;

             end
         end
     end

     % cost of taking a photo
     if( ii ~= gg )
         G(k, 5) = 1*P(ii,ii,5);                 % if fail but not caught
         G(k, 5) = G(k, 5) + 7*P(ii,gg,5);                % if go to gate
         G(k, 5) = G(k, 5) + 1*(1-P(ii,ii,5)-P(ii,gg,5)); % if success
     else
         % if take photo at gate, special treatment
         % assume you cannot see the mansion from the gate
         Pcam = findPCamera(cameras, FOV, n, m);
         Ps = findPSuccess(mansion, MV, n, m);
         G(k, 5) = (1-Ps)*(1*(1-Pcam) + 7*Pcam) + 1*Ps;
     end

end

end


function FOV = findFOV(cameras, map)
%
% construct field-of-view of cameras
%
%       map:
%           A (M x N)-matrix describing the terrain of the estate map.
%           Positive values indicate cells that are inaccessible (e.g.
%           trees, bushes or the mansion) and negative values indicate
%           ponds or pools.
%
%    	cameras:
%          	A (H x 3)-matrix indicating the positions and quality of the 
%           cameras.
%
%       FOV: field of view of cameras. Rows of FOV{i} contains (n, m) coordinates
%            of points in the field of view of the i-th camera
% 

    FOV{1} = [0,0];
    M = size(map, 1);
    N = size(map, 2);
    Ncam = size(cameras, 1);

    for c = 1:Ncam
        m = cameras(c, 2);
        n = cameras(c, 1);
        
        FOV_tmp = [0, 0];

        % north
        ii = 1;
        while( (m+ii <= M) && (map(m+ii, n)<=0))
            FOV_tmp = [FOV_tmp; [n, m+ii]];
            ii = ii+1;
        end

        % south
        ii = 1;
        while( (m-ii >= 1) && (map(m-ii, n)<=0))
            FOV_tmp = [FOV_tmp; [n, m-ii]];
            ii = ii+1;
        end
    
        % west
        ii = 1;
        while( (n-ii >= 1) && (map(m, n-ii)<=0))
            FOV_tmp = [FOV_tmp; [n-ii, m]];
            ii = ii+1;
        end

        % east
        ii = 1;
        while( (n+ii <= N) && (map(m, n+ii)<=0))
            FOV_tmp = [FOV_tmp; [n+ii, m]];
            ii = ii+1;
        end

        FOV{c} = FOV_tmp(2:end, :);
    end
end


function MV = findMansionView(mansion, map)
%
% find states where you can see the mansion
%

    M = size(map, 1);
    N = size(map, 2);
    Nm = size(mansion, 1);

    MV_tmp = [0, 0];
    for c = 1:Nm
        m = mansion(c, 2);
        n = mansion(c, 1);

        % north
        ii = 1;
        while( (m+ii <= M) && (map(m+ii, n)<=0))
            MV_tmp = [MV_tmp; [n, m+ii]];
            ii = ii+1;
        end

        % south
        ii = 1;
        while( (m-ii >= 1) && (map(m-ii, n)<=0))
            MV_tmp = [MV_tmp; [n, m-ii]];
            ii = ii+1;
        end
    
        % west
        ii = 1;
        while( (n-ii >= 1) && (map(m, n-ii)<=0))
            MV_tmp = [MV_tmp; [n-ii, m]];
            ii = ii+1;
        end

        % east
        ii = 1;
        while( (n+ii <= N) && (map(m, n+ii)<=0))
            MV_tmp = [MV_tmp; [n+ii, m]];
            ii = ii+1;
        end
    end
    MV = MV_tmp(2:end, :);
end



function Pcam = findPCamera(cameras, FOV, n, m)
%
% find chance of geting caught by camera(s)
%
%    	cameras:
%          	A (H x 3)-matrix indicating the positions and quality of the 
%           cameras.
%
%       FOV: field of view of cameras, computed using findFOV function
%
%       n, m: coordinates of the state to check
% 
%       Pcam: probability of being caught by the camera(s)
%

     Ncam = size(cameras, 1);
     Pnot = 1;
     for c = 1:Ncam
         if( sum( ismember(FOV{c}, [n,m], 'rows'))>0 )
             dist = sqrt((cameras(c,1)-n)^2 + (cameras(c,2)-m)^2);
             Pnot *= (1-cameras(c,3)/dist);
         end
     end

     Pcam = 1-Pnot;
end




function Ps = findPSuccess(mansion, MV, n, m)
%
% find chance of taking a good picture
%
%     MV: mansion view, computed with findMansionView function
%

  Ps = 0.001;

    % if we can see mansion directly
    if(sum( ismember(MV, [n,m], 'rows'))>0)
        dist_vec = sqrt( (mansion(:,1)-n).^2 + (mansion(:,2)-m).^2);
        dist_min = min(dist_vec);
        Ps = max(0.001, 0.5/dist_min);
    end
end
