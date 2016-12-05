function G = ComputeStageCosts( stateSpace, controlSpace, map, gate, mansion, cameras, P )
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

% size of the map
M = size(map, 1);
N = size(map, 2);

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
         if( map(m+1, n)<0 )
             G(k,1) += 4*P(ii,jj,1); % if in pool
             if( jj != gg )
                 G(K,1) += 11*P(ii,gg,1); % if next state is gate, then cannot be escorted 
             end
         else
             G(k,1) += 1*P(ii,jj,1); % if not in pool
             if( jj != gg )
                 G(K,1) += 7*P(ii,gg,1);
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
         if( map(m, n-1)<0 )
             G(k,2) += 4*P(ii,jj,2); % if in pool
             if( jj != gg )
                 G(K,2) += 11*P(ii,gg,2);
             end
         else
             G(k,2) += 1*P(ii,jj,2); % if not in pool
             if( jj != gg )
                 G(K,2) += 7*P(ii,gg,2);
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
         if( map(m-1, n)<0 )
             G(k,3) += 4*P(ii,jj,3); % if in pool
             if( jj != gg )
                 G(K,3) += 11*P(ii,gg,3);
             end
         else
             G(k,3) += 1*P(ii,jj,3); % if not in pool
             if( jj != gg )
                 G(K,3) += 7*P(ii,gg,3);
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
         if( map(m, n+1)<0 )
             G(k,4) += 4*P(ii,jj,4); % if in pool
             if( jj != gg )
                 G(K,4) += 11*P(ii,gg,4); 
             end
         else
             G(k,4) += 1*P(ii,jj,4); % if not in pool
             if( jj != gg )
                 G(K,4) += 7*P(ii,gg,4);
             end
         end
     end

     % cost of taking a photo
     G(k, 5) = 1*P(ii,ii,5);                 % if fail but not caught
         G(k, 5) += 7*P(ii,gg,5);                % if go to gate
         G(k, 5) += 1*(1-P(ii,ii,5)-P(ii,gg,5)); % if success

endfor

end
