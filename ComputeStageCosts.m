function G = ComputeStageCosts( stateSpace, controlSpace, mazeSize, walls, targetCell, holes, resetCell, p_f, c_p, c_r )
%COMPUTESTAGECOSTS Compute stage costs.
% 	Compute the stage costs for all states in the state space for all
%   attainable control inputs.
%
%   G = ComputeStageCosts(stateSpace, controlSpace, disturbanceSpace,
%   mazeSize, walls, targetCell) computes the stage costs for all states in
%   the state space for all attainable control inputs.
%
%   Input arguments:
%
%       stateSpace:
%           A (MN x 2) matrix, where the i-th row represents the i-th
%           element of the state space. Note that the state space also
%           contains the target cell, in order to simplify state indexing.
%
%       controlSpace:
%           A (L x 2) matrix, where the l-th row represents the l-th
%           element of the control space.
%
%       mazeSize:
%           A (1 x 2) matrix containing the width and the height of the
%           maze in number of cells.
%
%   	walls:
%          	A (2K x 2) matrix containing the K wall segments, where the start
%        	and end point of the k-th segment are stored in row 2k-1
%         	and 2k, respectively.
%
%    	targetCell:
%          	A (2 x 1) matrix describing the position of the target cell in
%         	the maze.
%       holes:
%         	A (H x 2) matrix containg the H holes of the maze. Each row
%         	represents the position of a hole.
%
%   	resetCell:
%         	A (1 x 2) matrix describing the position of the reset cell in
%           the maze.
%
%       p_f:
%           The probability of falling into a hole when the ball is
%           traversing through or to a cell with a hole
%       
%       c_p:
%           Every time the ball bounces into a wall or boundary, we get this number 
%            of time steps as penalty.
%       c_r:
%           Every time the ball falls into a hole, the ball is set to the reset cell
%           at the beginning of the next stage and we get this number of time steps
%           as additional penalty.
%
%   Output arguments:
%
%       G:
%           A (MN x L) matrix containing the stage costs of all states in
%           the state space for all attainable control inputs. The entry
%           G(i, l) represents the cost if we are in state i and apply
%           control input l.

% put your code here
    hole_state = zeros(size(holes,1),1);
    for i = 1:size(holes,1)
        hole_state(i) = (holes(i,1) - 1) * mazeSize( 2 ) + holes(i,2);
    end
    reset_state = (resetCell(1) - 1) * mazeSize( 2 ) + resetCell(2);    %corresponding state of resetCell
    num_state = size(stateSpace,1);
    num_strat = size(controlSpace,1);
    G = Inf(num_state,num_strat);   %initialize state cost matrix
    
    for state = 1:num_state
        %if ball is at targetCell, no stage cost
        if (isequal(stateSpace(state,:),targetCell))
            G(state,:) = 0;
            continue;
        end
        fea_ctr = [];          %feasible control strategy for current state
        % decide available control strategies for the current state
        for i = 1:num_strat    
            feasibility = AvailableControl(stateSpace(state,:), controlSpace(i,:),walls, mazeSize);
            if (feasibility)
                fea_ctr = [fea_ctr; controlSpace(i,:), i];
            end
        end  
        raw_n = stateSpace(state,1);
        raw_m = stateSpace(state,2);% get present position
        for j = 1:size(fea_ctr,1)
            label = fea_ctr(j,3);   %label of strategy
            G(state,label) = 0;     %feasible solution must have stage cost
            mid_n = raw_n + fea_ctr(j,1);   %position after control input while before disturbance
            mid_m = raw_m + fea_ctr(j,2);
            mid_state = (mid_n - 1) * mazeSize( 2 ) + mid_m;
            num_hole = TraverseHole(stateSpace(state,:),fea_ctr(j,1:2), holes, mazeSize);   %number of holes along the control direction
            prob_trav = (1 - p_f)^num_hole; %probability of not returning to reset state after control input while before disturbance
            G(state,label) = (1 - prob_trav)*c_r;     %cost from falling into holes when traversing
            % calculating cost in disturbance process
            num_bounce = 0;
            for k = 2:9
                if (AvailableControl([mid_n,mid_m], controlSpace(k,:),walls, mazeSize)) %check whether disturbance hit walls
                    end_state = (mid_n + controlSpace(k,1) - 1) * mazeSize( 2 ) + mid_m + controlSpace(k,2);
                    if (find(end_state == hole_state(:)))   %check whether endstate is a hole
                        G(state,label) = G(state,label) + c_r*prob_trav/9*p_f;  %cost caused by falling into holes after disturbance, except bounced back to mid_state
                    end
                else
                    num_bounce = num_bounce + 1;
                end
            end
            G(state,label) = G(state,label) + c_p*num_bounce*prob_trav/9;       %cost caused by bouncing into walls
            %if mid_state is a hole
            if (find(mid_state == hole_state(:)))
                G(state,label) = G(state,label) + c_r*num_bounce*prob_trav*p_f/9;
            end
            G(state,label) = G(state,label) + 1;        %time cost of every step
        end
    end
end


function feasibility = AvailableControl(state, control, walls, mazeSize) % check whether this strategy is available
    feasibility = 1;    %initialization
    num_walls = size(walls,1)/2;
    n = state(1);
    m = state(2);%get position
    step_n = control(1);
    step_m = control(2);%get step
    
	if (step_n == 0 && step_m == 0)
    elseif (step_n == 1 && step_m == 0)%step(1,0)
        if (n+1 > mazeSize(1))
        	feasibility = 0;    %present position is on the boundary, not feasible
        else     
            for k = 1:num_walls
                if (isequal([n,m-1;n,m],walls(2*k-1:2*k,:)) || isequal([n,m;n,m-1],walls(2*k-1:2*k,:)))
                    feasibility = 0;    %moving direction has wall, not feasible
                	break;     
                end
            end
        end
	elseif (step_n == 0 && step_m == 1)
        if (m+1 > mazeSize(2))
        	feasibility = 0;    %present position is on the boundary, not feasible
        else     
            for k = 1:num_walls
                if (isequal([n-1,m;n,m],walls(2*k-1:2*k,:)) || isequal([n,m;n-1,m],walls(2*k-1:2*k,:)))
                    feasibility = 0;    %moving direction has wall, not feasible
                	break;    
                end
            end
        end
	elseif (step_n == -1 && step_m == 0)
        if (n-1 < 1)
        	feasibility = 0;    %present position is on the boundary, not feasible
        else     
            for k = 1:num_walls
                if (isequal([n-1,m-1;n-1,m],walls(2*k-1:2*k,:)) || isequal([n-1,m-1;n-1,m],walls(2*k-1:2*k,:)))
                    feasibility = 0;    %moving direction has wall, not feasible
                	break;    
                end
            end
        end
	elseif (step_n == 0 && step_m == -1)
        if (m-1 < 1)
        	feasibility = 0;    %present position is on the boundary, not feasible
        else     
            for k = 1:num_walls
                if (isequal([n-1,m-1;n,m-1],walls(2*k-1:2*k,:)) || isequal([n,m-1;n-1,m-1],walls(2*k-1:2*k,:)))
                    feasibility = 0;    %moving direction has wall, not feasible
                	break;    
                end
            end
        end
	elseif (step_n == 2 && step_m == 0)
    	if (n+2 > mazeSize(1))
        	feasibility = 0;    %present position is on the boundary, not feasible
        else
            for k = 1:num_walls
                if (isequal([n,m-1;n,m],walls(2*k-1:2*k,:)) || isequal([n,m;n,m-1],walls(2*k-1:2*k,:))...
                    || isequal([n+1,m-1;n+1,m],walls(2*k-1:2*k,:)) || isequal([n+1,m;n+1,m-1],walls(2*k-1:2*k,:)))
                    feasibility = 0;    %moving direction has wall, not feasible
                	break;     
                end
            end
        end
	elseif (step_n == 0 && step_m == 2)
    	if (m+2 > mazeSize(2))
        	feasibility = 0;    %present position is on the boundary, not feasible
        else
            for k = 1:num_walls
                if (isequal([n-1,m;n,m],walls(2*k-1:2*k,:)) || isequal([n,m;n-1,m],walls(2*k-1:2*k,:))...
                    || isequal([n-1,m+1;n,m+1],walls(2*k-1:2*k,:)) || isequal([n,m+1;n-1,m+1],walls(2*k-1:2*k,:)))
                    feasibility = 0;    %moving direction has wall, not feasible
                	break;    
                end
            end
        end
	elseif (step_n == -2 && step_m == 0)
    	if (n-2 <1)
        	feasibility = 0;    %present position is on the boundary, not feasible
        else     
            for k = 1:num_walls
                if (isequal([n-1,m-1;n-1,m],walls(2*k-1:2*k,:)) || isequal([n-1,m-1;n-1,m],walls(2*k-1:2*k,:))...
                    || isequal([n-2,m-1;n-2,m],walls(2*k-1:2*k,:)) || isequal([n-2,m;n-2,m-1],walls(2*k-1:2*k,:)))
                    feasibility = 0;    %moving direction has wall, not feasible
                	break;    
                end
            end
        end        
	elseif (step_n == 0 && step_m == -2)
        if (m-2 < 1)
        	feasibility = 0;    %present position is on the boundary, not feasible
        else     
            for k = 1:num_walls
                if (isequal([n-1,m-1;n,m-1],walls(2*k-1:2*k,:)) || isequal([n,m-1;n-1,m-1],walls(2*k-1:2*k,:))...
                    || isequal([n-1,m-2;n,m-2],walls(2*k-1:2*k,:)) || isequal([n,m-2;n-1,m-2],walls(2*k-1:2*k,:)))
                    feasibility = 0;    %moving direction has wall, not feasible
                	break;    
                end
            end
        end
    elseif (step_n == 1 && step_m == 1)
        if (n+1 > mazeSize(1) || m+1 > mazeSize(2))
            feasibility = 0;
        else
            for k = 1:2*num_walls
                if(isequal([n,m],walls(k,:)))
                    feasibility = 0;
                end
            end
        end
    elseif (step_n == 1 && step_m == -1)
        if (n+1 > mazeSize(1) || m-1 < 1)
            feasibility = 0;
        else
            for k = 1:2*num_walls
                if(isequal([n,m-1],walls(k,:)))
                    feasibility = 0;
                end
            end
        end
    elseif (step_n == -1 && step_m == 1)
        if (n-1 < 1 || m+1 > mazeSize(2))
            feasibility = 0;
        else
            for k = 1:2*num_walls
                if(isequal([n-1,m],walls(k,:)))
                    feasibility = 0;
                end
            end
        end
    elseif (step_n == -1 && step_m == -1)
        if (n-1 < 1 || m-1 < 1)
            feasibility = 0;
        else
            for k = 1:2*num_walls
                if(isequal([n-1,m-1],walls(k,:)))
                    feasibility = 0;
                end
            end
        end
    elseif (step_n == 2 && step_m == 2)
        if (n+2 > mazeSize(1) || m+2 > mazeSize(2))
            feasibility = 0;
        else
            for k = 1:2*num_walls
                if(isequal([n,m],walls(k,:)) || isequal([n+1,m+1],walls(k,:)))
                    feasibility = 0;
                end
            end
        end
    elseif (step_n == 2 && step_m == -2)  
        if (n+2 > mazeSize(1) || m-2 < 1)
            feasibility = 0;
        else
            for k = 1:2*num_walls
                if(isequal([n,m-1],walls(k,:)) || isequal([n+1,m-2],walls(k,:)))
                    feasibility = 0;
                end
            end
        end
    elseif (step_n == -2 && step_m == 2)
        if (n-2 < 1 || m+2 > mazeSize(2))
            feasibility = 0;
        else
            for k = 1:2*num_walls
                if(isequal([n-1,m],walls(k,:)) || isequal([n-2,m+1],walls(k,:)))
                    feasibility = 0;
                end
            end
        end
    elseif (step_n == -2 && step_m == -2)
        if (n-2 < 1 || m-2 < 1)
            feasibility = 0;
        else
            for k = 1:2*num_walls
                if(isequal([n-1,m-1],walls(k,:)) || isequal([n-2,m-2],walls(k,:)))
                    feasibility = 0;
                end
            end
        end
    end
end


function num_hole =TraverseHole(state, control, holes, mazeSize) %check number of holes along the control direction
    num_hole = 0;   %initialization
    if (control(1) == 0 && control(2) == 0) %stay at previous position
    elseif (control(1) == 0)    %move along y axis
        for i = 1:abs(control(2))
            for j = 1:size(holes,1) %check whether traverse a hole
                if (isequal([state(1),state(2)+i*sign(control(2))],holes(j,:)))
                    num_hole = num_hole + 1;
                    break;
                end
            end
        end
    elseif (control(2) == 0)    %move along x axis
        for i = 1:abs(control(1))
            for j = 1:size(holes,1) %check whether end at a hole
                if (isequal([state(1)+i*sign(control(1)),state(2)],holes(j,:)))
                    num_hole = num_hole + 1;
                    break;
                end
            end
        end        
    else    %move along diag
        for i = 1:abs(control(1))
            for j = 1:size(holes,1) %check whether end at a hole
                if (isequal([state(1)+i*sign(control(1)),state(2)+i*sign(control(2))],holes(j,:)))
                    num_hole = num_hole + 1;
                    break;
                end
            end
        end         
    end  
end