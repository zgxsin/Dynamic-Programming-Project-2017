function [ J_opt, u_opt_ind ] = ValueIteration( P, G )
%VALUEITERATION Value iteration
%   Solve a stochastic shortest path problem by value iteration.
%
%   [J_opt, u_opt_ind] = ValueIteration(P, G) computes the optimal cost and
%   the optimal control input for each state of the state space.
%
%   Input arguments:
%
%       P:
%           A (MN x MN x L) matrix containing the transition probabilities
%           between all states in the state space for all attainable
%           control inputs. The entry P(i, j, l) represents the transition
%           probability from state i to state j if control input l is
%           applied.
%
%       G:
%           A (MN x L) matrix containing the stage costs of all states in
%           the state space for all attainable control inputs. The entry
%           G(i, l) represents the cost if we are in state i and apply
%           control input l.
%
%   Output arguments:
%
%       J_opt:
%       	A (1 x MN) matrix containing the optimal cost-to-go for each
%       	element of the state space.
%
%       u_opt_ind:
%       	A (1 x MN) matrix containing the indices of the optimal control
%       	inputs for each element of the state space.

% put your code here
    num_iteration = 1000;   % 1000 times of updating
    num_state = size(G,1);
    J_opt = zeros(1,num_state);
    J_tmp = zeros(1,num_state);
    u_opt_ind = zeros(1,num_state);
    
    for j = 1:num_state
    	if (P(j,j,:) == 1)  %target state
        	u_opt_ind(j) = 1;
            continue;
        end
        index = [];
        for k = 1:size(G,2)
        	if (sum(P(j,:,k)) > 0)
            	index = [index, k];
            end
        end
        [J_tmp(j), ctr_opt] = min( G(j,index) + J_opt*squeeze(P(j,:,index)) );
        u_opt_ind(j) = index(ctr_opt);
    end
    
    while(sum(abs(J_tmp - J_opt))> 0.000001)
        J_opt = J_tmp;
        for j = 1:num_state
            if (P(j,j,:) == 1)  %target state
                u_opt_ind(j) = 1;
                continue;
            end
            index = [];
            for k = 1:size(G,2)
                if (sum(P(j,:,k)) > 0)
                    index = [index, k];
                end
            end
            [J_tmp(j), ctr_opt] = min( G(j,index) + J_opt*squeeze(P(j,:,index)) );
            u_opt_ind(j) = index(ctr_opt);
        end
    end
end

