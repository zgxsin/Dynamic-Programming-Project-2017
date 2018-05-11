function [ J_opt, u_opt_ind ] = PolicyIteration( P, G )
%POLICYITERATION Value iteration
%   Solve a stochastic shortest path problem by policy iteration.
%
%   [J_opt, u_opt_ind] = PolicyIteration(P, G) computes the optimal cost and
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
    num_state = size(G,1);
	for i = 1:num_state
        if (P(i,i,:) == 1)  %target state
        	target_state = i;
        	break;
        end
    end
    P(target_state,:,:) = [];     %in policy iteration, terminal state are not included
    P(:,target_state,:) = [];
    G(target_state,:) = [];
    J_opt = zeros(num_state - 1,1);
    u_opt_ind = ones(num_state - 1,1);
    u_opt_pre = zeros(num_state - 1,1);  %initialization
    
    while not(isequal(u_opt_ind, u_opt_pre))
        u_opt_pre = u_opt_ind;  %update policy
        tmp_G = zeros(num_state-1,1);
        tmp_P = zeros(num_state-1,num_state-1);
        %policy evaluation
        for i = 1:(num_state-1)
            tmp_G(i) = G(i,u_opt_pre(i));
            tmp_P(i,:) = P(i,:,u_opt_pre(i));
        end
        J_opt = inv(eye(num_state-1) - tmp_P)*tmp_G; %policy evaluation
        for i = 1:(num_state-1)
            index = [];
            for k = 1:size(G,2)
                if (sum(P(i,:,k)) > 0)  %feasible control input
                    index = [index, k];
                end
            end
            [M, ctr_opt] = min(G(i,index)' + squeeze(P(i,:,index))'*J_opt);
            u_opt_ind(i) = index(ctr_opt);
        end
    end
    disp(target_state);
    J_opt = [J_opt(1:target_state-1)', 0, J_opt(target_state:end)'];
    u_opt_ind = [u_opt_ind(1:target_state-1)', 1, u_opt_ind(target_state:end)'];
end

