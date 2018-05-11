function [ J_opt, u_opt_ind ] = LinearProgramming( P, G )
%LINEARPROGRAMMING Value iteration
%   Solve a stochastic shortest path problem by linear programming.
%
%   [J_opt, u_opt_ind] = LinearProgramming(P, G) computes the optimal cost
%   and the optimal control input for each state of the state space.
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

% find the target index
[target_index1, target_index2]=find(P(:,:,1)==1)

% When we compute the P matrix, we set the P(target, target,:)=1, However,
% it may cause some trouble in linear programming, since J(target) will be
% unbounded. So we revise P matrix as follows, so that J(target) will be set to zero in linprog
% function
P(target_index1,:, :)=0;
P(target_index1,target_index1, 1)=1;


states_number = size(G,1);
L = size(P,3);

% set negative linear programming coefficients since linprog is a
% minimazition problem
f = -ones(states_number,1)*3;

% construct A and b for the linprog function
stack_A = [];
stack_b = [];
for i = 1 : L
    stack_A = [stack_A; eye(states_number)- P(:,:,i)];
    stack_b = [stack_b;G(:,i)];
end

% set the inf elements in stack_b to a large number so that the linprog can
% work
stack_b(find(stack_b==inf))=1e10;

% set upper and lower bound for the Js
lb = zeros(states_number,1);
ub = Inf(states_number,1);

J_opt = linprog(f,stack_A,stack_b,[],[],lb,ub);
u_opt_ind = zeros(states_number, 1);

% J_opt_ind is the optimal cost to go, it satifies the Bellman Equation
for index = 1 : states_number
    for u = 1 : L
        % for a certain policy at a certain state, if the difference between 
        % left-hand side and right-hand side of the Bellman Equation is very small, 
        % we can say this policy is the best policy at this state 
        if abs(J_opt(index) - G(index,u) - P(index, :, u)* J_opt) <= 0.01
            u_opt_ind(index)=u;
            break;
        end
    end
end

end


