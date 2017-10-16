function [time_allocated,individial_prob,possibility]=Solver2(PGW,nos,tau,prob,T)
    pp=interp1(tau,prob,'spline','pp');
    pp_val = fntlr(pp,2,tau);
    pp_val = pp_val(2,:);
function [prob_a_field,gamma]=func(X)
    num_var=numel(X);
    prob_a_field=zeros(1,num_var-1);
    for i=1:num_var-1
        prob_a_field(i)=PGW(i)*interp1(tau,prob,X(i), 'spline');
    end
    
    lambda=X(num_var);
    gamma=sum(prob_a_field)+lambda*(sum(X(1:end-1))-T);
end

function dlambda=dfunc(X)
    num_var=numel(X);
    dlambda=nan(size(X));
    
    for i=1:num_var-1
        dlambda(i)=PGW(i)*interp1(tau,pp_val,X(i))+X(end);
    end
        dlambda(end)=sum(X(1:end-1))-T;   
end

sum_prob_rl=sum(PGW(1:nos));
proposed_solution=zeros(1,length(PGW(1:nos)));

for j=1:length(PGW(1:nos))
    proposed_solution(j)=PGW(j)*T/sum_prob_rl;
end


[time_allocated,~,~]=fsolve(@dfunc,[ones(1,nos) 0],optimset('MaxFunEval',1e6,'MaxIter',8000,'Display','off','tolfun',1e-16,'tolx',1e-16));
% 
[individial_prob,possibility]=func(time_allocated);




end

