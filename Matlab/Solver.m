function [time_allocated,individual_prob,probability]=Solver(PGW,nos,tau,prob,T)
%[time_allocated,individual_prob,probability]=Solver(PGW,nos,tau,prob,T)
% Solver returns the detection probability, and the time allocations 
% 
% area_index is a vector containing the values of P_GW that is in.
% descending order of the values of P_GW.
%
% nos is the current number of fields being considered.
%
% tau, prob are the data for P_EM.
% 
% T is the total observation time excluding the slew time or readout time
% or both.
%
% time_allocated is the cell array that contains the time allocations.
% For example, time_allocated{10} is the time allocation for the highest 
% 10 fields in terms of P_GW when only these 10 fields are being observed.
%
% individual_prob is the detection prob. probability rendered by an
% individual field.
function [prob_a_field,probability]=func(X)
    prob_a_field=zeros(numel(X),1);
    for i=1:length(X)
        prob_a_field(i)=PGW(i)*interp1(tau,prob,X(i),'spline');
    end
   probability=sum(prob_a_field);
   
end

time_allocated=T/nos*ones(1,nos);
[individual_prob,probability]=func(time_allocated);
end

