function [tprob,time_allocation]= Main(T_tot,readout_time,slew_time,PGW,tau,prob)
%[tprob,time_allocation]= Main(T_tot,readout_time,slew_time,PGW,tau,prob)
% Main is the main body of the optimisation code
%
% T_tot: total available telescope time
%
% readout_time : the time needed to perform readout, scalar.
%
% slew time : the time needed to slew telescope , scalar.
%             If the time needed to perform both slew and readout is just the larger
%             of the two, just set one of these parameters equal to the larger one
%             and the other equal to zero.
%
%
% area_index : A matrix representing the fields and containing the GW
%             probabilities. an output of function 'fields'. 
%
% tau : A matrix containing the observation tau from 0.1 seconds to 1e5
%      seconds. An output of Pem.
%
% prob : A matrix containing the values of Pem at each of the values of
%       tau. An output of Pem.             
%
% T : The total Observation Time available, in second.

if isnumeric(T_tot) == 0
    error('please enter the total observation times, in hours.')
end
nof = length(PGW);
message1='The Number of Fields is ';
message1=[message1,num2str(nof)];
disp(message1)
fprintf('\n')

numberoffields=linspace(1,nof,nof);
tprob=zeros(1,nof);
prob_distri=cell(1,nof);
time_allocation=cell(1,nof);
message1='Start Generating the Optimal Strategy....';
disp(message1)
fprintf('\n')




for nos=1:nof
    T=T_tot-(readout_time+slew_time)*nos;                                   % This computes the total time available to observation only. 
    
    text1=['The number of fields being considered now is ', num2str(nos)];
    disp(text1)
    
	if length(nonzeros(PGW))>=nos                                    % In the case where some of the chosen fields actually contains zero GW probability, 
                                                                            % there would be therefore no need to perform further analysis. Therefore, only when
                                                                            % the number of fields whose GW probability are not zero is larger than the no. of
                                                                            % fields for which the analyis will be carried out, then the ananlysis will continue.    
                                                                            
        [time_allocated,individial_prob,probability]=...                    % calling solver to find the solution.
            Solver(PGW(1:nos),nos,tau,prob,T);
                                                      
        text1=['The exit condition of fsolve for ', num2str(nos), ...
            ' fields is:'];
        disp(text1)
        fprintf('\n')
    else
        text1=['Because the posterior probabilities of the fields from the '...
            ,num2str(nos), ' th field onwards are zero, this analysis aborts.'];
        disp(text1)
        break
    end
    
    tprob(nos)=probability;                                                 % the total detection prob. if the fields are observed.            
    
    prob_distri{nos}=individial_prob;                                       % the individual detection prob. from each of the observed fields.
    
    time_allocation{nos}=time_allocated;                                    % the time allocation. 

    
    text1=['Provided observing ', num2str(nos),...
        ' fields, the detection prob. of the kilonova is ',...
        num2str(tprob(nos))];
    disp(text1)
    fprintf('\n')
    
    text1='The time assigned to each of these fields is as follow:';
    disp(text1)
    disp(time_allocated(1:nos))
    fprintf('\n')
end
[m,i]=max(tprob);                                                           % finding the maximum detection prob. and the optimal number of observed fields. 
text1='Therefore, based on the analyses above';
disp(text1)
fprintf('\n')
text1=['The highest possibility could be obtained by observing ',...
    num2str(i), ' fields, and the possibility is '];
disp(text1)
disp(num2str(m))
fprintf('\n')
text1='The time should be allocated as follow: ';
disp(text1)
disp(time_allocation{i})

figure
plot(numberoffields,tprob,'b','linewidth',2)
xlabel('The Total Observed Fields k','fontsize',20)
ylabel('\it{P(D_{\rmEM}|k)}','fontsize',20)
grid
set(gca,'fontsize',20)
set(gca,'linewidth',2)



















