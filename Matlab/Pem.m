function [tau,prob]= Pem(r,lim_mag,lim_time,D_mu,D_sig,Loftau,N_ref,L_min,L_max)
%[tau,prob]= Pem(r,lim_mag,lim_time,D_mu,D_sig,Loftau,N_ref,L_min,L_max)
% Pem calculates the values of the P_EM integral 
%
% r : the raidus of the aperture, in meters.
%
% lim_mag : the limiting magnitude of the telescope. For kilonovae, the
%           magnitude should be in R-band.
%
% lim_time : the limiting time at which the lim_mag is achieved, in
%            seconds. For example, if a telescope can achieve apparent an R
%            band magnitude 21 in 3 minutes, then lim_time is 180 seconds,
%            and lim_mag is 21.
%
% D_mu : this method employs a Gaussian distribution to approximate the 
%        distance information. D_mu is the mean, in par sec.
%
% D_sig : The sigma on distance, in par sec.
%
% Loftau : the length of the tau vector. The the smallest value of tau is
%          0.1s, the largest is 1e(-1+(Loftau-1)*0.1)s. 
%          If the input is [], the value will be set to 61. 
%
% N_ref : the number of photons expected at apparent magnitude equal to
%         zero. If [] is input, the value would be equal to 9.7847*10^9.  
%
% L_min : the minimum peak Luminosity of kilonovae, in w. If the input is [],
%         the value will be set equal to 4.970e31w (M = -8) from  
%
% L_max : the maximum peak Luminosity of kilonovae, in w. If the input is [],
%         the value will be set equal to 4.970e33w (M = -13) from (Barnes &
%         Kasen, 2013).


if isempty(L_min) == 1
    L_min=4.9370*10^31;
elseif isnumeric(L_min) == 0
    error('L_min has to either be an empty array of a number.')
end

if isempty(L_max) == 1
    L_max=4.9370*10^33;
elseif isnumeric(L_max) == 0
    error('L_max has to either be an empty array of a number.')
end

if isnumeric(r) == 0
    error('please enter the raidus of the aperture, in meters.')
end

if isnumeric(lim_mag) == 0
    error('please enter the limiting magnitude of the interested telescope.')
end

if isnumeric(lim_time) == 0
    error('please enter the limiting time of the interested telescope that corresponds to the limiting magnitude.')
end

if isempty(D_mu) == 1
    D_mu=200.0*10^6;
elseif isnumeric(D_mu) == 0
    error('D_mu has to either be an empty array of a number.')
end

if isempty(D_sig) == 1
    D_sig=60.0*10^6;
elseif isnumeric(D_sig) == 0
    error('D_sig has to either be an empty array of a number.')
end

if isempty(Loftau) == 1
    Loftau=61;
elseif isnumeric(Loftau) == 0
    error('Loftau has to either be an empty array of a number.')
end


if isempty(N_ref) == 1
    N_ref=9.7847*10^9;
elseif isnumeric(N_ref) == 0
    error('please enter the no. of photons expected at m = 0, or enter [] for this entry.')
end
         
sample_length=1000;
L=linspace(L_min,L_max,sample_length);                                 
logL=log10(L);
dL=L(2)-L(1);
L_sun=3.85*10^26;
logL_sun=log10(L_sun);
k=1/(2*(sqrt(L_max)-sqrt(L_min)));
L_prior=k./sqrt(L);                                                         
M=4.77-2.5.*logL+2.5*logL_sun;
pMdM=L_prior*dL;                                                            % Performaning a changing variable from peak luminosity to peak absolute magnitude

D_min=10*10^6;                                                              % The minimal distance of kilonova from earth at which the source is sampled.
D_max=1000*10^6;                                                            % The maximal distance at which the source is sampled.
Rstepsize=(D_max-D_min)/sample_length;
Rsteps=(D_max-D_min)/Rstepsize;
R=linspace(D_min,D_max,Rsteps);
p_R=normpdf(R,D_mu,D_sig);



Loftau=61;
tau=zeros(Loftau,1);
prob=zeros(Loftau,1);
N_exp=zeros(length(R),sample_length);
                                                     
scaling_num=10^(log10(N_ref)-lim_mag/2.5);                                  % this is converting the limiting magnitude in limiting time to N*/A (eq. 2).     
scaling_flux=scaling_num*pi*r^2*lim_time;                       

message1='Computing the values of P_EM at tau equal to ';
message2='This calculation will be repeated ';
message3=' times';
message4='There are still ';
message5=' times to go.';

for i=1:Loftau                                                              % this for loop computes the values of P_EM at various values of observation tau. or equation 5b in the paper.
    tau(i)=10^(-1+(i-1)*0.1);
    display_message1=[message1,num2str(tau(i)),'s'];
    display_message2=[message2,num2str(Loftau),message3];
    disp(display_message1)
    disp(display_message2)
    display_message2=[message4,num2str(Loftau-i),message5];
    disp(display_message2)
    fprintf('\n')
    for c=1:length(R)
        for d=1:sample_length
            N_exp(c,d)=10^(log10(N_ref)-(M(d)+5*log10(R(c))-5)/2.5);
            int_value=(1-cdf('poiss',scaling_flux,N_exp(c,d)*tau(i)*pi*r^2))...
                *p_R(c)*pMdM(d)*Rstepsize;
            prob(i)=prob(i)+int_value;
        end
        
    end
end

logtau=log10(tau);
logprob=log10(prob);

figure 
 plot(tau,prob,'linewidth',2)
 xlabel('Observing Time \tau','fontsize',20)
 ylabel('\it{P_{\rmEM}|}','fontsize',20)
 grid
 set(gca,'linewidth',2)

figure 
plot(logtau,logprob,'linewidth',2)
 xlabel('Log of Observing Time \tau','fontsize',20)
 ylabel('Log of Detection Prob. Regardless of Skylocation','fontsize',20)
 grid
 set(gca,'linewidth',2)