%==========================================================================
% Date: 12/21/2021
%==========================================================================
% close all;
clear all;
disp('_________________IRS and Relay Combination________________________');
NS=input('Number of samples in 1000 ')*1000; % taking inputs
%==========================================================================
alpha_array = [0.5]; %[0 0.5 0.5]; TS factor array
Psi_ps_array = [0.5]; %[0.5 0 0.5]; %PS factor array
Tc = 1; %coherence time

eff = 0.8;%0.5; % Energy harvesting efficiency



%% system parameters
N=1;    % No. of IRSs
L=512;   % No. of elements
M=1;    % No. of BS antennas
K=1;    % No. of users 
Nr=1;   % No. of relays
sigma_N=1;  % noise power
eta=1;  % reflection coefficient

%--------------------------------------------------------------------------
% bandwidth = 20e6;       % channel bandwidth
% noiseFigure = 10;       % noise figure at BS (in dB)
% noiseVariancedBm = -174 + 10*log10(bandwidth) + noiseFigure;
% noiseVariance = 10^(0.1*noiseVariancedBm);
%--------------------------------------------------------------------------
for kk=1:length(alpha_array)
    %% variables
    alpha = alpha_array(kk); %TS factor
    Psi_ps = Psi_ps_array(kk); %power scaling factor
%     P=logspace(-1,3,20); 
P=linspace(0,100,30);
    gamma_bar=P;
%     gamma_bar=P/noiseVariance;
%% path-loss - correlation,large scale fading with shadowing
[Path_loss_S_IRS,Path_loss_IRS_R,Path_loss_R_D]=PathLoss_model(M,K,N,Nr);
%% simulation
Eh_sim_exct=[]; %--Energy harvest vector (like rate vector)
Eh_Analytical=[];

for pp=1:length(P)
    Ptx=P(pp);
    
    Eh_count_exct=0; %--EH count
    
    for ns=1:NS %number of realization = NS
% actual channels----------------------------------------------------------
        alpha_h=abs(sqrt(Path_loss_S_IRS).*(randn(M,L)+1j*randn(M,L))/sqrt(2));
        alpha_g=abs(sqrt(Path_loss_IRS_R).*(randn(Nr,L)+1j*randn(Nr,L))/sqrt(2));
        alpha_f=abs(sqrt(Path_loss_R_D).*(randn(K,Nr)+1j*randn(K,Nr))/sqrt(2));
% Assuming perfect phase cancelation---------------------------------------
         amp1=sum(eta.*alpha_h.*alpha_g,2);
         amp2=sqrt(Psi_ps)*amp1;
         
         %?????????????????????????????????????????????????????????????????
         % Optimum relay gain assuming noise power equal to 1 (Probably NEEDS some CHANGE)
         G_opt = sqrt(gamma_bar(pp)./(gamma_bar(pp).*(amp1).^2+1)); 
         %?????????????????????????????????????????????????????????????????
         
         power=eff*Ptx*alpha*Tc*(G_opt.*alpha_f.*amp1)^2+eff*Ptx*(1-alpha)*Tc*(G_opt.*alpha_f.*amp2)^2;
         
         Eh_count_exct=Eh_count_exct+power;

    end

% EH calculation---------------------------------------------------------

    Eh_sim_exct=[Eh_sim_exct (Eh_count_exct/NS)];
    
% % analysis
    Mu_X=0;
    Var_X=0;
    for ll=1:L
        lambda_l_sqr=eta^2.*Path_loss_S_IRS.*Path_loss_IRS_R./4;
        Mu_X=Mu_X+pi.*sqrt(lambda_l_sqr)./2;
        Var_X=Var_X+lambda_l_sqr.*(16-pi^2)./4;
    end
    Mu_R=sqrt(gamma_bar(pp))*Mu_X;
    Var_R=gamma_bar(pp)*Var_X;
   
    E_R2=Var_R+Mu_R.^2;
    tot_pwr=gamma_bar(pp)*(Path_loss_R_D)*E_R2/(gamma_bar(pp)*E_R2+1);
    %----------------------------------------------------------------------
    % EH analytical--------------------------------------------------------
    Eh=eff*alpha*Tc*Ptx*tot_pwr + eff*(1-alpha)*Psi_ps*Tc*Ptx*tot_pwr ;
    Eh_Analytical=[Eh_Analytical Eh];

%% output simulation progress-----------------------------------------------    
% display 
     disp(['Simulation:solve for ' num2str(pp) ' out of ' num2str(length(P))]);
    
end
%% plotting the results
P_db=10*log10(P);
figure(3)
plot(P,Eh_sim_exct,'--o','color',[1 0 0],'LineWidth',2,'MarkerSize',9);hold on;
plot(P,Eh_Analytical,'--','color',[0 0 1],'LineWidth',2,'MarkerSize',9);hold on;
end
% figure(3)
% plot(P_db,rate_sim_exct,'--o','color',[1 0 1],'LineWidth',2,'MarkerSize',9);hold on;
% plot(gamma_bar_db,rate_sim_ub,':s','color',[0 0 1],'LineWidth',2,'MarkerSize',10);hold on;
% plot(gamma_bar_db,rate_sim_min,':s','color',[0 1 1],'LineWidth',2,'MarkerSize',10);hold on;
% plot(gamma_bar_db,rate_Analytical_ub_int,'-','color',[1 0.5 1],'LineWidth',2);hold on;
% plot(gamma_bar_db,rate_Analytical_ub_apx,'--','color',[0 0.5 0],'LineWidth',2);