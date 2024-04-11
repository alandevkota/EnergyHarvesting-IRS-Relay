%==========================================================================
% Date: 12/21/2021
%==========================================================================
% close all;
clear all;
disp('_________________IRS and Relay Combination________________________');
NS=input('Number of samples in 1000 ')*1000; % taking inputs
%==========================================================================
alpha_array = [0.5]; %[0 0.5 0.5]; TS factor array
Psi_ps_array = [0.5]; %1 %[0.5 0 0.5]; %PS factor array
Tc = 1; %coherence time

eff = 0.9;%0.5; % Energy harvesting efficiency

bit_array = [2]; %[1 2 4];

delta= 9.079; omega= 2.9 ;psi=47.083*(10^(-3)); %psi =a, d=b=omega, delta = K
% delta= 20000; omega= 2.9 ;psi=6400;

nu=1/(1+exp(psi*omega));
PsiEH=@(p)  max(0,(delta/(1-nu))*((1./(1+exp(-psi*(p-omega))))-nu));


%% system parameters
N=1;    % No. of IRSs
L=256;   % No. of elements
M=1;    % No. of BS antennas
K=1;    % No. of users 
Nr=1;   % No. of relays
sigma_N=1;  % noise power
eta=0.9;  % reflection coefficient
%%
%**************************************************************************
% JUST ASSUMED THAT NOSE VARIANCE ARE EQUAL FOR THE RELAY AND THE DESTINATION
%**************************************************************************
% OTHERWISE WRITE DIFFERENT CODE BELOW FOR TWO NOISE VARIANCES
%--------------------------------------------------------------------------
bandwidth = 20e6;       % channel bandwidth
noiseFigure = 10;       % noise figure at BS (in dB)
noiseVariancedBm = -174 + 10*log10(bandwidth) + noiseFigure;
% noiseVariance = 10^(0.1*noiseVariancedBm);
noiseVariance = 0.1;
%--------------------------------------------------------------------------
%%
for qq = 1:length(bit_array)
    bit=bit_array(qq);
    tau=pi/(2^bit);
    
for kk=1:length(alpha_array)
    %% variables
    alpha = alpha_array(kk); %TS factor
    Psi_ps = Psi_ps_array(kk); %power scaling factor
%     P=logspace(-1,3,20); 
    P=linspace(0,70000,30); %transmitted  power
%     gamma_bar=P; % assuming noise variance equal to 1
    gamma_bar=P/noiseVariance;
%% path-loss - correlation,large scale fading with shadowing
%[Path_loss_S_IRS,Path_loss_IRS_R,Path_loss_R_D]=PathLoss_model(M,K,N,Nr);
Path_loss_S_IRS = 0.005;
Path_loss_IRS_R = 0.004;
Path_loss_R_D = 0.003;
%% simulation
Eh_sim_exct=[]; %--Energy harvest vector (like rate vector)
Eh_sim_exct_q=[]; %--Energy harvest vector with quantization
Eh_Analytical=[]; %--Energy harvest vector analytical
Eh_Analytical_q=[]; %--Energy harvest vector analytical with quantization
Ehnl_Analytical=[]; %--Non linear Energy harvest vector analytical
Ehnl_Analytical_q=[]; %--Non linear Energy harvest vector analytical with quantization

rate_sim_exct=[];
rate_sim_ub=[];
rate_sim_min=[]; 

rate_Analytical_ub_int=[];
rate_Analytical_ub_apx=[];
rate_final=[];

%==========================================================================
% GammaTH threshold
Rate_vec=[];

P1_vec=[];
P2_vec=[];

quantized_Rate_vec=[];

quantized_P1_vec=[];
quantized_P2_vec=[];
%==========================================================================

for pp=1:length(P)
    Ptx=P(pp);
    
    Eh_count_exct=0; %--EH count
    Eh_count_exct_q=0; %--EH quantized count
    
    r_gamma_out_count_exct=0;
    gamma_out_count_ub=0;
    gamma_star_count_min=0;
    
%==========================================================================
Rate_count=0;
Eh_count=0;
power1_count=0;
power2_count=0;

quantized_Rate_count=0;
quantized_power1_count=0;
quantized_power2_count=0;
%==========================================================================
%% Added by Dulaj
    E_WR = 0;
    E_WR2 = 0;
    E_WI2 = 0;
    
    for ll = 1:L % I have used zeta here
        E_WR = E_WR +  sqrt(Path_loss_S_IRS.*Path_loss_IRS_R);
        E_WR2 = E_WR2 +  Path_loss_S_IRS.*Path_loss_IRS_R*(0.5 + sin(2*tau)/(4*tau)) - pi*pi*Path_loss_S_IRS.*Path_loss_IRS_R*sin(tau)*sin(tau)/(16*tau*tau); % this is not E[(W_R)^2], this is variance
        E_WI2 = E_WI2 +  Path_loss_S_IRS.*Path_loss_IRS_R*(0.5 - sin(2*tau)/(4*tau));
    end
    
    E_WR = E_WR * eta*pi*sin(tau)/(4*tau);
    E_WR2 = E_WR2*eta*eta + (E_WR^2); % this is E[(W_R)^2]
    E_WI2 = E_WI2*eta*eta;
    G_opt_q = sqrt(gamma_bar(pp)./(gamma_bar(pp).*(E_WR2+E_WI2).^2+1));
    
    %%
    
    for ns=1:NS %number of realization = NS
% actual channels----------------------------------------------------------
        alpha_h=abs(sqrt(Path_loss_S_IRS).*(randn(M,L)+1j*randn(M,L))/sqrt(2));
        alpha_g=abs(sqrt(Path_loss_IRS_R).*(randn(Nr,L)+1j*randn(Nr,L))/sqrt(2));
        alpha_f=abs(sqrt(Path_loss_R_D).*(randn(K,Nr)+1j*randn(K,Nr))/sqrt(2));
        
% SNRs calculation---------------------------------------------------------
%         gamma_R=gamma_bar(pp).*sum(eta.*alpha_h.*alpha_g,2).^2;
%         gamma_D=gamma_bar(pp).*alpha_f.^2;
%         gamma_star=(gamma_R.*gamma_D)./(gamma_R+gamma_D+1);
%         r_gamma_out_count_exct=r_gamma_out_count_exct+log2(1+gamma_star);
%         gamma_out_count_ub=gamma_out_count_ub+gamma_star;
%         gamma_star_count_min=gamma_star_count_min+min(gamma_R,gamma_D);        
        
% Assuming perfect phase cancelation---------------------------------------

%          amp1=sum(eta.*alpha_h.*alpha_g,2); % breaking the expression into two parts and 1st part has amplitude of amp1
%          amp2=sqrt(Psi_ps)*amp1; % second amplitude has powerscaling factor in it
         
         %?????????????????????????????????????????????????????????????????
         % Optimum relay gain assuming noise power equal are equal (NEEDS CHANGE IF WE USE DIFFERENT NOISE POWER FOR USER AND RELAY SIDE)
%          G_opt = sqrt(gamma_bar(pp)./(gamma_bar(pp).*(amp1).^2+1)); 
         %?????????????????????????????????????????????????????????????????
         
%          power=eff*Ptx*alpha*Tc*(G_opt.*alpha_f.*amp1)^2+eff*Ptx*(1-alpha)*Tc*(G_opt.*alpha_f.*amp2)^2;
         
%          Eh_count_exct=Eh_count_exct+power;
 %==========================================================================        
%         Rate_count=Rate_count+power;
%         power1=Ptx*(sum(amp1).*G_opt.*alpha_f)^2 ;
%         power2=Ptx*(sum(amp2).*G_opt.*alpha_f)^2;
        
%         power1_count = power1_count+power1;
%         power2_count = power2_count+power2;
 %==========================================================================
% quantization
        q_error=2*tau*rand(1,L)-tau;
        amp1_q=sum(eta.*alpha_h.*alpha_g.*exp(1j*q_error),2);
%         amp2_q=sqrt(Psi_ps)*amp1_q;
        
        %%%%%%%%%% by Dulaj %%%%%%%%%%%%% G_opt_q = sqrt(gamma_bar(pp)./(gamma_bar(pp).*(amp1_q).^2+1)); 
        
%         power_q=eff*Ptx*alpha*Tc*(abs(G_opt_q.*alpha_f).*amp1_q)^2+eff*Ptx*(1-alpha)*Tc*(G_opt_q.*alpha_f.*amp2_q)^2;
%         Eh_count_exct_q=Eh_count_exct_q+power_q;
%     end
    
    
    
    %========================================================================== 
     
    quantized_amplitudes1 = sum(eta.*alpha_h.*alpha_g.*exp(1j*q_error),2);
    quantized_amplitudes2 = sqrt(Psi_ps)*amp1_q;
%     G_opt_q = sqrt(gamma_bar(pp)./(gamma_bar(pp).*(quantized_amplitudes1).^2+1));
    
    %% Added by Dulaj
%     E_WR = 0;
%     E_WR2 = 0;
%     E_WI2 = 0;
%     
%     for ll = 1:L % I have used zeta here
%         E_WR = E_WR +  sqrt(Path_loss_S_IRS.*Path_loss_IRS_R);
%         E_WR2 = E_WR2 +  Path_loss_S_IRS.*Path_loss_IRS_R*(0.5 + sin(2*tau)/(4*tau)) - pi*pi*Path_loss_S_IRS.*Path_loss_IRS_R*sin(tau)*sin(tau)/(16*tau*tau); % this is not E[(W_R)^2], this is variance
%         E_WI2 = E_WI2 +  Path_loss_S_IRS.*Path_loss_IRS_R*(0.5 - sin(2*tau)/(4*tau));
%     end
%     
%     E_WR = E_WR * eta*pi*sin(tau)/(4*tau);
%     E_WR2 = E_WR2*eta*eta + (E_WR^2); % this is E[(W_R)^2]
%     E_WI2 = E_WI2*eta*eta;
%     G_opt_q = sqrt(gamma_bar(pp)./(gamma_bar(pp).*(E_WR2+E_WI2).^2+1));
    
    %%
    
    quantized_power = eff*Ptx*alpha*Tc*(abs(G_opt_q.*quantized_amplitudes1*alpha_f))^2+eff*Ptx*(1-alpha)*Tc*(abs(G_opt_q.*alpha_f.*quantized_amplitudes2))^2;
    
    quantized_Rate_count=quantized_Rate_count+ quantized_power;
    quantized_power1 = Ptx*(sum(quantized_amplitudes1).*G_opt_q.*alpha_f)^2;
    quantized_power2 = Ptx*(sum(quantized_amplitudes2).*G_opt_q.*alpha_f)^2;
    
    quantized_power1_count = quantized_power1_count+quantized_power1;
    quantized_power2_count = quantized_power2_count+quantized_power2;
    %========================================================================== 
end
    
    
    
    
%==========================================================================  
%     Rate_vec=[Rate_vec (Rate_count/NS)];
%     P1_vec=[P1_vec (power1_count/NS)];
%     P2_vec=[P2_vec (power2_count/NS)];
    
    quantized_Rate_vec=[quantized_Rate_vec (quantized_Rate_count/NS)];
    quantized_P1_vec=[quantized_P1_vec (quantized_power1_count/NS)];
    quantized_P2_vec=[quantized_P2_vec (quantized_power2_count/NS)];
%==========================================================================  

%     average_snr_inst_ub=gamma_out_count_ub/NS;
%     average_snr_min=gamma_star_count_min/NS;    
    
% EH calculation for continuous and also for quantized --------------------

%     Eh_sim_exct=[Eh_sim_exct (Eh_count_exct/NS)];
%     Eh_sim_exct_q=[Eh_sim_exct_q (Eh_count_exct_q/NS)];
    
% rate calculation---------------------------------------------------------
%     rate_sim_exct=[rate_sim_exct r_gamma_out_count_exct/NS];
%     rate_sim_ub=[rate_sim_ub log2(1+average_snr_inst_ub)];
%     rate_sim_min=[rate_sim_min log2(1+average_snr_min)];
    
        
%% EH analysis
%     Mu_X=0;
%     Var_X=0;
%     for ll=1:L % I have used xi here instead of zeta, where zeta=xi/2
%         lambda_l_sqr=eta^2.*Path_loss_S_IRS.*Path_loss_IRS_R./4;
%         Mu_X=Mu_X+pi.*sqrt(lambda_l_sqr)./2;
%         Var_X=Var_X+lambda_l_sqr.*(16-pi^2)./4;
%     end
%     Mu_R=sqrt(gamma_bar(pp))*Mu_X;
%     Var_R=gamma_bar(pp)*Var_X;
    
%     sigma_D=gamma_bar(pp)*Path_loss_R_D;
%     psi=1/(qfunc(-Mu_R/sqrt(Var_R)));

%     E_R2=Var_R+Mu_R.^2;
%     tot_pwr=gamma_bar(pp)*(Path_loss_R_D)*E_R2/(gamma_bar(pp)*E_R2+1);
   %==========================================================================   
    
%     tmp1 = PsiEH(P1_vec);
%     tmp2 = PsiEH(P2_vec);
    
    quantized_tmp1 = PsiEH(quantized_P1_vec);
    quantized_tmp2 = PsiEH(quantized_P2_vec);
    %==========================================================================  
    %----------------------------------------------------------------------
    % EH analytical--------------------------------------------------------
    
    
    % Rate analytical--------------------------------------------------------
%     Eh=eff*alpha*Tc*Ptx*tot_pwr + eff*(1-alpha)*Psi_ps*Tc*Ptx*tot_pwr ;
%     Eh_Analytical=[Eh_Analytical Eh];
    
%     Eh_vec = alpha*Tc*tmp1 + (1-alpha)*Tc*tmp2;
    quantized_Eh_vec = alpha*Tc*quantized_tmp1 + (1-alpha)*Tc*quantized_tmp2;
    
    %==========================================================================  
%     Eh1=Ptx*tot_pwr ;
%     Eh2=Psi_ps*Ptx*tot_pwr ;
    
    
%     Eh_nl = alpha*Tc*PsiEH(Eh1) + (1-alpha)*Tc*PsiEH(Eh2);
%     Ehnl_Analytical = [Ehnl_Analytical Eh_nl];
%     Eh_Analytical_nl=[Eh_Analytical_nl Eh_nl];
    %==========================================================================  
     %==========================================================================     
% % %     % integral calculation-----------------------------------------------------
% % %     fun_int=@(u) (sqrt(Var_R).*u+Mu_R).*qfunc(u).*exp(-(sqrt(Var_R).*u+Mu_R).^2./sigma_D);
% % %     I_value_int=integral(fun_int,-Mu_R/sqrt(Var_R),Inf);
% % %     average_snr_ana_ub_int=2*sqrt(Var_R)*psi*I_value_int;
% % %     rate_Analytical_ub_int=[rate_Analytical_ub_int log2(1+average_snr_ana_ub_int)]; %try
    % analytical calculation---------------------------------------------------
%     average_snr_ana_ub_apx=1/2*sigma_D*psi*(1-erf(-Mu_R/(sqrt(2*Var_R))));
%     rate_Analytical_ub_apx=[rate_Analytical_ub_apx log2(1+average_snr_ana_ub_apx)];  
    
%% EH quantized analysis
%     E_WR = 0;
%     E_WR2 = 0;
%     E_WI2 = 0;
%     
%     for ll = 1:L % I have used zeta here
%         E_WR = E_WR +  sqrt(Path_loss_S_IRS.*Path_loss_IRS_R);
%         E_WR2 = E_WR2 +  Path_loss_S_IRS.*Path_loss_IRS_R*(0.5 + sin(2*tau)/(4*tau)) - pi*pi*Path_loss_S_IRS.*Path_loss_IRS_R*sin(tau)*sin(tau)/(16*tau*tau); % this is not E[(W_R)^2], this is variance
%         E_WI2 = E_WI2 +  Path_loss_S_IRS.*Path_loss_IRS_R*(0.5 - sin(2*tau)/(4*tau));
%     end
%     
%     E_WR = E_WR * eta*pi*sin(tau)/(4*tau);
%     E_WR2 = E_WR2*eta*eta + (E_WR^2); % this is E[(W_R)^2]
%     E_WI2 = E_WI2*eta*eta;
    
    tot_pwr_q=(G_opt_q^2)*(Path_loss_R_D)*(E_WR2+E_WI2);
    % EH quantized analytical----------------------------------------------
    
    Eh_q=(eff*alpha*Tc*Ptx+ eff*(1-alpha)*Psi_ps*Tc*Ptx)*tot_pwr_q;
    Eh_Analytical_q=[Eh_Analytical_q Eh_q];
    
    %==========================================================================   

    quantized_Eh1 = Ptx*G_opt_q^2*(Path_loss_R_D)*(E_WR2+E_WI2);
    quantized_Eh2 = Psi_ps*quantized_Eh1;
    
    quantized_temp1 = PsiEH(quantized_Eh1);
    quantized_temp2 = PsiEH(quantized_Eh2);
    
    quantized_Eh_nl = alpha*Tc*quantized_temp1 + (1-alpha)*Tc*quantized_temp2;
    Ehnl_Analytical_q=[Ehnl_Analytical_q quantized_Eh_nl];
    
%     quantized_ratio_linear = 100.*quantized_Rate_vec./Rate_vec;
%     quantized_ratio_non_linear = 100.*quantized_Eh_vec./Eh_vec;
    %========================================================================== 
    
%     SNR_avg = Ptx*G_opt_q^2*Path_loss_R_D*(E_WR2+E_WI2)/(G_opt_q^2*Path_loss_R_D*noiseVariance+noiseVariance);
    
%     SNR_test = SNR_avg/NS;
%     rate_final=[rate_final log2(1+SNR_avg)]; 
    %--------------quantization ratio--------------------------------------
%     quantizn_ratio_linear = 100.*Eh_sim_exct_q/Eh_sim_exct;
    
    

%% output simulation progress-----------------------------------------------    
% display 
     %disp(['Simulation:solve for ' num2str(pp) ' out of ' num2str(length(P))]);
    
end
%% plotting the results
% P_db=10*log10(P);
gamma_bar_db=10*log10(gamma_bar);
gammabardB=10*log10(P);

% figure(3)
% plot(P,Eh_sim_exct,'--o','color',[1 0 0],'LineWidth',2,'MarkerSize',9);hold on;
% plot(P,Eh_Analytical,'-','color',[0 0 1],'LineWidth',2,'MarkerSize',9);hold on;
% plot(P,Eh_sim_exct_q,'--o','color',[1 0 1],'LineWidth',2,'MarkerSize',9);hold on;
% plot(P,Eh_Analytical_q,'--s','color',[0 1 1],'LineWidth',2,'MarkerSize',9);hold on;

% figure(4)
% plot(gamma_bar_db,rate_sim_exct,'--o','color',[1 0 1],'LineWidth',2,'MarkerSize',9);hold on;
% plot(gamma_bar_db,rate_sim_ub,':s','color',[0 0 1],'LineWidth',2,'MarkerSize',10);hold on;
% plot(gamma_bar_db,rate_sim_min,':s','color',[0 1 1],'LineWidth',2,'MarkerSize',10);hold on;
% plot(gamma_bar_db,rate_Analytical_ub_int,'-','color',[1 0.5 1],'LineWidth',2);hold on;
% plot(gamma_bar_db,rate_Analytical_ub_apx,'--','color',[0 0.5 0],'LineWidth',2);
% plot(gamma_bar_db,rate_final,'--','color',[0 0.5 0],'LineWidth',2);

% figure(5)
% plot(gamma_bar,rate_sim_exct,'--o','color',[1 0 1],'LineWidth',2,'MarkerSize',9);hold on;
% plot(gamma_bar,rate_sim_ub,':s','color',[0 0 1],'LineWidth',2,'MarkerSize',10);hold on;
% plot(gamma_bar,rate_sim_min,':s','color',[0 1 1],'LineWidth',2,'MarkerSize',10);hold on;
% plot(gamma_bar,rate_Analytical_ub_int,'-','color',[1 0.5 1],'LineWidth',2);hold on;
% plot(gamma_bar,rate_Analytical_ub_apx,'--','color',[0 0.5 0],'LineWidth',2);
% plot(gamma_bar,rate_final,'--','color',[0 0.5 0],'LineWidth',2);

% figure(6)
% plot(gamma_bar_db,10.^(0.1*rate_sim_exct),'--o','color',[1 0 1],'LineWidth',2,'MarkerSize',9);hold on;
% plot(gamma_bar_db,10.^(0.1*rate_sim_ub),':s','color',[0 0 1],'LineWidth',2,'MarkerSize',10);hold on;
% plot(gamma_bar_db,10.^(0.1*rate_sim_min),':s','color',[0 1 1],'LineWidth',2,'MarkerSize',10);hold on;
% plot(gamma_bar_db,10.^(0.1*rate_Analytical_ub_int),'-','color',[1 0.5 1],'LineWidth',2);hold on;
% plot(gamma_bar_db,10.^(0.1*rate_Analytical_ub_apx),'--','color',[0 0.5 0],'LineWidth',2);
% plot(gamma_bar_db,10.^(0.1*rate_final),'--','color',[0 0.5 0],'LineWidth',2);
%========================================================================== 
    
% figure(7)
% % plot(Eh_sim_exct,rate_sim_exct,'--o','color',[1 0 0],'LineWidth',2,'MarkerSize',9);hold on;
% 
% plot(Eh_Analytical_q,rate_Analytical_ub_apx,'-','color',[0 0 1],'LineWidth',2);hold on;
% plot(Eh_nl,rate_Analytical_ub_apx,'o','color',[0 0 1],'LineWidth',2);
% plot(Eh_Analytical,rate_Analytical_ub_apx,'--','color',[0 0 1],'LineWidth',2);
% %========================================================================== 
    
% figure(8)
% 
% plot(P,Rate_vec,'o','color',[0 0 1],'LineWidth',2,'MarkerSize',8); hold on;
% plot(P,Eh_Analytical,'-','color',[1 0 0],'LineWidth',2,'MarkerSize',8); hold on;
% plot(P,Eh_vec,'o','color',[1 0 1],'LineWidth',2,'MarkerSize',8); hold on;
% plot(P,Ehnl_Analytical,'-','color',[0 0 0],'LineWidth',2,'MarkerSize',8); hold on;

% 
figure(9);
plot(P,quantized_Rate_vec,'o','color',[0 0 1],'LineWidth',2,'MarkerSize',8); hold on;
plot(P,Eh_Analytical_q,'-','color',[1 0 0],'LineWidth',2,'MarkerSize',8); hold on;
plot(P,quantized_Eh_vec,'+','color',[1 0 1],'LineWidth',2,'MarkerSize',8); hold on;
plot(P,Ehnl_Analytical_q,'-','color',[0 0 0],'LineWidth',2,'MarkerSize',8); hold on;
% plot(gammabardB,lower_bound,'--','color',[0 0 0],'LineWidth',2); hold on;

% figure(10);
% plot(Pt,quantized_ratio_non_linear,'--*','color',[1 0 0],'LineWidth',2)
% hold on;  
% plot(Pt,quantized_ratio_linear,'--*','color',[0 0 1],'LineWidth',2)
end
end