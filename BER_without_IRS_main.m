%==========================================================================

%==========================================================================
% clear all;
% clc;
disp('_________________IRS and Relay Combination________________________');
% taking input for the number of iterations
NS=input('input number of samples for estimation (1000s) : ')*1000;
%--------------------------------------------------------------------------
%==========================================================================
%% system parameters
N=1;    % No. of IRSs
% L=64;   % No. of elements
M=1;    % No. of BS antennas
K=1;    % No. of users 
Nr=1;   % No. of relays
sigma_N=1;  % noise power
eta=1;  % reflection coefficient
omega=1;
vartheta=2;

%% variables
gamma_bar=logspace(-3,3,20); 

%% path-loss - correlation,large scale fading with shadowing
[Path_loss_S_IRS,Path_loss_IRS_R,Path_loss_R_D,Path_loss_S_R]=PathLoss_model(M,K,N,Nr);


%% simulation
ber_sim=[];
ber_sim_sr=[];
for pp=1:length(gamma_bar)
    gamma_out_count=0;
    gamma_out_count_sr=0;
    for ns=1:NS
% actual channels----------------------------------------------------------
        alpha_h=abs(sqrt(Path_loss_S_IRS).*(randn(M,L)+1j*randn(M,L))/sqrt(2));
        alpha_g=abs(sqrt(Path_loss_IRS_R).*(randn(Nr,L)+1j*randn(Nr,L))/sqrt(2));
        alpha_f=abs(sqrt(Path_loss_R_D).*(randn(K,Nr)+1j*randn(K,Nr))/sqrt(2));
        alpha_sr=abs(sqrt(Path_loss_S_R).*(randn(K,Nr)+1j*randn(K,Nr))/sqrt(2));
% SNRs calculation---------------------------------------------------------
        gamma_R=gamma_bar(pp).*sum(eta.*alpha_h.*alpha_g,2).^2;
        gamma_D=gamma_bar(pp).*alpha_f.^2;
        gamma_star=(gamma_R.*gamma_D)./(gamma_R+gamma_D+1);
        gamma_out_count=gamma_out_count+omega*qfunc(sqrt(vartheta.*gamma_star));
        
        gamma_R_sr=gamma_bar(pp).*alpha_sr.^2;
        gamma_D_sr=gamma_bar(pp).*alpha_f.^2;
        gamma_star_sr=(gamma_R_sr.*gamma_D_sr)./(gamma_R_sr+gamma_D_sr+1);
        gamma_out_count_sr=gamma_out_count_sr+omega*qfunc(sqrt(vartheta.*gamma_star_sr));
    end
    ber_inst=gamma_out_count/NS;
    ber_sim=[ber_sim ber_inst];
    
    ber_inst_sr=gamma_out_count_sr/NS;
    ber_sim_sr=[ber_sim_sr ber_inst_sr];
    
% output simulation progress-----------------------------------------------
    disp(['Simulation:solve for ' num2str(pp) ' out of ' num2str(length(gamma_bar))]);
end

%% plotting the results
gamma_bar_db=10*log10(gamma_bar);
figure(1)
semilogy(gamma_bar_db,ber_sim,'--o','color',[1 0 0],'LineWidth',2,'MarkerSize',9);hold on;
semilogy(gamma_bar_db,ber_sim_sr,'-','color',[0 0 1],'LineWidth',2);hold on;
xlim([-10 30])

