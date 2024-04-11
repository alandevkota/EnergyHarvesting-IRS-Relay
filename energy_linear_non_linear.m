%==========================================================================
% Date: 12/01/2019

%==========================================================================
clear all;
clc;
disp('______________________Rate_IRS___________________________________');
N=input('Number of samples in 1000 ')*1000; 
%--------------------------------------------------------------------------
M=256; %# of elements
betam = 0.9; %refln coeff
alpha_array = [0.5]; %[0 0.5 0.5];
theta_array = [0]; %[0.5 0 0.5];
Tc = 1;

eff = 0.8;%0.5; % Energy harvesting efficiency

delta= 9.079; omega= 2.9 ;psi=47.083*1e-3; % delta is K in paper, and nu is ohm

nu=1/(1+exp(psi*omega));

PsiEH=@(p) max(0, delta/(1-nu)*(1./(1+exp(-psi*(p-omega)))-nu));

for kk=1:length(alpha_array)
    alpha = alpha_array(kk);
    theta = theta_array(kk);


%function M_4 = Eng(M,N,betam,alpha,theta,Tc)

%% System parameters
Ntx=1; %# transmit antennas
Nrx=1; %# receive antennas
%M=256; % IRS elements
%zetah=0.1; % PL parameter
%zetag=0.1;%PL parameter
%betam=0.9; %Refection coefficient
eta = 3.5; %Pathloss exponent
frequency = 3*(10^9);
d0 = 1;
d_closest_tx = 50;
d_closest_rx = 70;

aj = 64;
bj = 0.003;
Mj = 0.02;
Oj = 1/(1 + exp(aj*bj));

PL_BS_IRS = Path_Loss_Matrix(M, Ntx, d0, eta, frequency,d_closest_tx);
PL_IRS_MS = Path_Loss_Matrix(M, Nrx, d0, eta, frequency,d_closest_rx);
PL_BS_MS  = Path_Loss_Matrix(1, 1, d0, eta, frequency,d_closest_tx);

zetah = eye(M);
zetag = eye(M);

for di=1:M
    zetah(di,di)=PL_BS_IRS(di);
    zetag(di,di)=PL_IRS_MS(di);
end
 
Pt=linspace(0,70,30);
Number_of_realizations=N;
%==========================================================================
% GammaTH threshold
Rate_vec=[];

P1_vec=[];
P2_vec=[];
 for sinr=1:length(Pt)
%% ==========================================================================
Ptx=Pt(sinr);
  

Rate_count=0;
Eh_count=0;
power1_count=0;
power2_count=0;

    for nn=1:Number_of_realizations
        %generate channel from Tx to IRS
        h_tilde=(randn(M,1)+1j*randn(M,1))/sqrt(2);
        h=sqrt(zetah)*h_tilde;
        %generate channel from IRS to Rx
        g_tilde=(randn(1,M)+1j*randn(1,M))/sqrt(2);
        g=g_tilde*sqrt(zetag);
        %generate channel from Tx to MS
        f_tilde=(randn(1,1)+1j*randn(1,1))/sqrt(2);
        f=sqrt(PL_BS_MS)*f_tilde;
        
        %Assume IRS perfectly cancel the phases 
        % therefore take only aplitudes
        amplitudes1=betam*abs(h).*abs(g.');
        amplitudes2=sqrt(theta)*betam*abs(h).*abs(g.');
        power=eff*Ptx*alpha*Tc*(sum(amplitudes1)+abs(f))^2 + eff*Ptx*(1-alpha)*Tc*(sum(amplitudes2)+sqrt(theta)*abs(f))^2;
        
        Rate_count=Rate_count+power;
        
        power1=Ptx*(sum(amplitudes1)+abs(f))^2 ;
        power2=Ptx*(sum(amplitudes2)+sqrt(theta)*abs(f))^2;
        
        power1_count = power1_count+power1;
        power2_count = power2_count+power2;

    end
    Rate_vec=[Rate_vec (Rate_count/Number_of_realizations)];
    P1_vec=[P1_vec (power1_count/Number_of_realizations)];
    P2_vec=[P2_vec (power2_count/Number_of_realizations)];
    
    %% display 
    B = sprintf('solve for %d from %d' ,sinr, length(Pt));
    disp (B) 
    
 end
%% ==========================================================================

%tmp1 = Mj./(1+exp(-aj.*(P1_vec-bj)));
% tmp2 = Mj./(1+exp(-aj.*(P2_vec-bj)));
 
 tmp1 = PsiEH(P1_vec);
  tmp2 = PsiEH(P2_vec);
  
 %Eh_vec = (tmp1-Mj*Oj)*alpha*Tc/(1-Oj) + (tmp2-Mj*Oj)*(1-alpha)*Tc/(1-Oj);

 Eh_vec = alpha*Tc*tmp1 + (1-alpha)*Tc*tmp2;
 
%% CLT approximation
sigma=(diag(zetah).'/2).*(diag(zetag).'/2)*betam^2;
Mean_CLT=sum(pi*sqrt(sigma)/2);
Var_CLT=sum(sigma*(4-(pi/2)^2));
E_Y2 = Var_CLT + Mean_CLT.^2;
tot_pwr = E_Y2 + sqrt(PL_BS_MS*pi)*Mean_CLT+PL_BS_MS;
%--------------------------------------------------------------------------
% Rate analytical
Eh=eff*alpha*Tc*Pt*tot_pwr + eff*(1-alpha)*theta*Tc*Pt*tot_pwr ;

Eh1=Pt*tot_pwr ;
Eh2=theta*Pt*tot_pwr ;

%temp1 = Mj./(1+exp(-aj*(Eh1-bj)));
%temp2 = Mj./(1+exp(-aj*(Eh2-bj)));

temp1 = PsiEH(Eh1);
  temp2 = PsiEH(Eh2);

%Eh_nl = alpha*Tc*(temp1-Mj*Oj)/(1-Oj) + (1-alpha)*Tc*(temp2-Mj*Oj)/(1-Oj);

Eh_nl = alpha*Tc*temp1 + (1-alpha)*Tc*temp2;
%==========================================================================


%% plot results
figure(7);
gammabardB=10*log10(Pt);

plot(Pt,Rate_vec,'o','color',[0 0 1],'LineWidth',2); hold on;
plot(Pt,Eh,'--','color',[1 0 0],'LineWidth',2); hold on;
%plot(Pt,Eh_vec,'o','color',[0 0 1],'LineWidth',2); hold on;
%plot(Pt,Eh_nl,'-','color',[1 0 0],'LineWidth',2); hold on;

%plot(gammabardB,lower_bound,'--','color',[0 0 0],'LineWidth',2); hold on;

end

legend('SIM','Ana_UB');
ylabel('Average rate [bits/s/Hz]');
xlabel('Average transmit SNR [$\bar{\gamma}$\,dB]');

%xlim([-10 30])
%ylim([1e-6 1])
