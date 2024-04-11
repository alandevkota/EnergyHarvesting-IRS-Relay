function [Path_loss_S_IRS,Path_loss_IRS_R,Path_loss_R_D,Path_loss_S_R]=PathLoss_model(Num_BS_antennas,Num_US,Num_IRS,Num_R)
% Num_BS_antennas=1;   % No. of BS antennas 
% Num_US=1;     % No. of total users
% Num_IRS=1;      % No. of IRS elements
% Num_R=1;

Distance_scale = 1000;       % scale the distance

%% creating the area
% area cordinates
a = 1;      % cordinate scale
x_area = [2*a 1*0 0 2*a 2*a];   % x cordinates of the area
y_area = [1*a 1*a 0 1*0 1*a];   % y cordinates of the area

f1=figure(15);
% ploting the boundary
plot(x_area*Distance_scale, y_area*Distance_scale, 'k-', 'LineWidth', 2);
% axis square;
grid on
hold on;

%% BS placenemt in the area 
x_BS = 0.2.*ones(1, Num_BS_antennas);
y_BS = 0.1.*ones(1, Num_BS_antennas);
% ploting the APs
plot((x_BS)*Distance_scale, (y_BS)*Distance_scale, 'rp', 'LineWidth', 4);
hold on;

%% IRS placenemt in the area 
x_IRS = 0.6;
y_IRS = 0.5;
% ploting the APs
plot((x_IRS)*Distance_scale, (y_IRS)*Distance_scale, 'ks', 'LineWidth', 10);
hold on;

%% R placenemt in the area 
x_R = 1.4;
y_R = 0.6;
% ploting the APs
plot((x_R)*Distance_scale, (y_R)*Distance_scale, 'bs', 'LineWidth', 10);
hold on;

%% user pacement
x_US = 1.8;
y_US = 0.9;
plot(x_US*Distance_scale, y_US*Distance_scale, 'go', 'LineWidth', 2);


%% path loss model
% path gain calculation
% Constant term in pathloss model (assuming distances in meters)
%% parameters
d0 = 1;     % reference distance  
sigma_sf = 8;    % shadow fading standared diviation in dB
G_tx = 0;      % trasmit antenna gain
G_rx = 0;       % receiver antenna gain
K0 = -41.96;       % pathloss at reference distance (5m = -55.96 / 1m = -41.96 )
bandwidth = 20e6;       % channel bandwidth
noiseFigure = 10;       % noise figure at BS (in dB)
noiseVariancedBm = -174 + 10*log10(bandwidth) + noiseFigure;        % compute total noise power
eta = 3.4;       %path loss exponent

%% path gain calculation BS to IRS
Path_loss_S_IRS = [];
for ii = 1:Num_IRS
    % selecting ii_th user
    x_IRS_ii = x_IRS(ii);     % assign ii_th user-x cordinates
    y_IRS_ii = y_IRS(ii);     % assign ii_th user-y cordinates
    
    zeta_H_ll = [];
    jj = 0;
    while jj < Num_BS_antennas
        jj = jj + 1;
        % selecting jj_th AP
        x_US_ii = x_BS(jj);     % assign jj_th AP-x cordinates
        y_US_ii = y_BS(jj);     % assign jj_th AP-y cordinates
        d_AP_US = sqrt((x_US_ii - x_IRS_ii)^2 + (y_US_ii - y_IRS_ii)^2)*Distance_scale;       % distance between ii_th user and the jj_th AP
        X_sf = sigma_sf*randn(1,1);      % shadow fading in dB for the ii_th user 
        P_L_dB = -K0 + eta*10*log10(d_AP_US/d0) + X_sf + noiseVariancedBm - G_tx - G_rx;      % pathloss in dB for the ii_th user 
        P_L_linear = 10^(P_L_dB/10);       % pathloss in linear scale for the ii_th user 
        P_G_linear = 1/P_L_linear;      % calcutate path gain
        
        if P_G_linear >= 0.1 || P_G_linear <= 0.01
            jj = jj -1;
        else
            zeta_H_ll = [zeta_H_ll P_G_linear];        % storing path gain per user
        end
    end
    
    Path_loss_S_IRS = [Path_loss_S_IRS  zeta_H_ll]; 
end

%% path gain calculation IRS to R
Path_loss_IRS_R = [];
for ii = 1:Num_IRS
    % selecting ii_th user
    x_IRS_ii = x_IRS(ii);     % assign ii_th user-x cordinates
    y_IRS_ii = y_IRS(ii);     % assign ii_th user-y cordinates
    
    zeta_F_ll = [];
    jj = 0;
    while jj < Num_R
        jj = jj + 1;
        % selecting jj_th AP
        x_US_ii = x_R(jj);     % assign jj_th AP-x cordinates
        y_US_ii = y_R(jj);     % assign jj_th AP-y cordinates
        d_BS_IRS = sqrt((x_US_ii - x_IRS_ii)^2 + (y_US_ii - y_IRS_ii)^2)*Distance_scale;       % distance between ii_th user and the jj_th AP
        X_sf = sigma_sf*randn(1,1);      % shadow fading in dB for the ii_th user 
        P_L_dB = -K0 + eta*10*log10(d_BS_IRS/d0) + X_sf + noiseVariancedBm - G_tx - G_rx;      % pathloss in dB for the ii_th user 
        P_L_linear = 10^(P_L_dB/10);       % pathloss in linear scale for the ii_th user 
        P_G_linear = 1/P_L_linear;      % calcutate path gain
        
        if P_G_linear >= 0.1 || P_G_linear <= 0.01
            jj = jj -1;
        else
            zeta_F_ll = [zeta_F_ll P_G_linear];        % storing path gain per user
        end
    end
    
    Path_loss_IRS_R = [Path_loss_IRS_R zeta_F_ll]; 
end

%% path gain calculation R to users
Path_loss_R_D = [];
for ii = 1:Num_US
    % selecting ii_th user
    x_US_ii = x_US(ii);     % assign ii_th user-x cordinates
    y_US_ii = y_US(ii);     % assign ii_th user-y cordinates
    
    zeta_G_ll = [];
    jj = 0;
    while jj < Num_R
        jj = jj + 1;
        % selecting jj_th AP
        x_BS_jj = x_R(jj);     % assign jj_th AP-x cordinates
        y_BS_jj = y_R(jj);     % assign jj_th AP-y cordinates
        d_IRS_US = sqrt((x_BS_jj - x_US_ii)^2 + (y_BS_jj - y_US_ii)^2)*Distance_scale;       % distance between ii_th user and the jj_th AP
        X_sf = sigma_sf*randn(1,1);      % shadow fading in dB for the ii_th user 
        P_L_dB = -K0 + eta*10*log10(d_IRS_US/d0) + X_sf + noiseVariancedBm - G_tx - G_rx;      % pathloss in dB for the ii_th user 
        P_L_linear = 10^(P_L_dB/10);       % pathloss in linear scale for the ii_th user 
        P_G_linear = 1/P_L_linear;      % calcutate path gain
        
        if P_G_linear >= 0.1 || P_G_linear <= 0.01
            jj = jj -1;
        else
            zeta_G_ll = [zeta_G_ll P_G_linear];        % storing path gain per user
        end
    end
    
    Path_loss_R_D = [Path_loss_R_D  zeta_G_ll]; 
end

%% path gain calculation S to R
Path_loss_S_R = [];
for ii = 1:Num_R
    % selecting ii_th user
    x_R_ii = x_R(ii);     % assign ii_th user-x cordinates
    y_R_ii = y_R(ii);     % assign ii_th user-y cordinates
    
    zeta_R_ll = [];
    jj = 0;
    while jj < Num_BS_antennas
        jj = jj + 1;
        % selecting jj_th AP
        x_US_ii = x_BS(jj);     % assign jj_th AP-x cordinates
        y_US_ii = y_BS(jj);     % assign jj_th AP-y cordinates
        d_AP_US = sqrt((x_US_ii - x_R_ii)^2 + (y_US_ii - y_R_ii)^2)*Distance_scale;       % distance between ii_th user and the jj_th AP
        X_sf = sigma_sf*randn(1,1);      % shadow fading in dB for the ii_th user 
        P_L_dB = -K0 + eta*10*log10(d_AP_US/d0) + X_sf + noiseVariancedBm - G_tx - G_rx;      % pathloss in dB for the ii_th user 
        P_L_linear = 10^(P_L_dB/10);       % pathloss in linear scale for the ii_th user 
        P_G_linear = 1/P_L_linear;      % calcutate path gain
        
        if P_G_linear >= 0.1 || P_G_linear <= 0.003
            jj = jj -1;
        else
            zeta_R_ll = [zeta_R_ll P_G_linear];        % storing path gain per user
        end
    end
    
    Path_loss_S_R = [Path_loss_S_R  zeta_R_ll]; 
end

close(f1)
