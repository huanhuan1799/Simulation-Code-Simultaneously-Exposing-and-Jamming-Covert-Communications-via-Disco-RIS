%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program: Simulation of Detection Error Probabilities and Communication 
% Performance with and without Disco RIS (DRIS)
% 
% - Generates random transmission/silence sequences (Alice -> Bob).
% - Simulates Covert Communications for various SNR levels.
% - Compares Monte Carlo simulation results with Theoretical derivations.
%
% Reference: H. Huang, H. Zhang, Y. Cai, D. Niyato, A. L. Swindlehurst, and Z. Han, 
% Simultaneously Exposing and Jamming Covert Communications via Disco Reconfigurable Intelligent Surfaces,
%IEEE J. Sel. Areas Commun., 2025 12. DOI: 10.1109/JSAC.2025.3646955
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;  
clc;
clear all 
rng(1799)

%% ------------------------ Parameter Settings ------------------------
K = 1;                % Number of users (Bob)
N = 1;                % Number of antennas at Alice
Simulation_Num = 2^17;
result = randi([0, 1], 1, Simulation_Num);  % Random sequence: 1 (Transmission), 0 (Silence)
WW = 64;              % Width of DRIS
HH = 32;              % Height of DRIS
NN = WW * HH;         % Total number of DRIS elements (N_D)

% SNR Range (dB) and Linear scale
Pmin = -15;   % Minimum SNR (dB)
Pmax = 9;     % Maximum SNR (dB)
PdBm = Pmin : 4 : Pmax;   % Step size
P = 10.^(PdBm/10);        % Convert to linear scale (P_0)

% Initialize Detection Parameters
% MM corresponds to 'M' in PDF (Number of samples within coherence time)
MM = 5;       
% N_The corresponds to 'N' in PDF (Detection threshold parameter in Eq. (9))
N_The = 2;
PC = 1/2;     % Prior probability of Transmission (H1)
PS = 1 - PC;  % Prior probability of Silence (H0)

% Arrays to store Simulation and Theoretical results for Error Probabilities
p_MDR            = zeros(1, length(PdBm)); % Missed Detection Rate (With DRIS)
p_FAR            = zeros(1, length(PdBm)); % False Alarm Rate (With DRIS)
p_MDR_The        = zeros(1, length(PdBm)); % Theoretical MDR (With DRIS)
p_FAR_The        = zeros(1, length(PdBm)); % Theoretical FAR (With DRIS)
p_MDR_w          = zeros(1, length(PdBm)); % MDR (Without DRIS)
p_FAR_w          = zeros(1, length(PdBm)); % FAR (Without DRIS)
p_FAR_w_The      = zeros(1, length(PdBm)); % Theoretical FAR (Without DRIS)
p_MDR_w_The      = zeros(1, length(PdBm)); % Theoretical MDR (Without DRIS)

% Arrays to store SNR/SJNR Simulation and Theoretical results
SNR         = zeros(1, length(PdBm));  
SNR_lp      = zeros(1, length(PdBm));  % Low Power Active Jamming
SNR_hp      = zeros(1, length(PdBm));  % High Power Active Jamming
SNR_Theo    = zeros(1, length(PdBm));
SNR_w       = zeros(1, length(PdBm));  
SNR_w_Theo  = zeros(1, length(PdBm));

%% -------------------- DRIS Random Reflection Parameters --------------------
% Modeled based on Section IV Simulation Setup 
Amp = [0.8, 1];        % Discrete amplitudes (Omega)
aph = (0.8^2 + 1)/2;   % alpha_bar: Expectation of squared amplitude (E[|beta|^2])
                       % [Proposition 1, below Eq. (15)]
Omg = [pi/9, 7*pi/6];  % Discrete phase shifts (Psi)
 
% Noise Power (dBm to Linear)
% [Section IV, Simulation Results Setup]
noise = -170 + 10*log10(180*1e3);  
detaw = 10^(noise/10); % sigma_w^2

%% ---------------------- Main Loop: Iterate over SNR ----------------------
for pn = 1:length(PdBm)
    fprintf('Processing SNR index: %d / %d\n', pn, length(PdBm));
    result = randi([0, 1], 1, Simulation_Num);
    
    % Counters for correct/error decisions (With/Without DRIS)
    Dec_ture   = 0;   
    Dec_Err    = 0;
    Dec_w_ture = 0;   
    Dec_w_Err  = 0;
    
    Sig_power = 0;
    AJ_lowpower = 0;
    AJ_highpower =0;
    ACA_power = 0;
    N_power   = 0;

    %% ------------------ Iterate over each bit ------------------
    for m = 1:length(result)
        
        % Channel large-scale fading parameters (LoS + NLoS)
        % [Section II-B, Eq. (7) & Eq. (8)]
        eb = 10; % Rician factor epsilon
        eb22 = 1 / (1 + eb);
        eb11 = eb / (1 + eb);
        eb1 = sqrt(eb11);
        eb2 = sqrt(eb22);
        
        % Generate Channels (Alice-Willie, Alice-Bob, etc.)
        % [Eq. (5)~Eq. (8)]
        [Haw,Gw,Hrw,Hab,Hrb,Hjb,Hjw,lar, lrw, law,lab,lrb,ljb,ljw] = Covert_Communication_GenerateAJchannel(WW, HH, N, eb1, eb2);
            
        % Initialize received signals and thresholds
        y_w  = zeros(MM,1);  % Received signal with DRIS
        yw_w = zeros(MM,1);  % Received signal without DRIS
        HD_w = zeros(MM,1);
        thevar = zeros(MM,1);
        sigma1 =  zeros(MM,1);
        
        for t = 1:MM
            % Generate Random DRIS Phase Configuration phi(t)
            index = randi([0,1], 1, WW*HH) + 1;
            theta_inbar = Amp(index) .* exp(1j * Omg(index));      
            
            % Synthesize Cascaded Channel 
            % [Section II-A, Eq. (2)]
            HD_w(t) = Hrw * diag(theta_inbar) * Gw;
            
            % Calculate Variance under H1 (With DRIS)
            % [Eq. (17) Variance term]
            sigma1(t)   = P(pn)*norm(Haw+HD_w(t),'fro')^2 + detaw; 
            
            % Calculate Detection Threshold epsilon(m) (With DRIS)
            % [Proposition 2, Eq. (27)]
            thevar(t)   = sigma1(t) * log(sigma1(t)/ detaw) / (sigma1(t) * 10^(-noise/10) - 1);        
        end
        
        % Calculate Variance and Threshold (Without DRIS)
        % [Eq. (28) Fixed threshold epsilon]
        sigma1_w = P(pn)*norm(Haw,'fro')^2 + detaw; 
        thevar_w = sigma1_w * log(sigma1_w / detaw) / (sigma1_w * 10^(-noise/10) - 1);
        
        if result(m) == 1
            %% --------------- Case: Transmission (H1) ---------------
            % Generate QPSK Symbols (Power P(pn)) - Alice's signal s(m)
            % [Eq. (1), s(m) ~ CN(0, P0)]
            Txsig = sqrt(P(pn)/2) * (randn(MM,1) + 1j*randn(MM,1));
            
            % Active Jammer (AJ) Signals (Optional comparison in Section IV)
            AJlpsig = sqrt( (10.^(-7/10))/2) * (randn(MM,1) + 1j*randn(MM,1));
            AJhpsig = sqrt( (10.^(5/10))/2) * (randn(MM,1) + 1j*randn(MM,1));
            
            % Generate Noise Vector (Complex Gaussian)
            % [Eq. (1), n_w(m)]
            n_2 = sqrt(detaw/2) * (randn(MM,1) + 1j*randn(MM,1));
            
            % Compute Received Signals and Perform Detection
            Dec_Willie_Ture = 0;
            Dec_Willie_wo_Ture = 0;
            
            for t = 1:MM
                clear Hcom
                % Total Channel at Willie = Direct + DRIS
                Hcom = HD_w(t) + Haw;
                        
                %% Simulate SJNR (Signal-to-Jamming-plus-Noise Ratio)
                % [Eq. (14) Components]
                Sig_power = Sig_power + ((Hab*Txsig(t))'*(Hab*Txsig(t))); 
                % Active Channel Aging (ACA) interference power from DRIS
                ACA_power = ACA_power + ((Hrb * diag(theta_inbar) * Gw)*Txsig(t))'*((Hrb * diag(theta_inbar) * Gw)*Txsig(t));
                AJ_lowpower = AJ_lowpower + ((Hjb*AJlpsig(t))'*(Hjb*AJlpsig(t)));
                AJ_highpower = AJ_highpower + ((Hjb*AJhpsig(t))'*(Hjb*AJhpsig(t)));
                N_power =  N_power + n_2(t)'*n_2(t);
                
                % Received Signal at Willie (With DRIS)
                % [Eq. (16)]
                y_w(t) = Hcom * Txsig(t) + n_2(t);
                
                %% Detection (With DRIS)
                % [Eq. (9) Decision Region S]
                if norm(y_w(t),'fro')^2 >= thevar(t)
                    Dec_Willie_Ture = Dec_Willie_Ture + 1; % H1 is true, Detection is correct
                end
                
                % Received Signal at Willie (Without DRIS)
                yw_w(t) = Haw * Txsig(t) + n_2(t);
                
                %% Detection (Without DRIS)
                if norm(yw_w(t),'fro')^2 >= thevar_w
                    Dec_Willie_wo_Ture = Dec_Willie_wo_Ture + 1;
                end                
            end 
            
            %% %%%%%%%%%%%%%%%%%%%%% Theoretical SNR %%%%%%%%%%%%%%%%%%%%%%%%
            % [Theorem 2, Eq. (31)]
            % Note: This implements the asymptotic convergence limit
            SNR_Theo(pn) = SNR_Theo(pn) + (P(pn)*lab* 10^(-noise/10))/(P(pn)*WW*HH*aph*(lar*lrb)* 10^(-noise/10) + 1);
            SNR_w_Theo(pn) = SNR_w_Theo(pn) + (P(pn)*lab* 10^(-noise/10));
            
            % Calculate Missed Detection (MDR) for DRIS
            % If detections < N_The, Willie decides H0 (Missed Detection)
            if Dec_Willie_Ture < N_The 
               p_MDR(pn) = p_MDR(pn) + 1;
            end 
            
            % Calculate Missed Detection (MDR) for No DRIS
            if Dec_Willie_wo_Ture < N_The 
               p_MDR_w(pn) = p_MDR_w(pn) + 1;
            end
            
            %% %%%%%%%%%% Theoretical MDR Calculation (With DRIS) %%%%%%%%%%
            % [Theorem 1, Eq. (30)]
            % Summation over combinations of probabilities where detection fails
            for i = 1 : N_The
                clear comb
                comb = nchoosek(1:MM, i-1);
                
                if length(comb) == 0 %% No samples cross the threshold
                    clear comb1
                    comb1 = nchoosek(1:MM, MM);
                    Pro = 1;
                    for rr = 1:length( comb1 ) %% Iterate all samples
                        % CDF of |y_w|^2 under H1
                        % [Eq. (48)]
                        Pro = Pro*( 1-exp( -( thevar(comb1(rr) ))/(sigma1(comb1(rr))) ) );
                    end
                    p_MDR_The(pn) = p_MDR_The(pn) + Pro ;
                else
                    for tt = 1:length(comb(:,1))  
                        Pro = 1;
                        for rr = 1:length( comb(tt,:) )
                            Pro = Pro*exp( -( thevar(comb(tt,rr) ))/(sigma1(comb(tt,rr))) );
                        end
                        clear mis_ele
                        mis_ele = setdiff(1:MM, comb(tt,:));
                        if length(mis_ele) ~= 0
                            for jj = 1:length(mis_ele)
                                Pro = Pro* ( 1-exp( -( thevar( mis_ele(jj) ))/(sigma1( mis_ele(jj) )) ) );
                            end
                        end
                        p_MDR_The(pn) = p_MDR_The(pn) + Pro ;
                    end
                end
                
            end
            
           %% %%%%%%%%%% Theoretical MDR Calculation (No DRIS) %%%%%%%%%%
           % Standard derivation without time-varying channel effects
           for i = 1 : N_The
               clear comb
               comb = nchoosek(1:MM, i-1);
               
               if length(comb) == 0
                   clear comb1
                    comb1 = nchoosek(1:MM, MM);
                    Pro = 1;
                    for rr = 1:length( comb1 )
                        Pro = Pro*( 1-exp( -( thevar_w )/(sigma1_w) ) );
                    end
                    p_MDR_w_The(pn) = p_MDR_w_The(pn) + Pro ;
                else
                    for tt = 1:length(comb(:,1))
                        Pro = 1;
                        for rr = 1:length( comb(tt,:) )
                            Pro = Pro*exp( -( thevar_w)/(sigma1_w ) );
                        end
                        clear mis_ele
                        mis_ele = setdiff(1:MM, comb(tt,:));
                        if length(mis_ele) ~= 0
                            for jj = 1:length(mis_ele)
                                Pro = Pro* ( 1-exp( -( thevar_w )/(sigma1_w ) ) );
                            end
                        end
                        p_MDR_w_The(pn) = p_MDR_w_The(pn) + Pro ;
                    end
                end
            end 
            
        else 
            %% --------------- Case: Silence (H0) --------------- 
            % Received signal is noise only
            % [Eq. (26)]
            y_w = sqrt(detaw/2) * (randn(MM,1) + 1j*randn(MM,1));
            yw_w =sqrt(detaw/2) * (randn(MM,1) + 1j*randn(MM,1));
            
            Dec_Willie_False = 0;
            Dec_Willie_wo_False = 0;
            %% %%%%%%%%%%%%%%%%%%%%%%%% Detection Process %%%%%%%%%%%%%%%%%%%%%%
            for t = 1:MM
                % False Alarm Check (With DRIS)
                if norm(y_w(t),'fro')^2 >= thevar(t)
                    Dec_Willie_False = Dec_Willie_False + 1;
                end
                % False Alarm Check (No DRIS)
                if norm(yw_w(t),'fro')^2 >= thevar_w
                    Dec_Willie_wo_False = Dec_Willie_wo_False + 1;
                end        
            end
            %% %%%%%%%%%%%%%%%%%%%%%%%%% Decision %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % If detections >= N_The, Willie decides H1 (False Alarm)
            % [Eq. (9), (10)]
            if Dec_Willie_False >= N_The
                p_FAR(pn) = p_FAR(pn) + 1;
            end
            
            if Dec_Willie_wo_False >= N_The
               p_FAR_w(pn) = p_FAR_w(pn) + 1;
            end
            
            %% Theoretical FAR (With DRIS)    
            % [PDF: Theorem 1, Eq. (29)]
            % Probability that T >= N samples exceed threshold under H0
            for i = N_The : MM
                clear comb
                comb = nchoosek(1:MM, i);
                for tt = 1:length(comb(:,1))
                    Pro = 1;
                    for rr = 1:length( comb(tt,:) )
                        % exp(-epsilon/sigma_w^2)
                        % [Eq. (43) related terms]
                        Pro = Pro*exp( -thevar(comb(tt,rr))/detaw );
                    end
                    clear mis_ele
                    mis_ele = setdiff(1:MM, comb(tt,:));
                    if length(mis_ele) ~= 0
                        for jj = 1:length(mis_ele)
                            Pro = Pro* (1-exp( -thevar(mis_ele(jj))/detaw ) );
                        end
                    end
                    p_FAR_The(pn) = p_FAR_The(pn) + Pro;
                end
            end
            
            %% Theoretical FAR (No DRIS)
            for i = N_The : MM
               clear comb
               comb = nchoosek(1:MM, i);
               
               for tt = 1:length(comb(:,1))
                   Pro = 1;
                   for rr = 1:length( comb(tt,:) )
                       Pro = Pro*exp( -( thevar_w)/(detaw ) );
                   end
                   clear mis_ele
                   mis_ele = setdiff(1:MM, comb(tt,:));
                   if length(mis_ele) ~= 0
                       for jj = 1:length(mis_ele)
                           Pro = Pro* ( 1-exp( -( thevar_w )/(detaw ) ) );
                       end
                   end
                   p_FAR_w_The(pn) = p_FAR_w_The(pn) + Pro;
               end
            end
            
        end
    end
    mm = length( result(find(result == 1)) ); %% Count of Transmission events
    mn = length( result(find(result ~= 1)) ); %% Count of Silence events
    %% --------------- Normalize Error Rates ---------------
    p_MDR(pn) =  p_MDR(pn) / length(result);
    p_MDR_w(pn) = p_MDR_w(pn) / length(result);
    
    p_FAR(pn) =  p_FAR(pn) / length(result);
    p_FAR_w(pn) = p_FAR_w(pn) / length(result);
    %% Apply Priors to Theoretical Results
    p_MDR_The(pn) = PC*p_MDR_The(pn)/mm;
    p_MDR_w_The(pn) = PC*p_MDR_w_The(pn)/mm;
    
    p_FAR_The(pn) = PS*p_FAR_The(pn)/mn;
    p_FAR_w_The(pn) = PS*p_FAR_w_The(pn)/mn;
    %% --------------- Normalize SNR/SJNR --------------
    SNR_Theo(pn) = SNR_Theo(pn) / ( mm );
    SNR_w_Theo(pn) = SNR_w_Theo(pn) / ( mm );
    
    % SJNR Calculation (Signal / (Interference + Noise))
    % [PDF: Eq. (14)]
    SNR(pn) =  Sig_power / ( ACA_power + N_power );                
    SNR_w(pn) =  Sig_power / ( N_power ); 

    SNR_lp(pn) = Sig_power / ( AJ_lowpower + N_power );
    SNR_hp(pn) = Sig_power / ( AJ_highpower + N_power );
end

%% ------------------------ Plot Results ------------------------
figure;
%% Plot Probability Curves (Left Axis)
yyaxis left
semilogy(PdBm, p_FAR_w, 'mh', 'LineWidth', 1.8,'MarkerSize', 7);           
hold on 
semilogy(PdBm, p_FAR_w_The, 'm-', 'LineWidth', 1.8,'MarkerSize', 7);
hold on
semilogy(PdBm, p_FAR, 'cd', 'LineWidth', 1.8,'MarkerSize', 7);           
hold on 
semilogy(PdBm, p_FAR_The, 'c-', 'LineWidth', 1.8,'MarkerSize', 7); 
hold on
semilogy(PdBm, p_MDR_w , 'ro', 'LineWidth', 1.8,'MarkerSize', 7);           
hold on 
semilogy(PdBm, p_MDR_w_The, 'r-', 'LineWidth', 1.8,'MarkerSize', 7);  
hold on
semilogy(PdBm, p_MDR, 'bs', 'LineWidth', 1.8,'MarkerSize', 7);           
hold on 
semilogy(PdBm, p_MDR_The, 'b-', 'LineWidth', 1.8,'MarkerSize', 7); 
ylabel('Error Probabilities at Willie');
hold on

%% Plot SNR/Rate Curves (Right Axis)
yyaxis right
plot(PdBm, log2(1+SNR_w),'kp', 'LineWidth', 1.8,'MarkerSize', 7); 
hold on 
plot(PdBm, log2(1+SNR_w_Theo),'k-.','LineWidth', 1.8,'MarkerSize', 7); 
hold on
plot(PdBm, log2(1+SNR),'gv','LineWidth', 1.8,'MarkerSize', 7); 
hold on 
plot(PdBm, log2(1+SNR_Theo),'g-.','LineWidth', 1.8,'MarkerSize', 7); 
hold on 
plot(PdBm, log2(1+SNR_lp),'>','LineWidth', 1.8,'MarkerSize', 7); 
hold on 
plot(PdBm, log2(1+SNR_hp),'+ ','LineWidth', 1.8,'MarkerSize', 7); 
xlabel('Transmit Power at Alice [dBm]');
ylabel('Achievable Rate at Bob [bits/symbol/user]');

grid on;
legend('FAR W/O DRIS','Theoretical FAR','FAR W/ DRIS','Theoretical FAR',...
    'MDR W/O DRIS','Theoretical MDR','MDR W/ DRIS','Theoretical MDR',...
    'SJNR W/O DRIS','Theoretical SJNR','SJNR W/ DRIS','Theoretical SJNR','AJ @ -7dBm','AJ @ 3dBm' );