% =========================================================================   
% (c) 2016 Ronald Nissel, ronald.nissel@gmail.com
% ========================================================================= 
% This script simulates an FBMC and OFDM transmission over a doubly-flat
% channel, including channel estimation. The pilot symbol aided channel 
% estimation in FBMC is based on R. Nissel, M. Rupp, "On Pilot-Symbol Aided
% Channel Estimation in FBMC-OQAM", IEEE ICASSP, 2016.

clear; close all;
addpath('./Theory');

M_SNR_OFDM_dB = [0:5:30];           % Signal-to-Noise Ratio in dB
NrRepetitions = 1000;               % Number of Monte Carlo repetition (different channel realizations)      
QAM_ModulationOrder = 16;           % QAM signal constellation order, 4, 16, 64, 256, 1024,...

%% FBMC Object
FBMC = Modulation.FBMC(...
    12,...                          % Number subcarriers
    30,...                          % Number FBMC symbols
    15e3,...                        % Subcarrier spacing (Hz)
    15e3*14*12,...                  % Sampling rate (Samples/s)
    15e3*20,...                     % Intermediate frequency first subcarrier (Hz)
    false,...                       % Transmit real valued signal 
    'Hermite-OQAM',...              % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
    8, ...                          % Overlapping factor (corresponding to the prototype filter length)
    0, ...                          % Initial phase shift
    true ...                        % Polyphase implementation
    );

%% OFDM Object
OFDM = Modulation.OFDM(...
    12,...                          % Number subcarriers
    15,...                          % Number OFDM Symbols
    15e3,...                        % Subcarrier spacing (Hz)
    15e3*14*12,...                  % Sampling rate (Samples/s)
    15e3*20,...                     % Intermediate frequency first subcarrier (Hz)
    false,...                       % Transmit real valued signal
    0, ...                          % Cyclic prefix length (s), LTE: 1/15e3/14
    (8-1/2)*1/15e3*1/2 ...          % Zero guard length (s)
    );

%% PAM and QAM Object
PAM = Modulation.SignalConstellation(sqrt(QAM_ModulationOrder),'PAM');
QAM = Modulation.SignalConstellation(QAM_ModulationOrder,'QAM');

%% Channel Estimation Objects
ChannelEstimation_OFDM = ChannelEstimation.PilotSymbolAidedChannelEstimation(...
    'Diamond',...                           % Pilot pattern
    [...                                    % Matrix that represents the pilot pattern parameters
    OFDM.Nr.Subcarriers,...                 % Number of subcarriers
    6; ...                                  % Pilot spacing in the frequency domain
    OFDM.Nr.MCSymbols,...                   % Number of FBMC/OFDM Symbols
    4 ...                                   % Pilot spacing in the time domain
    ],...                                   
    'linear'...                             % Interpolation(Extrapolation) method 'linear','spline','FullAverage,'MovingBlockAverage'
    );
ChannelEstimation_FBMC = ChannelEstimation.PilotSymbolAidedChannelEstimation(...
    'Diamond',...                           % Pilot pattern
    [...                                    % Matrix that represents the pilot pattern parameters
    FBMC.Nr.Subcarriers,...                 % Number of subcarriers
    6; ...                                  % Pilot spacing in the frequency domain
    FBMC.Nr.MCSymbols,...                   % Number of FBMC/OFDM Symbols
    8 ...                                   % Pilot spacing in the time domain
    ],...                                   
    'linear'...                             % Interpolation(Extrapolation) method 'linear','spline','FullAverage,'MovingBlockAverage',...
    );

%% Imaginary Interference Cancellation Objects
AuxiliaryMethod = ChannelEstimation.ImaginaryInterferenceCancellationAtPilotPosition(...
    'Auxiliary', ...                                    % Cancellation method
    ChannelEstimation_FBMC.GetAuxiliaryMatrix(1), ...   % PilotMatrix
    FBMC.GetFBMCMatrix, ...                             % Imaginary interference matrix
    16, ...                                             % Cancel 16 closest interferers
    2 ...                                               % Pilot to data power offset
    );
CodingMethod = ChannelEstimation.ImaginaryInterferenceCancellationAtPilotPosition(...
    'Coding', ...                                       % Cancellation method
    ChannelEstimation_FBMC.PilotMatrix, ...             % PilotMatrix
    FBMC.GetFBMCMatrix, ...                             % Imaginary interference matrix
    16, ...                                             % Cancel 16 closest interferers
    2 ...                                               % Pilot to data power offset
    );

BER_FBMC_Aux = nan(length(M_SNR_OFDM_dB),NrRepetitions);
BER_FBMC_Cod = nan(length(M_SNR_OFDM_dB),NrRepetitions);
BER_FBMC_perfect = nan(length(M_SNR_OFDM_dB),NrRepetitions);
BER_OFDM = nan(length(M_SNR_OFDM_dB),NrRepetitions);
BER_OFDM_perfect = nan(length(M_SNR_OFDM_dB),NrRepetitions);
for i_rep = 1:NrRepetitions
    for i_SNR = 1:length(M_SNR_OFDM_dB)
        SNR_OFDM_dB = M_SNR_OFDM_dB(i_SNR);
        Pn_time = OFDM.PHY.SamplingRate/(OFDM.PHY.SubcarrierSpacing*OFDM.Nr.Subcarriers)*10^(-SNR_OFDM_dB/10); 

        %% Generate Random BitStream
        BinaryDataStream_FBMC_Aux = randi([0 1],AuxiliaryMethod.NrDataSymbols*log2(PAM.ModulationOrder),1);
        BinaryDataStream_FBMC_Cod = randi([0 1],CodingMethod.NrDataSymbols*log2(PAM.ModulationOrder),1);
        BinaryDataStream_OFDM     = randi([0 1],(OFDM.Nr.Subcarriers*OFDM.Nr.MCSymbols-ChannelEstimation_OFDM.NrPilotSymbols)*log2(QAM.ModulationOrder),1);

        %% Transmitted Data Symbols
        xD_FBMC_Aux = PAM.Bit2Symbol(BinaryDataStream_FBMC_Aux);
        xD_FBMC_Cod = PAM.Bit2Symbol(BinaryDataStream_FBMC_Cod);
        xD_OFDM     = QAM.Bit2Symbol(BinaryDataStream_OFDM);

        %% Transmitted Pilot Symbols
        xP_FBMC = PAM.SymbolMapping(randi(PAM.ModulationOrder,[ChannelEstimation_FBMC.NrPilotSymbols 1]));
        xP_FBMC = xP_FBMC./abs(xP_FBMC);
        xP_OFDM = QAM.SymbolMapping(randi(QAM.ModulationOrder,[ChannelEstimation_OFDM.NrPilotSymbols 1]));
        xP_OFDM = xP_OFDM./abs(xP_OFDM); 

        %% Transmitted Symbols
        x_FBMC_Aux = reshape(AuxiliaryMethod.PrecodingMatrix*[xP_FBMC;xD_FBMC_Aux],[FBMC.Nr.Subcarriers FBMC.Nr.MCSymbols]);
        x_FBMC_Cod = reshape(CodingMethod.PrecodingMatrix*[xP_FBMC;xD_FBMC_Cod],[FBMC.Nr.Subcarriers FBMC.Nr.MCSymbols]);
        x_OFDM = nan(OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols);
        x_OFDM(ChannelEstimation_OFDM.PilotMatrix==1) = xP_OFDM;
        x_OFDM(ChannelEstimation_OFDM.PilotMatrix==0) = xD_OFDM;

        %% Transmitted FBMC Signal (time domain)
        s_FBMC_Aux = FBMC.Modulation(x_FBMC_Aux); 
        s_FBMC_Cod = FBMC.Modulation(x_FBMC_Cod);     
        s_OFDM     = OFDM.Modulation(x_OFDM);

        %% Channel (doubly flat fading and AWGN)     
        h = sqrt(1/2)*(randn+1j*randn);
    %     h = 1; % Pure AWGN
        n_FBMC = sqrt(Pn_time/2)*(randn(size(s_FBMC_Cod))+1j*randn(size(s_FBMC_Cod)));
        n_OFDM = sqrt(Pn_time/2)*(randn(size(s_OFDM))+1j*randn(size(s_OFDM)));

        r_FBMC_Aux = h*s_FBMC_Aux + n_FBMC;
        r_FBMC_Cod = h*s_FBMC_Cod + n_FBMC;   
        r_OFDM     = h*s_OFDM + n_OFDM;

        %% Demodulate OFDM and FBMC signal
        y_FBMC_Aux = FBMC.Demodulation(r_FBMC_Aux);
        y_FBMC_Cod = FBMC.Demodulation(r_FBMC_Cod);
        y_OFDM     = OFDM.Demodulation(r_OFDM);

        %% LS channel estimates at pilot positions
        hP_LS_FBMC_Aux = y_FBMC_Aux(ChannelEstimation_FBMC.PilotMatrix==1)./xP_FBMC/sqrt(AuxiliaryMethod.PilotToDataPowerOffset*AuxiliaryMethod.DataPowerReduction);
        hP_LS_FBMC_Cod = y_FBMC_Cod(ChannelEstimation_FBMC.PilotMatrix==1)./xP_FBMC/sqrt(CodingMethod.PilotToDataPowerOffset);
        hP_LS_OFDM     = y_OFDM(ChannelEstimation_OFDM.PilotMatrix==1)./xP_OFDM;

        %% Channel Estimation using Interpolation
        h_FBMC_Aux = ChannelEstimation_FBMC.ChannelInterpolation(hP_LS_FBMC_Aux);
        h_FBMC_Cod = ChannelEstimation_FBMC.ChannelInterpolation(hP_LS_FBMC_Cod);
        h_OFDM     = ChannelEstimation_OFDM.ChannelInterpolation(hP_LS_OFDM);

        %% Equalized received symbols at data position
        y_EQ_FBMC_Aux = real(y_FBMC_Aux(AuxiliaryMethod.PilotMatrix==0)./h_FBMC_Aux(AuxiliaryMethod.PilotMatrix==0)/sqrt(AuxiliaryMethod.DataPowerReduction));
        y_EQ_FBMC_Cod = real(CodingMethod.PrecodingMatrix(:,CodingMethod.NrPilotSymbols+1:end)'*(y_FBMC_Cod(:)./h_FBMC_Cod(:)));
        y_EQ_FBMC_perfect = real(CodingMethod.PrecodingMatrix(:,CodingMethod.NrPilotSymbols+1:end)'*(y_FBMC_Cod(:)./h));

        y_EQ_OFDM = y_OFDM(ChannelEstimation_OFDM.PilotMatrix==0)./h_OFDM(ChannelEstimation_OFDM.PilotMatrix==0);
        y_EQ_OFDM_perfect = y_OFDM(ChannelEstimation_OFDM.PilotMatrix==0)./h;

        %% Detect BitStream
        DetectedBitStream_FBMC_Aux = PAM.Symbol2Bit(real(y_EQ_FBMC_Aux(:)));
        DetectedBitStream_FBMC_Cod = PAM.Symbol2Bit(real(y_EQ_FBMC_Cod(:)));
        DetectedBitStream_FBMC_perfect = PAM.Symbol2Bit(real(y_EQ_FBMC_perfect(:)));

        DetectedBitStream_OFDM = QAM.Symbol2Bit(y_EQ_OFDM(:));
        DetectedBitStream_OFDM_perfect = QAM.Symbol2Bit(y_EQ_OFDM_perfect(:));

        %% Calculate BER
        BER_FBMC_Aux(i_SNR,i_rep) = mean(BinaryDataStream_FBMC_Aux~=DetectedBitStream_FBMC_Aux);
        BER_FBMC_Cod(i_SNR,i_rep) = mean(BinaryDataStream_FBMC_Cod~=DetectedBitStream_FBMC_Cod);
        BER_FBMC_perfect(i_SNR,i_rep) = mean(BinaryDataStream_FBMC_Cod~=DetectedBitStream_FBMC_perfect);

        BER_OFDM(i_SNR,i_rep) = mean(BinaryDataStream_OFDM~=DetectedBitStream_OFDM);   
        BER_OFDM_perfect(i_SNR,i_rep) = mean(BinaryDataStream_OFDM~=DetectedBitStream_OFDM_perfect);   

    end
    
    if mod(i_rep,100)==0
       disp([int2str(i_rep/NrRepetitions*100) '%']);
    end
end

%% Theoretical BEP for perfect channel knowledge 
% BEP_4QAM = 1/2-1./(2*sqrt(2*(1+10.^(-M_SNR_OFDM_dB/10))-1));
M_SNR_OFDM_dB_morePoints = min(M_SNR_OFDM_dB):0.5:max(M_SNR_OFDM_dB);
BEP_perfect = BitErrorProbabilityDoublyFlatRayleigh(M_SNR_OFDM_dB_morePoints,QAM.SymbolMapping,QAM.BitMapping);


%% Plot BER and BEP
figure();
semilogy(M_SNR_OFDM_dB,mean(BER_FBMC_Aux,2),'red -o');
hold on;
semilogy(M_SNR_OFDM_dB,mean(BER_FBMC_Cod,2),'blue -o');
semilogy(M_SNR_OFDM_dB,mean(BER_OFDM,2),'black -o'); 
semilogy(M_SNR_OFDM_dB,mean(BER_FBMC_perfect,2),'blue -x');
semilogy(M_SNR_OFDM_dB,mean(BER_OFDM_perfect,2),'black -x');
semilogy(M_SNR_OFDM_dB_morePoints,BEP_perfect','black');
xlabel('SNR for OFDM (dB)'); 
ylabel('BER, BEP');
legend('Simulation: FBMC Auxiliary','Simulation: FBMC Coding','Simulation: OFDM', 'Simulation FBMC perfect CSI', 'Simulation OFDM perfect CSI','Theory perfect CSI','Location','SouthWest');

%% Plot Pilot Pattern
figure();
ChannelEstimation_OFDM.PlotPilotPattern;
title('OFDM');
figure();
ChannelEstimation_FBMC.PlotPilotPattern(AuxiliaryMethod.PilotMatrix)
title('FBMC Auxiliary');
figure();
ChannelEstimation_FBMC.PlotPilotPattern(-(CodingMethod.ConsideredInterferenceMatrix<0)+(CodingMethod.ConsideredInterferenceMatrix>0))
title('FBMC Coding');

%% Calculate and Plot Expected Transmit Power Over Time
[Power_FBMC_Aux,t_FBMC] = FBMC.PlotTransmitPower(AuxiliaryMethod.PrecodingMatrix*AuxiliaryMethod.PrecodingMatrix');
[Power_FBMC_Cod,~] = FBMC.PlotTransmitPower(CodingMethod.PrecodingMatrix*CodingMethod.PrecodingMatrix');
[Power_OFDM,t_OFDM] = OFDM.PlotTransmitPower;
figure();
plot(t_FBMC,Power_FBMC_Aux,'red');
hold on;
plot(t_FBMC,Power_FBMC_Cod,'blue');
plot(t_OFDM,Power_OFDM,'black ');
legend({'FBMC Auxiliary','FBMC Coding','OFDM'});
ylabel('Transmit Power');
xlabel('Time(s)');

%% Calculate Power Spectral Density
[PSD_FBMC_Aux,t_FBMC] = FBMC.PlotPowerSpectralDensity(AuxiliaryMethod.PrecodingMatrix*AuxiliaryMethod.PrecodingMatrix');
[PSD_FBMC_Cod,~] = FBMC.PlotPowerSpectralDensity(CodingMethod.PrecodingMatrix*CodingMethod.PrecodingMatrix');
[PSD_OFDM,t_OFDM] = OFDM.PlotPowerSpectralDensity;
figure();
plot(t_FBMC,10*log10(PSD_FBMC_Aux),'red');
hold on;
plot(t_FBMC,10*log10(PSD_FBMC_Cod),'blue');
plot(t_OFDM,10*log10(PSD_OFDM),'black ');
legend({'FBMC Auxiliary','FBMC Coding','OFDM'});
ylabel('Power Spectral Density (dB)');
xlabel('Frequency (Hz)');

