% =========================================================================   
% (c) 2018 Ronald Nissel, ronald.nissel@gmail.com
% =========================================================================  
% Allows to reproduce Figure 2, 3, 4 and 5 of "Doubly-Selective Channel 
% Estimation in FBMC-OQAM and OFDM Systems", R. Nissel, et al, IEEE VTC
% Fall, 2018. In particular, this script simulates an FBMC and OFDM
% transmission over a doubly-selective channel, including doubly-selective
% MMSE channel estimation and interference cancellation. Note that, to 
% truly reproduce the figures, the lines 41-46 must be uncommented.


clear; close all;
addpath('./Theory');


%% Parameters
% Simulation
M_SNR_dB                  = [10:5:40];              % Signal-to-Noise Ratio in dB
NrRepetitions             = 25;                     % Number of Monte Carlo repetition (different channel realizations)             
ZeroThresholdSparse       = 8;                      % Set some matrix values, smaller than "10^(-ZeroThresholdSparse)", to zero.
PlotIterationStepsSNRdB   = 35;                     % Plot the BER over iteration step for an SNR of 35dB. 

% FBMC and OFDM parameters
L                         = 12*2;                   % Number of Subcarriers, one resource block consists of 12 subcarriers (and 0.5ms in time)
F                         = 15e3;                   % Subcarrier Spacing in Hz
SamplingRate              = F*12*2;                 % Sampling rate in (Samples/s)
NrSubframes               = 1;                      % Number of subframes. One subframe requires 1ms for F=15kHz.                             
QAM_ModulationOrder       = 256;                    % QAM signal constellation order, 4, 16, 64, 256, 1024,...

% Channel estimation parameters
PilotToDataPowerOffset    = 2;                      % Pilot to data power offset for OFDM. In FBMC data spreading, the power offset is twice this number. 
PilotToDataPowerOffsetAux = 4.685;                  % Pilot to data power offset for FBMC, auxiliary method.
NrIterations              = 4;                      % Number of iterations for interference cancellation scheme.

% Channel
Velocity_kmh              = 500;                    % Velocity in km/h. Note that [mph]*1.6=[kmh] and [m/s]*3.6=[kmh]        
PowerDelayProfile         = 'VehicularA';           % Channel model, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate) 


% ###########################################################################
% % In the paper:
% M_SNR_dB                = [10:2:40];
% PlotIterationStepsSNRdB = 32;
% NrRepetitions           = 1000;
% SamplingRate            = F*14*14;
% NrSubframes             = 2;
% ###########################################################################


%% FBMC Object
FBMC = Modulation.FBMC(...
    L,...                               % Number subcarriers
    30*NrSubframes,...                  % Number FBMC symbols
    F,...                               % Subcarrier spacing (Hz)
    SamplingRate,...                    % Sampling rate (Samples/s)
    0,...                               % Intermediate frequency first subcarrier (Hz)
    false,...                           % Transmit real valued signal
    'Hermite-OQAM',...                  % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
    8, ...                              % Overlapping factor (also determines oversampling in the frequency domain)
    0, ...                              % Initial phase shift
    true ...                            % Polyphase implementation
    );


%% OFDM Object (Add zeroes to the OFDM signal so that it fits the FBMC signal)
ZeroGuardTimeLength = ((FBMC.Nr.SamplesTotal-(round((1/15e3/14)*SamplingRate)+round(SamplingRate/15e3))*14*NrSubframes)/2)/SamplingRate;
OFDM = Modulation.OFDM(...
    L,...                           % Number Subcarriers
    14*NrSubframes,...              % Number OFDM Symbols
    F,...                           % Subcarrier spacing (Hz)
    SamplingRate,...                % Sampling rate (Samples/s)
    0,...                           % Intermediate frequency first subcarrier (Hz)
    false,...                       % Transmitreal valued signal
    1/15e3/14, ...                  % Cyclic prefix length (s)
    ZeroGuardTimeLength ...         % Zero guard length (s)
    );

%% Check Number of Samples
if  OFDM.Nr.SamplesTotal~=FBMC.Nr.SamplesTotal
   error('Total number of samples must be the same for OFDM and FBMC.');
end
N = OFDM.Nr.SamplesTotal;


%% PAM and QAM Object
PAM = Modulation.SignalConstellation(sqrt(QAM_ModulationOrder),'PAM');
QAM = Modulation.SignalConstellation(QAM_ModulationOrder,'QAM');


%% Pilot Matrices, 0="Data", 1=Pilot
PilotMatrix_OFDM = zeros(OFDM.Nr.Subcarriers,14);
PilotMatrix_OFDM(2:2*6:end,2:2*3.5:end)=1;
PilotMatrix_OFDM(5:2*6:end,6:2*3.5:end)=1;
PilotMatrix_OFDM(8:2*6:end,2:2*3.5:end)=1;
PilotMatrix_OFDM(11:2*6:end,6:2*3.5:end)=1;
PilotMatrix_OFDM = repmat(PilotMatrix_OFDM,[1 NrSubframes]);

PilotMatrix_FBMC =  zeros(FBMC.Nr.Subcarriers,30);
PilotMatrix_FBMC(2:12:end,3:16:end) = 1;
PilotMatrix_FBMC(5:12:end,11:16:end) = 1;
PilotMatrix_FBMC(8:12:end,3+1:16:end) = 1;
PilotMatrix_FBMC(11:12:end,11+1:16:end) = 1;
PilotMatrix_FBMC = repmat(PilotMatrix_FBMC,[1 NrSubframes]);


AuxilaryPilotMatrix_FBMC = PilotMatrix_FBMC;
[a,b] = find(PilotMatrix_FBMC);
for i_pilot = 1:length(a)
    AuxilaryPilotMatrix_FBMC(a(i_pilot)+1,b(i_pilot))=-1;
    AuxilaryPilotMatrix_FBMC(a(i_pilot)-1,b(i_pilot))=-1;
    AuxilaryPilotMatrix_FBMC(a(i_pilot),b(i_pilot)+1)=-1;
    AuxilaryPilotMatrix_FBMC(a(i_pilot),b(i_pilot)-1)=-1;
end

%% Cancel Imaginary Interference At Pilot Position Object (Precoding Matrix)
AuxiliaryMethod = ChannelEstimation.ImaginaryInterferenceCancellationAtPilotPosition(...
    'Auxiliary', ...                                    % Cancellation method
    AuxilaryPilotMatrix_FBMC, ...                       % PilotAndAuxiliaryMatrix
    FBMC.GetFBMCMatrix, ...                             % Imaginary interference matrix
    28, ...                                             % Cancel 28 closest interferers
    PilotToDataPowerOffsetAux ...                       % Pilot to data power offset
    );
CodingMethod = ChannelEstimation.ImaginaryInterferenceCancellationAtPilotPosition(...
    'Coding', ...                                       % Cancellation method
    PilotMatrix_FBMC, ...                               % PilotMatrix
    FBMC.GetFBMCMatrix, ...                             % Imaginary interference matrix
    20, ...                                             % Cancel 28 closest interferers
    2*PilotToDataPowerOffset ...                        % Pilot to data power offset
    );

NrPilotSymbols = sum(PilotMatrix_OFDM(:)==1);
NrDataSymbols_OFDM = sum(PilotMatrix_OFDM(:)==0);

PilotMapping_OFDM = zeros(numel(PilotMatrix_OFDM));
PilotMapping_OFDM(PilotMatrix_OFDM(:)==1,1:NrPilotSymbols) = sqrt(PilotToDataPowerOffset)*eye(NrPilotSymbols);
PilotMapping_OFDM(PilotMatrix_OFDM(:)==0,NrPilotSymbols+1:end) = eye(NrDataSymbols_OFDM);
PilotMapping_OFDM = PilotMapping_OFDM/sqrt(mean(diag(PilotMapping_OFDM*PilotMapping_OFDM')));
DataPowerReduction_OFDM = numel(PilotMatrix_OFDM)/(NrPilotSymbols*PilotToDataPowerOffset+NrDataSymbols_OFDM);

Kappa_Aux  = AuxiliaryMethod.PilotToDataPowerOffset*AuxiliaryMethod.DataPowerReduction;
Kappa_Cod  = CodingMethod.PilotToDataPowerOffset*CodingMethod.DataPowerReduction;
Kappa_OFDM = PilotToDataPowerOffset*DataPowerReduction_OFDM;

%% Evaluate only the BER near the center. i.e., ignore edges 
ConsideredTimeFrequencyPositions_FBMC = zeros(size(PilotMatrix_FBMC));
ConsideredTimeFrequencyPositions_FBMC(5:end-4,11:end-10)= 1;

ConsideredTimeFrequencyPositions_OFDM = zeros(size(PilotMatrix_OFDM));
ConsideredTimeFrequencyPositions_OFDM(5:end-4,6:end-5)= 1;

for i_lk = 1: AuxiliaryMethod.NrDataSymbols
    xD_Temp = zeros(AuxiliaryMethod.NrDataSymbols,1);
    xD_Temp(i_lk) = 1;
    x_FBMC_Aux = AuxiliaryMethod.PrecodingMatrix*[zeros(NrPilotSymbols,1);xD_Temp];    
    ConsideredDataPositions_FBMC_Aux(i_lk) = sum(abs(x_FBMC_Aux(1==ConsideredTimeFrequencyPositions_FBMC.*not(AuxilaryPilotMatrix_FBMC))))>AuxiliaryMethod.DataPowerReduction*0.9;
end
for i_lk = 1: CodingMethod.NrDataSymbols
    xD_Temp = zeros(CodingMethod.NrDataSymbols,1);
    xD_Temp(i_lk) = 1;
    x_FBMC_Cod = CodingMethod.PrecodingMatrix*[zeros(NrPilotSymbols,1);xD_Temp];      
    ConsideredDataPositions_FBMC_Cod(i_lk) =not(any(x_FBMC_Cod(not(ConsideredTimeFrequencyPositions_FBMC))));
end
for i_lk = 1: NrDataSymbols_OFDM
    xD_Temp = zeros(NrDataSymbols_OFDM,1);
    xD_Temp(i_lk) = 1;
    x_OFDM = PilotMapping_OFDM*[zeros(NrPilotSymbols,1);xD_Temp];
    ConsideredDataPositions_OFDM(i_lk) = sum(abs(x_OFDM(1==ConsideredTimeFrequencyPositions_OFDM.*not(PilotMatrix_OFDM))))>DataPowerReduction_OFDM*0.9;
end

ConsideredBits_FBMC_Aux = reshape(repmat(ConsideredDataPositions_FBMC_Aux,log2(QAM_ModulationOrder)/2,1),[],1);
ConsideredBits_FBMC_Cod = reshape(repmat(ConsideredDataPositions_FBMC_Cod,log2(QAM_ModulationOrder)/2,1),[],1);
ConsideredBits_OFDM     = reshape(repmat(ConsideredDataPositions_OFDM,log2(QAM_ModulationOrder),1),[],1);


%% Channel Model Object
ChannelModel = Channel.FastFading(...
    SamplingRate,...                                   % Sampling rate (Samples/s)
    PowerDelayProfile,...                              % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate)
    N,...                                              % Number of total samples
    Velocity_kmh/3.6*2.5e9/2.998e8,...                 % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8
    'Jakes',...                                        % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large
    200,...                                            % Number of paths for the WSSUS process. Only relevant for a 'Jakes' and 'Uniform' Doppler spectrum
    1,...                                              % Number of transmit antennas
    1,...                                              % Number of receive antennas
    1 ...                                              % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
    );
R_vecH = ChannelModel.GetCorrelationMatrix;


%% Precalculate Transmit and Receive Matrices
G_FBMC = FBMC.GetTXMatrix;
Q_FBMC = (FBMC.GetRXMatrix)';

G_OFDM = OFDM.GetTXMatrix;
Q_OFDM = (OFDM.GetRXMatrix)';

GP_FBMC = G_FBMC(:,PilotMatrix_FBMC(:)==1);
GP_OFDM = G_OFDM(:,PilotMatrix_OFDM(:)==1);

QP_FBMC = Q_FBMC(:,PilotMatrix_FBMC(:)==1);
QP_OFDM = Q_OFDM(:,PilotMatrix_OFDM(:)==1);

G_Aux = G_FBMC*AuxiliaryMethod.PrecodingMatrix;
G_Cod = G_FBMC*CodingMethod.PrecodingMatrix;
G_OFDM_PilotMapping = G_OFDM*PilotMapping_OFDM;


%% Calculate Correlation Matrices
disp('Calculate Correlation Matrix of Pilot Estimates (no Noise, no Interference)...');
R_hP_FBMC  = nan(NrPilotSymbols,NrPilotSymbols);
R_hP_OFDM = nan(NrPilotSymbols,NrPilotSymbols);
for j_pilot = 1:NrPilotSymbols
    R_hP_FBMC(:,j_pilot) = sum((QP_FBMC'*reshape(R_vecH*kron(GP_FBMC(:,j_pilot).',QP_FBMC(:,j_pilot)')',N,N)).*(GP_FBMC.'),2);    
    R_hP_OFDM(:,j_pilot) = sum((QP_OFDM'*reshape(R_vecH*kron(GP_OFDM(:,j_pilot).',QP_OFDM(:,j_pilot)')',N,N)).*(GP_OFDM.'),2);
end


disp('Calculate Correlation Matrix of Pilot Estimates (no Interference)...');
R_hP_est_noNoise_FBMC_Aux = R_hP_FBMC;
R_hP_est_noNoise_FBMC_Cod = R_hP_FBMC;
R_hP_est_noNoise_OFDM = R_hP_OFDM;
for i_pilots = 1: NrPilotSymbols
    % FBMC Auxiliary Method, Similar as Equation (13) but computationally more efficient
    Temp = kron(sparse(eye(N)),QP_FBMC(:,i_pilots)')/sqrt(Kappa_Aux);
    R_hP_est_noNoise_FBMC_Aux(i_pilots,i_pilots)=abs(sum(sum((G_Aux.'*(Temp*R_vecH*Temp')).*G_Aux',2)));     

    % FBMC Coding
    Temp = kron(sparse(eye(N)),QP_FBMC(:,i_pilots)')/sqrt(Kappa_Cod);
    R_hP_est_noNoise_FBMC_Cod(i_pilots,i_pilots)=abs(sum(sum((G_Cod.'*(Temp*R_vecH*Temp')).*G_Cod',2)));     

    % OFDM
    Temp = kron(sparse(eye(N)),QP_OFDM(:,i_pilots)')/sqrt(Kappa_OFDM);
    R_hP_est_noNoise_OFDM(i_pilots,i_pilots)=abs(sum(sum((G_OFDM_PilotMapping.'*(Temp*R_vecH*Temp')).*G_OFDM_PilotMapping',2)));           
end 


disp('Calculate Correlation Matrix of Pilot Estimates...');
R_hP_est_FBMC_Aux = repmat(R_hP_est_noNoise_FBMC_Aux,[1 1 length(M_SNR_dB)]);
R_hP_est_FBMC_Cod = repmat(R_hP_est_noNoise_FBMC_Cod,[1 1 length(M_SNR_dB)]);
R_hP_est_OFDM = repmat(R_hP_est_noNoise_OFDM,[1 1 length(M_SNR_dB)]);
for i_SNR = 1:length(M_SNR_dB)
    SNR_dB = M_SNR_dB(i_SNR);
    Pn_time = SamplingRate/(F*L)*10^(-SNR_dB/10);   
    for i_pilots = 1: NrPilotSymbols
        R_hP_est_FBMC_Aux(i_pilots,i_pilots,i_SNR)=R_hP_est_noNoise_FBMC_Aux(i_pilots,i_pilots)+Pn_time*QP_FBMC(:,i_pilots)'*QP_FBMC(:,i_pilots)/(Kappa_Aux);     
        R_hP_est_FBMC_Cod(i_pilots,i_pilots,i_SNR)=R_hP_est_noNoise_FBMC_Cod(i_pilots,i_pilots)+Pn_time*QP_FBMC(:,i_pilots)'*QP_FBMC(:,i_pilots)/(Kappa_Cod);     
        R_hP_est_OFDM(i_pilots,i_pilots,i_SNR)=R_hP_est_noNoise_OFDM(i_pilots,i_pilots)+Pn_time*QP_OFDM(:,i_pilots)'*QP_OFDM(:,i_pilots)/Kappa_OFDM;           
    end    
end

R_hP_est_noInterference_FBMC_Aux = R_hP_est_FBMC_Aux-repmat((R_hP_est_noNoise_FBMC_Aux-R_hP_FBMC),[1 1 length(M_SNR_dB)]);
R_hP_est_noInterference_FBMC_Cod = R_hP_est_FBMC_Cod-repmat((R_hP_est_noNoise_FBMC_Cod-R_hP_FBMC),[1 1 length(M_SNR_dB)]);
R_hP_est_noInterference_OFDM = R_hP_est_OFDM-repmat((R_hP_est_noNoise_OFDM-R_hP_OFDM),[1 1 length(M_SNR_dB)]);


disp('Calculate Correlation Matrix between Transmission-Matrix D and Pilot Estimates...');
R_Dij_hP_FBMC = sparse(size(G_FBMC,2)^2,NrPilotSymbols);
R_Dij_hP_OFDM = sparse(size(G_OFDM,2)^2,NrPilotSymbols);
for i_pilot = 1: NrPilotSymbols
    R_Dij_hP_FBMC_Temp = reshape(Q_FBMC'*reshape(R_vecH*kron(GP_FBMC(:,i_pilot).',QP_FBMC(:,i_pilot)')',N,N)*G_FBMC,[],1);
    R_Dij_hP_OFDM_Temp = reshape(Q_OFDM'*reshape(R_vecH*kron(GP_OFDM(:,i_pilot).',QP_OFDM(:,i_pilot)')',N,N)*G_OFDM,[],1);
    
    R_Dij_hP_FBMC_Temp(abs(R_Dij_hP_FBMC_Temp)<10^(-ZeroThresholdSparse))=0;
    R_Dij_hP_OFDM_Temp(abs(R_Dij_hP_OFDM_Temp)<10^(-ZeroThresholdSparse))=0;
     
    R_Dij_hP_FBMC(:,i_pilot) = R_Dij_hP_FBMC_Temp;
    R_Dij_hP_OFDM(:,i_pilot) = R_Dij_hP_OFDM_Temp;
end
clear R_vecH;

%% Calculate SIR at Pilot Positions
SIR_P_FBMC_Aux_dB = 10*log10(trace(abs(R_hP_FBMC))./trace(abs(R_hP_est_noNoise_FBMC_Aux-R_hP_FBMC)));
SIR_P_FBMC_Cod_dB = 10*log10(trace(abs(R_hP_FBMC))./trace(abs(R_hP_est_noNoise_FBMC_Cod-R_hP_FBMC)));
SIR_P_OFDM_dB     = 10*log10(trace(abs(R_hP_OFDM))./trace(abs(R_hP_est_noNoise_OFDM-R_hP_OFDM)));


%% Calculate MMSE Estimation Matrix Using Correlation
disp('Calculate MMSE Solution ...');
W_MMSE_FBMC_Aux = sparse(size(G_FBMC,2)*size(G_FBMC,2)*NrPilotSymbols,length(M_SNR_dB));
W_MMSE_FBMC_Cod = sparse(size(G_FBMC,2)*size(G_FBMC,2)*NrPilotSymbols,length(M_SNR_dB));
W_MMSE_OFDM     = sparse(size(G_OFDM,2)*size(G_OFDM,2)*NrPilotSymbols,length(M_SNR_dB));
for i_SNR = 1:length(M_SNR_dB)
    W_MMSE_FBMC_Aux_Temp = reshape(R_Dij_hP_FBMC*pinv(R_hP_est_FBMC_Aux(:,:,i_SNR)),size(G_FBMC,2)*size(G_FBMC,2)*NrPilotSymbols,1);
    W_MMSE_FBMC_Cod_Temp = reshape(R_Dij_hP_FBMC*pinv(R_hP_est_FBMC_Cod(:,:,i_SNR)),size(G_FBMC,2)*size(G_FBMC,2)*NrPilotSymbols,1);
    W_MMSE_OFDM_Temp     = reshape(R_Dij_hP_OFDM*pinv(R_hP_est_OFDM(:,:,i_SNR)),size(G_OFDM,2)*size(G_OFDM,2)*NrPilotSymbols,1);
 
    W_MMSE_FBMC_Aux_Temp(abs(W_MMSE_FBMC_Aux_Temp)<10^(-ZeroThresholdSparse))=0;
    W_MMSE_FBMC_Cod_Temp(abs(W_MMSE_FBMC_Cod_Temp)<10^(-ZeroThresholdSparse))=0;
    W_MMSE_OFDM_Temp(abs(W_MMSE_OFDM_Temp)<10^(-ZeroThresholdSparse))=0;
    
    W_MMSE_FBMC_Aux(:,i_SNR) = W_MMSE_FBMC_Aux_Temp;
    W_MMSE_FBMC_Cod(:,i_SNR) = W_MMSE_FBMC_Cod_Temp;
    W_MMSE_OFDM(:,i_SNR)     = W_MMSE_OFDM_Temp;
end

%% Calculate MMSE Estimation if no Interference is present at the pilot positions
disp('Calculate MMSE Solution (no Interference)...');
W_MMSE_noInterference_FBMC_Aux = sparse(size(G_FBMC,2)*size(G_FBMC,2)*NrPilotSymbols,length(M_SNR_dB));
W_MMSE_noInterference_FBMC_Cod = sparse(size(G_FBMC,2)*size(G_FBMC,2)*NrPilotSymbols,length(M_SNR_dB));
W_MMSE_noInterference_OFDM     = sparse(size(G_OFDM,2)*size(G_OFDM,2)*NrPilotSymbols,length(M_SNR_dB));
for i_SNR = 1:length(M_SNR_dB)
    W_MMSE_noInterference_FBMC_Aux_Temp = reshape(R_Dij_hP_FBMC*pinv(R_hP_est_noInterference_FBMC_Aux(:,:,i_SNR)),size(G_FBMC,2)*size(G_FBMC,2)*NrPilotSymbols,1);
    W_MMSE_noInterference_FBMC_Cod_Temp = reshape(R_Dij_hP_FBMC*pinv(R_hP_est_noInterference_FBMC_Cod(:,:,i_SNR)),size(G_FBMC,2)*size(G_FBMC,2)*NrPilotSymbols,1);
    W_MMSE_noInterference_OFDM_Temp = reshape(R_Dij_hP_OFDM*pinv(R_hP_est_noInterference_OFDM(:,:,i_SNR)),size(G_OFDM,2)*size(G_OFDM,2)*NrPilotSymbols,1);

    W_MMSE_noInterference_FBMC_Aux_Temp(abs(W_MMSE_noInterference_FBMC_Aux_Temp)<10^(-ZeroThresholdSparse))=0;
    W_MMSE_noInterference_FBMC_Cod_Temp(abs(W_MMSE_noInterference_FBMC_Cod_Temp)<10^(-ZeroThresholdSparse))=0;
    W_MMSE_noInterference_OFDM_Temp(abs(W_MMSE_noInterference_OFDM_Temp)<10^(-ZeroThresholdSparse))=0;
    
    W_MMSE_noInterference_FBMC_Aux(:,i_SNR) = W_MMSE_noInterference_FBMC_Aux_Temp;
    W_MMSE_noInterference_FBMC_Cod(:,i_SNR) = W_MMSE_noInterference_FBMC_Cod_Temp;
    W_MMSE_noInterference_OFDM(:,i_SNR)     = W_MMSE_noInterference_OFDM_Temp;    
end


%% Theoretical BEP for a Doubly Flat Rayleigh Channel
M_SNR_dB_morePoints = min(M_SNR_dB):1:max(M_SNR_dB);
BitErrorProbability = BitErrorProbabilityDoublyFlatRayleigh(M_SNR_dB_morePoints,QAM.SymbolMapping,QAM.BitMapping);


%% Preallocate for Parfor
BER_FBMC_Aux_PerfectCSI_InterferenceCancellation        = nan(length(M_SNR_dB),NrRepetitions,NrIterations);
BER_FBMC_Aux_PerfectCSI_InterferenceCancellation_NoEdge = nan(length(M_SNR_dB),NrRepetitions,NrIterations);
BER_FBMC_Cod_PerfectCSI_InterferenceCancellation        = nan(length(M_SNR_dB),NrRepetitions,NrIterations);
BER_FBMC_Cod_PerfectCSI_InterferenceCancellation_NoEdge = nan(length(M_SNR_dB),NrRepetitions,NrIterations);
BER_OFDM_PerfectCSI_InterferenceCancellation            = nan(length(M_SNR_dB),NrRepetitions,NrIterations);
BER_OFDM_PerfectCSI_InterferenceCancellation_NoEdge     = nan(length(M_SNR_dB),NrRepetitions,NrIterations);
BER_FBMC_Aux_OneTapEqualizer_PerfectCSI                 = nan(length(M_SNR_dB),NrRepetitions);
BER_FBMC_Aux_OneTapEqualizer_PerfectCSI_NoEdge          = nan(length(M_SNR_dB),NrRepetitions);
BER_FBMC_Cod_OneTapEqualizer_PerfectCSI                 = nan(length(M_SNR_dB),NrRepetitions);
BER_FBMC_Cod_OneTapEqualizer_PerfectCSI_NoEdge          = nan(length(M_SNR_dB),NrRepetitions);
BER_OFDM_OneTapEqualizer_PerfectCSI                     = nan(length(M_SNR_dB),NrRepetitions);
BER_OFDM_OneTapEqualizer_PerfectCSI_NoEdge              = nan(length(M_SNR_dB),NrRepetitions);
BER_FBMC_Aux_InterferenceCancellation                   = nan(length(M_SNR_dB),NrRepetitions,NrIterations);
BER_FBMC_Cod_InterferenceCancellation                   = nan(length(M_SNR_dB),NrRepetitions,NrIterations);
BER_OFDM_InterferenceCancellation                       = nan(length(M_SNR_dB),NrRepetitions,NrIterations);
BER_FBMC_Aux_InterferenceCancellation_NoEdge            = nan(length(M_SNR_dB),NrRepetitions,NrIterations);
BER_FBMC_Cod_InterferenceCancellation_NoEdge            = nan(length(M_SNR_dB),NrRepetitions,NrIterations);
BER_OFDM_InterferenceCancellation_NoEdge                = nan(length(M_SNR_dB),NrRepetitions,NrIterations);
BER_FBMC_Aux_OneTapEqualizer                            = nan(length(M_SNR_dB),NrRepetitions);
BER_FBMC_Aux_OneTapEqualizer_NoEdge                     = nan(length(M_SNR_dB),NrRepetitions);
BER_FBMC_Cod_OneTapEqualizer                            = nan(length(M_SNR_dB),NrRepetitions);
BER_FBMC_Cod_OneTapEqualizer_NoEdge                     = nan(length(M_SNR_dB),NrRepetitions);
BER_OFDM_OneTapEqualizer                                = nan(length(M_SNR_dB),NrRepetitions);
BER_OFDM_OneTapEqualizer_NoEdge                         = nan(length(M_SNR_dB),NrRepetitions);

%% Start Simulation
tic
disp('Monte Carlo Simulation ...');
for i_rep = 1:NrRepetitions
    %% Update Channel
    ChannelModel.NewRealization;
    
    %% Binary Data
    BinaryDataStream_FBMC_Aux = randi([0 1],AuxiliaryMethod.NrDataSymbols*log2(PAM.ModulationOrder),1);
    BinaryDataStream_FBMC_Cod = randi([0 1],CodingMethod.NrDataSymbols*log2(PAM.ModulationOrder),1);    
    BinaryDataStream_OFDM     = randi([0 1],NrDataSymbols_OFDM*log2(QAM.ModulationOrder),1);
       
    %% Data Symbols
    xD_FBMC_Aux = PAM.Bit2Symbol(BinaryDataStream_FBMC_Aux);
    xD_FBMC_Cod = PAM.Bit2Symbol(BinaryDataStream_FBMC_Cod);
    xD_OFDM     = QAM.Bit2Symbol(BinaryDataStream_OFDM);
       
    %% Pilot Symbols
    xP_FBMC = PAM.SymbolMapping(randi(PAM.ModulationOrder,AuxiliaryMethod.NrPilotSymbols,1));
    xP_FBMC = xP_FBMC./abs(xP_FBMC);
    xP_OFDM = QAM.SymbolMapping(randi(QAM.ModulationOrder,NrPilotSymbols,1));
    xP_OFDM = xP_OFDM./abs(xP_OFDM);
    
    %% Transmitted Data Symbols (Map bin to symbol)
    x_FBMC_Aux = AuxiliaryMethod.PrecodingMatrix*[xP_FBMC;xD_FBMC_Aux];
    x_FBMC_Cod = CodingMethod.PrecodingMatrix*[xP_FBMC;xD_FBMC_Cod];
    x_OFDM     = PilotMapping_OFDM*[xP_OFDM;xD_OFDM];
                
    %% Transmitted Signal (time domain)
    s_FBMC_Aux = G_FBMC*x_FBMC_Aux(:); % Same as "FBMC.Modulation(x_FBMC_Aux)" which is computationally more efficient. But G_FBMC is consistent with the paper.
    s_FBMC_Cod = G_FBMC*x_FBMC_Cod(:); 
    s_OFDM     = G_OFDM*x_OFDM(:);
       
    %% Channel
    ConvolutionMatrix = ChannelModel.GetConvolutionMatrix{1};
   
    r_FBMC_Aux_noNoise = ConvolutionMatrix*s_FBMC_Aux;
    r_FBMC_Cod_noNoise = ConvolutionMatrix*s_FBMC_Cod;
    r_OFDM_noNoise     = ConvolutionMatrix*s_OFDM;
    
    %% Transmission Matrix
    D_FBMC = Q_FBMC'*ConvolutionMatrix*G_FBMC;
    D_OFDM = Q_OFDM'*ConvolutionMatrix*G_OFDM;
        
    %% One-Tap Channel (for perfect channel knowledge)
    h_FBMC = diag(D_FBMC);
    h_OFDM = diag(D_OFDM);
           
    for i_SNR = 1:length(M_SNR_dB)
        %% Add Noise
        SNR_dB  = M_SNR_dB(i_SNR);
        Pn_time = SamplingRate/(F*L)*10^(-SNR_dB/10);
        noise   = sqrt(Pn_time/2)*(randn(size(s_OFDM))+1j*randn(size(s_OFDM)));

        r_FBMC_Aux = r_FBMC_Aux_noNoise+noise; 
        r_FBMC_Cod = r_FBMC_Cod_noNoise+noise;
        r_OFDM     = r_OFDM_noNoise+noise;

        %% Demodulate FBMC signal
        y_FBMC_Aux         = Q_FBMC'*r_FBMC_Aux; % Same as "FBMC.Demodulation(r_FBMC_Aux)" 
        y_FBMC_Cod         = Q_FBMC'*r_FBMC_Cod; % Same as "FBMC.Demodulation(r_FBMC_Cod)" 
        y_FBMC_Cod_PostCod = CodingMethod.PrecodingMatrix'*y_FBMC_Cod;
        y_OFDM             = Q_OFDM'*r_OFDM; % Same as "OFDM.Demodulation(r_OFDM)" 

        %% Channel Estimation at Pilot Position
        hP_est_FBMC_Aux = y_FBMC_Aux(PilotMatrix_FBMC==1)./xP_FBMC/sqrt(Kappa_Aux);
        hP_est_FBMC_Cod = y_FBMC_Cod(PilotMatrix_FBMC==1)./xP_FBMC/sqrt(Kappa_Cod);
        hP_est_OFDM     = y_OFDM(PilotMatrix_OFDM==1)./xP_OFDM/sqrt(Kappa_OFDM);

        %% Estimate Transmit Matix
        D_FBMC_est_Aux = sum(bsxfun(@times,...
            reshape(full(W_MMSE_FBMC_Aux(:,i_SNR)),size(G_FBMC,2),size(G_FBMC,2),NrPilotSymbols),...
            reshape(hP_est_FBMC_Aux,1,1,[])),3);
        D_FBMC_est_Cod = sum(bsxfun(@times,...
            reshape(full(W_MMSE_FBMC_Cod(:,i_SNR)),size(G_FBMC,2),size(G_FBMC,2),NrPilotSymbols),...
            reshape(hP_est_FBMC_Cod,1,1,[])),3);     
        D_OFDM_est = sum(bsxfun(@times,...
            reshape(full(W_MMSE_OFDM(:,i_SNR)),size(G_OFDM,2),size(G_OFDM,2),NrPilotSymbols),...
            reshape(hP_est_OFDM,1,1,[])),3);     

        %% One-Tap Equalizer
        h_est_FBMC_Aux = diag(D_FBMC_est_Aux);  
        x_est_OneTapEqualizer_FBMC_Aux = y_FBMC_Aux./h_est_FBMC_Aux;
        xD_est_OneTapEqualizer_FBMC_Aux = real(x_est_OneTapEqualizer_FBMC_Aux(AuxilaryPilotMatrix_FBMC(:)==0)./sqrt(AuxiliaryMethod.DataPowerReduction));
        DetectedBitStream_OneTapEqualizer_FBMC_Aux = PAM.Symbol2Bit(xD_est_OneTapEqualizer_FBMC_Aux);   
        BER_FBMC_Aux_OneTapEqualizer(i_SNR,i_rep) = mean(BinaryDataStream_FBMC_Aux~=DetectedBitStream_OneTapEqualizer_FBMC_Aux);    
        BER_FBMC_Aux_OneTapEqualizer_NoEdge(i_SNR,i_rep) = mean(BinaryDataStream_FBMC_Aux(ConsideredBits_FBMC_Aux)~=DetectedBitStream_OneTapEqualizer_FBMC_Aux(ConsideredBits_FBMC_Aux));

        h_est_FBMC_Cod = diag(D_FBMC_est_Cod); 
        x_est_OneTapEqualizer_FBMC_Cod = CodingMethod.PrecodingMatrix'*(y_FBMC_Cod./h_est_FBMC_Cod);
        xD_est_OneTapEqualizer_FBMC_Cod = real(x_est_OneTapEqualizer_FBMC_Cod(NrPilotSymbols+1:end))/CodingMethod.DataPowerReduction;
        DetectedBitStream_OneTapEqualizer_FBMC_Cod = PAM.Symbol2Bit(xD_est_OneTapEqualizer_FBMC_Cod);   
        BER_FBMC_Cod_OneTapEqualizer(i_SNR,i_rep) = mean(BinaryDataStream_FBMC_Cod~=DetectedBitStream_OneTapEqualizer_FBMC_Cod);    
        BER_FBMC_Cod_OneTapEqualizer_NoEdge(i_SNR,i_rep) = mean(BinaryDataStream_FBMC_Cod(ConsideredBits_FBMC_Cod)~=DetectedBitStream_OneTapEqualizer_FBMC_Cod(ConsideredBits_FBMC_Cod));

        h_est_OFDM = diag(D_OFDM_est);
        x_est_OneTapEqualizer_OFDM = y_OFDM./h_est_OFDM;
        xD_est_OneTapEqualizer_OFDM = x_est_OneTapEqualizer_OFDM(PilotMatrix_OFDM(:)==0)./sqrt(DataPowerReduction_OFDM);   
        DetectedBitStream_OneTapEqualizer_OFDM = QAM.Symbol2Bit(xD_est_OneTapEqualizer_OFDM);   
        BER_OFDM_OneTapEqualizer(i_SNR,i_rep) = mean(BinaryDataStream_OFDM~=DetectedBitStream_OneTapEqualizer_OFDM);
        BER_OFDM_OneTapEqualizer_NoEdge(i_SNR,i_rep) = mean(BinaryDataStream_OFDM(ConsideredBits_OFDM)~=DetectedBitStream_OneTapEqualizer_OFDM(ConsideredBits_OFDM));

        %% One-Tap Equalizer, Perfect Channel Knowledge
        x_est_OneTapEqualizer_FBMC_Aux_PerfectCSI = y_FBMC_Aux./h_FBMC;
        xD_est_OneTapEqualizer_FBMC_Aux_PerfectCSI = real(x_est_OneTapEqualizer_FBMC_Aux_PerfectCSI(AuxilaryPilotMatrix_FBMC(:)==0)./sqrt(AuxiliaryMethod.DataPowerReduction));
        DetectedBitStream_OneTapEqualizer_FBMC_Aux_PerfectCSI = PAM.Symbol2Bit(xD_est_OneTapEqualizer_FBMC_Aux_PerfectCSI);   
        BER_FBMC_Aux_OneTapEqualizer_PerfectCSI(i_SNR,i_rep) = mean(BinaryDataStream_FBMC_Aux~=DetectedBitStream_OneTapEqualizer_FBMC_Aux_PerfectCSI);    
        BER_FBMC_Aux_OneTapEqualizer_PerfectCSI_NoEdge(i_SNR,i_rep) = mean(BinaryDataStream_FBMC_Aux(ConsideredBits_FBMC_Aux)~=DetectedBitStream_OneTapEqualizer_FBMC_Aux_PerfectCSI(ConsideredBits_FBMC_Aux));    

        x_est_OneTapEqualizer_FBMC_Cod_PerfectCSI = CodingMethod.PrecodingMatrix'*(y_FBMC_Cod./h_FBMC);
        xD_est_OneTapEqualizer_FBMC_Cod_PerfectCSI = real(x_est_OneTapEqualizer_FBMC_Cod_PerfectCSI(NrPilotSymbols+1:end))/CodingMethod.DataPowerReduction;
        DetectedBitStream_OneTapEqualizer_FBMC_Cod_PerfectCSI = PAM.Symbol2Bit(xD_est_OneTapEqualizer_FBMC_Cod_PerfectCSI);   
        BER_FBMC_Cod_OneTapEqualizer_PerfectCSI(i_SNR,i_rep) = mean(BinaryDataStream_FBMC_Cod~=DetectedBitStream_OneTapEqualizer_FBMC_Cod_PerfectCSI);    
        BER_FBMC_Cod_OneTapEqualizer_PerfectCSI_NoEdge(i_SNR,i_rep) = mean(BinaryDataStream_FBMC_Cod(ConsideredBits_FBMC_Cod)~=DetectedBitStream_OneTapEqualizer_FBMC_Cod_PerfectCSI(ConsideredBits_FBMC_Cod));    

        x_est_OneTapEqualizer_OFDM_PerfectCSI = y_OFDM./h_OFDM;
        xD_est_OneTapEqualizer_OFDM_PerfectCSI = x_est_OneTapEqualizer_OFDM_PerfectCSI(PilotMatrix_OFDM(:)==0)./sqrt(DataPowerReduction_OFDM);   
        DetectedBitStream_OneTapEqualizer_OFDM_PerfectCSI = QAM.Symbol2Bit(xD_est_OneTapEqualizer_OFDM_PerfectCSI);   
        BER_OFDM_OneTapEqualizer_PerfectCSI(i_SNR,i_rep) = mean(BinaryDataStream_OFDM~=DetectedBitStream_OneTapEqualizer_OFDM_PerfectCSI);
        BER_OFDM_OneTapEqualizer_PerfectCSI_NoEdge(i_SNR,i_rep) = mean(BinaryDataStream_OFDM(ConsideredBits_OFDM)~=DetectedBitStream_OneTapEqualizer_OFDM_PerfectCSI(ConsideredBits_OFDM));

        %% Improved Channel Estimation and Data Detection
        xD_est_FBMC_Aux_Temp = xD_est_OneTapEqualizer_FBMC_Aux; % initialize with one tap estimates    
        xD_est_FBMC_Cod_Temp = xD_est_OneTapEqualizer_FBMC_Cod; % initialize with one tap estimates    
        xD_est_OFDM_Temp     = xD_est_OneTapEqualizer_OFDM; % initialize with one tap estimates
        xD_est_FBMC_Aux_PerfectCSI_Temp = xD_est_OneTapEqualizer_FBMC_Aux_PerfectCSI; % initialize with one tap estimates    
        xD_est_FBMC_Cod_PerfectCSI_Temp = xD_est_OneTapEqualizer_FBMC_Cod_PerfectCSI; % initialize with one tap estimates    
        xD_est_OFDM_PerfectCSI_Temp     = xD_est_OneTapEqualizer_OFDM_PerfectCSI; % initialize with one tap estimates     
        D_FBMC_est_Aux_Temp  = D_FBMC_est_Aux;
        D_FBMC_est_Cod_Temp  = D_FBMC_est_Cod;
        D_OFDM_est_Temp      = D_OFDM_est;
        h_est_FBMC_Aux_Temp  = h_est_FBMC_Aux;
        h_est_FBMC_Cod_Temp  = h_est_FBMC_Cod;
        h_est_OFDM_Temp      = h_est_OFDM;
        for i_iteration = 1:NrIterations
            y_FBMC_Aux_InterferenceCancellation = (y_FBMC_Aux(:) - (D_FBMC_est_Aux_Temp-diag(h_est_FBMC_Aux_Temp))*AuxiliaryMethod.PrecodingMatrix*[xP_FBMC;PAM.SymbolQuantization(xD_est_FBMC_Aux_Temp)]);                
            y_FBMC_Cod_InterferenceCancellation = (y_FBMC_Cod - (D_FBMC_est_Cod_Temp-diag(h_est_FBMC_Cod_Temp))*CodingMethod.PrecodingMatrix*[xP_FBMC;PAM.SymbolQuantization(xD_est_FBMC_Cod_Temp)]);                   
            y_OFDM_InterferenceCancellation     = (y_OFDM - (D_OFDM_est_Temp-diag(h_est_OFDM_Temp))*PilotMapping_OFDM*[xP_OFDM;QAM.SymbolQuantization(xD_est_OFDM_Temp)]);        

            % New Channel Estimates at Pilot Positions            
            hP_est_FBMC_Aux_Temp = y_FBMC_Aux_InterferenceCancellation(PilotMatrix_FBMC==1)./xP_FBMC/sqrt(Kappa_Aux);
            hP_est_FBMC_Cod_Temp = y_FBMC_Cod_InterferenceCancellation(PilotMatrix_FBMC==1)./xP_FBMC/sqrt(Kappa_Cod);
            hP_est_OFDM_Temp     = y_OFDM_InterferenceCancellation(PilotMatrix_OFDM==1)./xP_OFDM/sqrt(Kappa_OFDM);

            % Improved Channel Estimation
            if i_iteration<=NrIterations/2
                D_FBMC_est_Aux_Temp = sum(bsxfun(@times,...
                    reshape(full(W_MMSE_FBMC_Aux(:,i_SNR)),size(G_FBMC,2),size(G_FBMC,2),NrPilotSymbols),...
                    reshape(hP_est_FBMC_Aux_Temp,1,1,[])),3);
                D_FBMC_est_Cod_Temp = sum(bsxfun(@times,...
                    reshape(full(W_MMSE_FBMC_Cod(:,i_SNR)),size(G_FBMC,2),size(G_FBMC,2),NrPilotSymbols),...
                    reshape(hP_est_FBMC_Cod_Temp,1,1,[])),3);     
                D_OFDM_est_Temp = sum(bsxfun(@times,...
                    reshape(full(W_MMSE_OFDM(:,i_SNR)),size(G_OFDM,2),size(G_OFDM,2),NrPilotSymbols),...
                    reshape(hP_est_OFDM_Temp,1,1,[])),3);          
            else
                D_FBMC_est_Aux_Temp = sum(bsxfun(@times,...
                    reshape(full(W_MMSE_noInterference_FBMC_Aux(:,i_SNR)),size(G_FBMC,2),size(G_FBMC,2),NrPilotSymbols),...
                    reshape(hP_est_FBMC_Aux_Temp,1,1,[])),3);
                D_FBMC_est_Cod_Temp = sum(bsxfun(@times,...
                    reshape(full(W_MMSE_noInterference_FBMC_Cod(:,i_SNR)),size(G_FBMC,2),size(G_FBMC,2),NrPilotSymbols),...
                    reshape(hP_est_FBMC_Cod_Temp,1,1,[])),3);     
                D_OFDM_est_Temp = sum(bsxfun(@times,...
                    reshape(full(W_MMSE_noInterference_OFDM(:,i_SNR)),size(G_OFDM,2),size(G_OFDM,2),NrPilotSymbols),...
                    reshape(hP_est_OFDM_Temp,1,1,[])),3);       
            end

            % One-Tap Channel
            h_est_FBMC_Aux_Temp = diag(D_FBMC_est_Aux_Temp);
            h_est_FBMC_Cod_Temp = diag(D_FBMC_est_Cod_Temp);
            h_est_OFDM_Temp     = diag(D_OFDM_est_Temp);

            x_est_FBMC_Aux_Temp = y_FBMC_Aux_InterferenceCancellation(:)./h_est_FBMC_Aux_Temp;
            x_est_FBMC_Cod_Temp = CodingMethod.PrecodingMatrix'*(y_FBMC_Cod_InterferenceCancellation./h_est_FBMC_Cod_Temp);
            x_est_OFDM_Temp     =  y_OFDM_InterferenceCancellation./h_est_OFDM_Temp;

            xD_est_FBMC_Aux_Temp = real(x_est_FBMC_Aux_Temp(AuxilaryPilotMatrix_FBMC(:)==0)./sqrt(AuxiliaryMethod.DataPowerReduction));
            xD_est_FBMC_Cod_Temp = real(x_est_FBMC_Cod_Temp(NrPilotSymbols+1:end))/CodingMethod.DataPowerReduction;
            xD_est_OFDM_Temp     = x_est_OFDM_Temp(PilotMatrix_OFDM(:)==0)./sqrt(DataPowerReduction_OFDM); 

            DetectedBitStream_FBMC_Aux_Temp = PAM.Symbol2Bit(xD_est_FBMC_Aux_Temp); 
            DetectedBitStream_FBMC_Cod_Temp = PAM.Symbol2Bit(xD_est_FBMC_Cod_Temp); 
            DetectedBitStream_OFDM_Temp     = QAM.Symbol2Bit(xD_est_OFDM_Temp);

            BER_FBMC_Aux_InterferenceCancellation(i_SNR,i_rep,i_iteration) = mean(BinaryDataStream_FBMC_Aux~=DetectedBitStream_FBMC_Aux_Temp);
            BER_FBMC_Cod_InterferenceCancellation(i_SNR,i_rep,i_iteration) = mean(BinaryDataStream_FBMC_Cod~=DetectedBitStream_FBMC_Cod_Temp);   
            BER_OFDM_InterferenceCancellation(i_SNR,i_rep,i_iteration)     = mean(BinaryDataStream_OFDM~=DetectedBitStream_OFDM_Temp);

            BER_FBMC_Aux_InterferenceCancellation_NoEdge(i_SNR,i_rep,i_iteration) = mean(BinaryDataStream_FBMC_Aux(ConsideredBits_FBMC_Aux)~=DetectedBitStream_FBMC_Aux_Temp(ConsideredBits_FBMC_Aux));
            BER_FBMC_Cod_InterferenceCancellation_NoEdge(i_SNR,i_rep,i_iteration) = mean(BinaryDataStream_FBMC_Cod(ConsideredBits_FBMC_Cod)~=DetectedBitStream_FBMC_Cod_Temp(ConsideredBits_FBMC_Cod));   
            BER_OFDM_InterferenceCancellation_NoEdge(i_SNR,i_rep,i_iteration)     = mean(BinaryDataStream_OFDM(ConsideredBits_OFDM)~=DetectedBitStream_OFDM_Temp(ConsideredBits_OFDM));

            
            % Perfect Channel Knoweledge
            y_FBMC_Aux_InterferenceCancellation_PerfectCSI = (y_FBMC_Aux(:) - (D_FBMC-diag(h_FBMC))*AuxiliaryMethod.PrecodingMatrix*[xP_FBMC;PAM.SymbolQuantization(xD_est_FBMC_Aux_PerfectCSI_Temp)]);                
            y_FBMC_Cod_InterferenceCancellation_PerfectCSI = (y_FBMC_Cod - (D_FBMC-diag(h_FBMC))*CodingMethod.PrecodingMatrix*[xP_FBMC;PAM.SymbolQuantization(xD_est_FBMC_Cod_PerfectCSI_Temp)]);                   
            y_OFDM_InterferenceCancellation_PerfectCSI     = (y_OFDM - (D_OFDM-diag(h_OFDM))*PilotMapping_OFDM*[xP_OFDM;QAM.SymbolQuantization(xD_est_OFDM_PerfectCSI_Temp)]);        

            x_est_FBMC_Aux_PerfectCSI_Temp = y_FBMC_Aux_InterferenceCancellation_PerfectCSI./h_FBMC;
            xD_est_FBMC_Aux_PerfectCSI_Temp = real(x_est_FBMC_Aux_PerfectCSI_Temp(AuxilaryPilotMatrix_FBMC(:)==0)./sqrt(AuxiliaryMethod.DataPowerReduction));
            DetectedBitStream_FBMC_Aux_PerfectCSI_Temp = PAM.Symbol2Bit(xD_est_FBMC_Aux_PerfectCSI_Temp);   
            BER_FBMC_Aux_PerfectCSI_InterferenceCancellation(i_SNR,i_rep,i_iteration) = mean(BinaryDataStream_FBMC_Aux~=DetectedBitStream_FBMC_Aux_PerfectCSI_Temp);    
            BER_FBMC_Aux_PerfectCSI_InterferenceCancellation_NoEdge(i_SNR,i_rep,i_iteration) = mean(BinaryDataStream_FBMC_Aux(ConsideredBits_FBMC_Aux)~=DetectedBitStream_FBMC_Aux_PerfectCSI_Temp(ConsideredBits_FBMC_Aux));    

            x_est_FBMC_Cod_PerfectCSI_Temp = CodingMethod.PrecodingMatrix'*(y_FBMC_Cod_InterferenceCancellation_PerfectCSI./h_FBMC);
            xD_est_FBMC_Cod_PerfectCSI_Temp = real(x_est_FBMC_Cod_PerfectCSI_Temp(NrPilotSymbols+1:end))/CodingMethod.DataPowerReduction;
            DetectedBitStream_FBMC_Cod_PerfectCSI_Temp = PAM.Symbol2Bit(xD_est_FBMC_Cod_PerfectCSI_Temp);   
            BER_FBMC_Cod_PerfectCSI_InterferenceCancellation(i_SNR,i_rep,i_iteration) = mean(BinaryDataStream_FBMC_Cod~=DetectedBitStream_FBMC_Cod_PerfectCSI_Temp);    
            BER_FBMC_Cod_PerfectCSI_InterferenceCancellation_NoEdge(i_SNR,i_rep,i_iteration) = mean(BinaryDataStream_FBMC_Cod(ConsideredBits_FBMC_Cod)~=DetectedBitStream_FBMC_Cod_PerfectCSI_Temp(ConsideredBits_FBMC_Cod));    

            x_est_OFDM_PerfectCSI_Temp = y_OFDM_InterferenceCancellation_PerfectCSI./h_OFDM;
            xD_est_OFDM_PerfectCSI_Temp = x_est_OFDM_PerfectCSI_Temp(PilotMatrix_OFDM(:)==0)./sqrt(DataPowerReduction_OFDM);   
            DetectedBitStream_OFDM_PerfectCSI_Temp = QAM.Symbol2Bit(xD_est_OFDM_PerfectCSI_Temp);   
            BER_OFDM_PerfectCSI_InterferenceCancellation(i_SNR,i_rep,i_iteration) = mean(BinaryDataStream_OFDM~=DetectedBitStream_OFDM_PerfectCSI_Temp);
            BER_OFDM_PerfectCSI_InterferenceCancellation_NoEdge(i_SNR,i_rep,i_iteration) = mean(BinaryDataStream_OFDM(ConsideredBits_OFDM)~=DetectedBitStream_OFDM_PerfectCSI_Temp(ConsideredBits_OFDM));

        end    

    end
    TimeNeededSoFar = toc;
    disp([int2str(i_rep/NrRepetitions*100) '% Completed! Time Left: ' int2str(TimeNeededSoFar/i_rep*(NrRepetitions-i_rep)/60) 'min, corresponding to approx. '  int2str(TimeNeededSoFar/i_rep*(NrRepetitions-i_rep)/3600) 'hour']);


    %% Plot results
    % OFDM
    figure(2);
    Markersize = 4;
    hold off;
    semilogy(M_SNR_dB_morePoints, BitErrorProbability,'Color',[1 1 1]*0.75);
    hold on;
    semilogy(M_SNR_dB, nanmean(BER_OFDM_PerfectCSI_InterferenceCancellation(:,:,end),2),'-x black','Markersize',Markersize);
    semilogy(M_SNR_dB, nanmean(BER_OFDM_InterferenceCancellation(:,:,end),2),'-s magenta','Markersize',Markersize);
    semilogy(M_SNR_dB, nanmean(BER_OFDM_InterferenceCancellation_NoEdge(:,:,end),2),'-o blue','Markersize',Markersize);
    semilogy(M_SNR_dB, nanmean(BER_OFDM_OneTapEqualizer_PerfectCSI,2),'-x','Color',[1 1 0]*0.7,'Markersize',Markersize);
    semilogy(M_SNR_dB, nanmean(BER_OFDM_OneTapEqualizer,2),'-s red','Markersize',Markersize);
    ylim([10^-2 0.5]);
    title(['OFDM, Realization ' int2str(i_rep) '/'  int2str(NrRepetitions)])
    legend({'Doubly-Flat Theory','Cancellation (Perfect CSI)','Cancellation','Cancellation (no Edges)','One-Tap (Perfect CSI)','One-Tap'});
    ylabel('Bit Error Ratio');
    xlabel('Signal-to-Noise Ratio [dB]');
    
    % FBMC, Auxiliary Method
    figure(3);
    hold off;
    semilogy(M_SNR_dB_morePoints, BitErrorProbability,'Color',[1 1 1]*0.75);
    hold on;
    semilogy(M_SNR_dB, nanmean(BER_FBMC_Aux_PerfectCSI_InterferenceCancellation(:,:,end),2),'-x black','Markersize',Markersize);
    semilogy(M_SNR_dB, nanmean(BER_FBMC_Aux_InterferenceCancellation(:,:,end),2),'-s magenta','Markersize',Markersize);
    semilogy(M_SNR_dB, nanmean(BER_FBMC_Aux_InterferenceCancellation_NoEdge(:,:,end),2),'-o blue','Markersize',Markersize);
    semilogy(M_SNR_dB, nanmean(BER_FBMC_Aux_OneTapEqualizer_PerfectCSI,2),'-x','Color',[1 1 0]*0.7,'Markersize',Markersize);
    semilogy(M_SNR_dB, nanmean(BER_FBMC_Aux_OneTapEqualizer,2),'-s red','Markersize',Markersize);
    ylim([10^-2 0.5]);
    semilogy([PlotIterationStepsSNRdB PlotIterationStepsSNRdB],[10^-2 10^-1],'Color',[1 1 1]*0.5,'Linewidth',1);
    title(['FBMC Auxiliary Symbols, Realization ' int2str(i_rep) '/'  int2str(NrRepetitions)])
    legend({'Doubly-Flat Theory','Cancellation (Perfect CSI)','Cancellation','Cancellation (no Edges)','One-Tap (Perfect CSI)','One-Tap'});
    ylabel('Bit Error Ratio');
    xlabel('Signal-to-Noise Ratio [dB]');
    
    % FBMC, Data Spreading Method
    figure(4);
    hold off;    
    semilogy(M_SNR_dB_morePoints, BitErrorProbability,'Color',[1 1 1]*0.75);
    hold on;
    semilogy(M_SNR_dB, nanmean(BER_FBMC_Cod_PerfectCSI_InterferenceCancellation(:,:,end),2),'-x black','Markersize',Markersize);
    semilogy(M_SNR_dB, nanmean(BER_FBMC_Cod_InterferenceCancellation(:,:,end),2),'-s magenta','Markersize',Markersize);
    semilogy(M_SNR_dB, nanmean(BER_FBMC_Cod_InterferenceCancellation_NoEdge(:,:,end),2),'-o blue','Markersize',Markersize);
    semilogy(M_SNR_dB, nanmean(BER_FBMC_Cod_OneTapEqualizer_PerfectCSI,2),'-x','Color',[1 1 0]*0.7,'Markersize',Markersize);
    semilogy(M_SNR_dB, nanmean(BER_FBMC_Cod_OneTapEqualizer,2),'-s red','Markersize',Markersize);
    ylim([10^-2 0.5]);
    title(['FBMC Data Spreading, Realization ' int2str(i_rep) '/'  int2str(NrRepetitions)])
    legend({'Doubly-Flat Theory','Cancellation (Perfect CSI)','Cancellation','Cancellation (no Edges)','One-Tap (Perfect CSI)','One-Tap'});
    ylabel('Bit Error Ratio');
    xlabel('Signal-to-Noise Ratio [dB]');
    
    % FBMC, Auxiliary Method, BER over Interation
    figure(5);
    hold off;
    semilogy(0:NrIterations, repmat(BitErrorProbability(find(PlotIterationStepsSNRdB==M_SNR_dB_morePoints)),NrIterations+1,1),'Color',[1 1 1]*0.75);
    hold on;
    Index = find(PlotIterationStepsSNRdB==M_SNR_dB);
    semilogy(0:NrIterations, [nanmean(BER_FBMC_Aux_OneTapEqualizer_PerfectCSI(Index,:),2);squeeze(nanmean(BER_FBMC_Aux_PerfectCSI_InterferenceCancellation(Index,:,:),2))],'-x black','Markersize',Markersize);
    semilogy(0:NrIterations, [nanmean(BER_FBMC_Aux_OneTapEqualizer(Index,:),2);squeeze(nanmean(BER_FBMC_Aux_InterferenceCancellation(Index,:,:),2))],'-s magenta','Markersize',Markersize);
    semilogy(0:NrIterations, [nanmean(BER_FBMC_Aux_OneTapEqualizer_NoEdge(Index,:),2);squeeze(nanmean(BER_FBMC_Aux_InterferenceCancellation_NoEdge(Index,:,:),2))],'-o blue','Markersize',Markersize);
    semilogy(0:NrIterations, repmat(nanmean(BER_FBMC_Aux_OneTapEqualizer_PerfectCSI(Index,:),2),NrIterations+1,1),'-x','Color',[1 1 0]*0.7,'Markersize',Markersize);
    semilogy(0:NrIterations, repmat(nanmean(BER_FBMC_Aux_OneTapEqualizer(Index,:),2),NrIterations+1,1),'-s red','Markersize',Markersize);
    title(['FBMC Auxiliary Symbols, Realization ' int2str(i_rep) '/'  int2str(NrRepetitions)])
    legend({'Doubly-Flat Theory','Cancellation (Perfect CSI)','Cancellation','Cancellation (no Edges)','One-Tap (Perfect CSI)','One-Tap'});
    set(gca, 'XTick',0:NrIterations);
    ylabel('Bit Error Ratio');
    xlabel('Iteration Step i');
    
    pause(0.01);

end


%% Plot Additional Information
fprintf('=============================\n');
fprintf('========= Data Rate =========\n');
fprintf('OFDM       |%7.2f Mbit/s  | \n', length(BinaryDataStream_OFDM)     / (OFDM.PHY.TimeSpacing*OFDM.Nr.MCSymbols)/1e6   );
fprintf('FBMC, Aux. |%7.2f Mbit/s  | \n', length(BinaryDataStream_FBMC_Aux) / (OFDM.PHY.TimeSpacing*OFDM.Nr.MCSymbols)/1e6   );
fprintf('FBMC, Cod. |%7.2f Mbit/s  | \n', length(BinaryDataStream_FBMC_Cod) / (OFDM.PHY.TimeSpacing*OFDM.Nr.MCSymbols)/1e6   );
fprintf('=============================\n');

% The power is normalized so that the average transmit power is one
fprintf('================================================\n');
fprintf('============== Relative SNR Shift ==============\n');
fprintf('================================================\n');
fprintf('           |    SNR    |  Data SNR  | Pilot SNR |\n');
fprintf('OFDM       |   %2.1fdB   |   %2.1fdB   |   %2.1fdB   |\n', 0, 10*log10(DataPowerReduction_OFDM), 10*log10(Kappa_OFDM))
fprintf('FBMC, Aux. |   %2.1fdB   |   %2.1fdB   |   %2.1fdB   |\n', 0, 10*log10(AuxiliaryMethod.DataPowerReduction), 10*log10(Kappa_Aux/2))
fprintf('FBMC, Cod. |   %2.1fdB   |   %2.1fdB   |   %2.1fdB   |\n', 0, 10*log10(CodingMethod.DataPowerReduction), 10*log10(Kappa_Cod/2))
fprintf('================================================\n');


