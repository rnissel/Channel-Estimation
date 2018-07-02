classdef OFDM < handle
    % =====================================================================        
    % This MATLAB class represents an implementation of OFDM. The 
    % modulation parameters are initialized by the class contructor. 
    % The modulation of data symbols x and the demodulation of the
    % received samples r are then performed by the methods ".Modulation(x)"
    % and ".Demodulation(r)".
    % Additionally, there exist some other useful methods, such as,
    % ".PlotPowerSpectralDensity" or ".GetTXMatrix".
    % =====================================================================    
    % Ronald Nissel, rnissel@nt.tuwien.ac.at
    % (c) 2017 by Institute of Telecommunications, TU Wien
    % www.nt.tuwien.ac.at
    % =====================================================================    
     
    properties (SetAccess = private)
        Nr                  % for dimensionless parameters
        PHY                 % for parameters with physical interpretation
        Implementation      % implmentation relevent parameters
    end
    
    methods
        function obj = OFDM(varargin)
            % Initialize parameters, set default values
            if numel(varargin) == 8
                obj.Nr.Subcarriers              = varargin{1};      % Number of subcarriers
                obj.Nr.MCSymbols                = varargin{2};      % Number OFDM symbols in time
                obj.PHY.SubcarrierSpacing       = varargin{3};      % Subcarrier spacing (Hz)
                obj.PHY.SamplingRate            = varargin{4};      % Sampling rate (Samples/s)
                obj.PHY.IntermediateFrequency   = varargin{5};      % Intermediate frequency of the first subcarrier (Hz)
                obj.PHY.TransmitRealSignal      = varargin{6};      % Transmit real valued signal (sampling theorem must be fulfilled!)
                obj.PHY.CyclicPrefixLength      = varargin{7};      % Length of the cyclic prefix (s)
                obj.PHY.ZeroGuardTimeLength     = varargin{8};      % Length of the guard time (s), that is, zeros at the beginning and at the end of the transmission
            elseif numel(varargin) == 0
                % LTE default values (can later be changed using OFDM.Set...)
                obj.Nr.Subcarriers              = 24;
                obj.Nr.MCSymbols                = 14;
                obj.PHY.SubcarrierSpacing       = 15e3;
                obj.PHY.SamplingRate            = 15e3*24*14;
                obj.PHY.IntermediateFrequency   = 0;
                obj.PHY.TransmitRealSignal      = false;
                obj.PHY.CyclicPrefixLength      = 1/(14*15e3);
                obj.PHY.ZeroGuardTimeLength     = 0;
            else
                error('Number of input variables must be either 0 (default values) or 8');
            end
            
            % calculate and set all dependent parameters
            obj.SetDependentParameters();            
            
        end
        
        function SetDependentParameters(obj)
            % method that sets all parameters which are dependent on other
            % parameters
            
            if mod(round(obj.PHY.SamplingRate/(obj.PHY.SubcarrierSpacing)*10^5)/10^5,1)~=0
                obj.PHY.SubcarrierSpacing=obj.PHY.SamplingRate/(round(obj.PHY.SamplingRate/(obj.PHY.SubcarrierSpacing)));
                disp('Sampling rate must be a multiple of the subcarrier spacing!');
                disp(['Therefore, the subcarrier spacing is set to: ' int2str(obj.PHY.SubcarrierSpacing) 'Hz']);
            end
            
            if mod(round(obj.PHY.IntermediateFrequency/obj.PHY.SubcarrierSpacing*10^5)/10^5,1)~=0
                obj.PHY.IntermediateFrequency = round(obj.PHY.IntermediateFrequency/obj.PHY.SubcarrierSpacing)*obj.PHY.SubcarrierSpacing;
                disp('The intermediate frequency must be a multiple of the subcarrier spacing!');
                disp(['Therefore, the intermediate frequency is set to ' int2str(obj.PHY.IntermediateFrequency) 'Hz']);
            end
            
            if (obj.PHY.SamplingRate<obj.Nr.Subcarriers*obj.PHY.SubcarrierSpacing)
                error('Sampling theorem is not fullfilled: sampling rate must be higher than the number of subcarriers times subcarrier spacing');
            end
            
            if abs(mod(round(obj.PHY.CyclicPrefixLength*obj.PHY.SamplingRate*10^5)/10^5,1))~=0
                obj.PHY.CyclicPrefixLength=round(obj.PHY.CyclicPrefixLength*obj.PHY.SamplingRate)/obj.PHY.SamplingRate;
                disp('The length of the cyclic prefix times the sampling rate must be an integer!');
                disp(['Therefore, the cyclic prefix length is set to ' num2str(obj.PHY.CyclicPrefixLength) 's']);
            end
            
            obj.Implementation.CyclicPrefix             = round( obj.PHY.CyclicPrefixLength * obj.PHY.SamplingRate );                                       % number of samples for the CP
            obj.Implementation.ZeroGuardSamples         = round( obj.PHY.ZeroGuardTimeLength * obj.PHY.SamplingRate );                                      % number of samples for the zero guard
            obj.Implementation.TimeSpacing              = round( obj.PHY.SamplingRate / obj.PHY.SubcarrierSpacing ) + obj.Implementation.CyclicPrefix;      % number of samples per sybol including CP
            obj.Implementation.FFTSize                  = round( obj.PHY.SamplingRate / obj.PHY.SubcarrierSpacing );
            obj.Implementation.IntermediateFrequency    = round( obj.PHY.IntermediateFrequency / obj.PHY.SubcarrierSpacing );
            obj.Implementation.NormalizationFactor      = sqrt( obj.PHY.SamplingRate^2 / obj.PHY.SubcarrierSpacing^2 / obj.Nr.Subcarriers);                 % Normalization factor so that power = 1
            obj.PHY.dt                                  = 1 / obj.PHY.SamplingRate;                                                                         % inter sample time
            obj.PHY.TimeSpacing                         = obj.Implementation.TimeSpacing * obj.PHY.dt;                                                      % symbol duration including CP
            obj.Nr.SamplesTotal                         = (obj.Nr.MCSymbols*obj.Implementation.TimeSpacing)+2*obj.Implementation.ZeroGuardSamples;          % total number of samples for all symbols including zero guards
        end
        
        % Set Functions
        function SetNrSubcarriers(varargin) 
            % set the number of subcarriers
            
            obj = varargin{1};
            % set specific property
            obj.Nr.Subcarriers = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function SetNrMCSymbols(varargin)
            % set the number of symbols
            
            obj = varargin{1};
            % set specific property
            obj.Nr.MCSymbols = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function SetSubcarrierSpacing(varargin)
            % set the subcarrier spacing
            
            obj = varargin{1};
            % set specific property
            obj.PHY.SubcarrierSpacing = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function SetSamplingRate(varargin)
            % set the sampling rate
            
            obj = varargin{1};
            % set specific property
            obj.PHY.SamplingRate = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function SetIntermediateFrequency(varargin)
            % set the intermediate frequency
            
            obj = varargin{1};
            % set specific property
            obj.PHY.IntermediateFrequency = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function SetTransmitRealSignal(varargin)
            % set real transmit signal indicator
            
            obj = varargin{1};
            % set specific property
            obj.PHY.TransmitRealSignal = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
          
        % Modulation and Demodulation
        function TransmitSignal = Modulation(obj,DataSymbols)
            % modulate the data symbols. The input argument is a matrix
            % of size "Number of subcarriers" \times "Number of OFDM symbols (time)"
            % which represents the transmit data symbols
            
            DataSymbolsTemp = zeros(obj.Implementation.FFTSize,obj.Nr.MCSymbols);
            DataSymbolsTemp(obj.Implementation.IntermediateFrequency+(1:obj.Nr.Subcarriers),:) = DataSymbols.*obj.Implementation.NormalizationFactor;
            if obj.PHY.TransmitRealSignal
                DataSymbolsTemp = (DataSymbolsTemp+conj(DataSymbolsTemp([1 end:-1:2],:)))/sqrt(2);
            end
            TransmitSignalNoCP = ifft(DataSymbolsTemp);
            TransmitSignal = [zeros(obj.Implementation.ZeroGuardSamples,1);reshape([TransmitSignalNoCP(end-obj.Implementation.CyclicPrefix+1:end,:);TransmitSignalNoCP],obj.Implementation.TimeSpacing*obj.Nr.MCSymbols,1);zeros(obj.Implementation.ZeroGuardSamples,1)];
        end
        
        function ReceivedSymbols = Demodulation(obj, ReceivedSignal)
            % demodulates the received time signal and returns a matrix of 
            % size "Number of subcarriers" \times "Number of OFDM symbols"
            % which represents the received symbols after demodulation but
            % before channel equalization
                        
            ReceivedSignal_reshape = reshape(ReceivedSignal((obj.Implementation.ZeroGuardSamples+1):(end-obj.Implementation.ZeroGuardSamples)),obj.Implementation.TimeSpacing,obj.Nr.MCSymbols);
            ReceivedSymbolsTemp = fft(ReceivedSignal_reshape(obj.Implementation.CyclicPrefix+1:end,:));
            if obj.PHY.TransmitRealSignal
                %             ReceivedSymbolsTemp = (ReceivedSymbolsTemp+conj(ReceivedSymbolsTemp([1 end:-1:2],:)))/sqrt(2);
                ReceivedSymbolsTemp = ReceivedSymbolsTemp*sqrt(2);
            end
            
            ReceivedSymbols = ReceivedSymbolsTemp(obj.Implementation.IntermediateFrequency+(1:obj.Nr.Subcarriers),:)/(obj.Implementation.NormalizationFactor);
        end
        
        % Matrix Description
        function TXMatrix = GetTXMatrix(obj)
            % returns a matrix G so that s=G*x(:) is the same as 
            % s=obj.Modulation(x)
            
            if obj.PHY.TransmitRealSignal
                % Does not work because conjungate complex is not a linear operation for complex symols. For FBMC it workes because real symbols are transmitted.
                error('GetTXMatrix is not supported for PHY.TransmitRealSignal == true!')
            end
            TXMatrix=zeros(obj.Nr.SamplesTotal,obj.Nr.Subcarriers*obj.Nr.MCSymbols);
            TXMatrixTemp=zeros(obj.Nr.SamplesTotal,obj.Nr.Subcarriers);
            x = zeros(obj.Nr.Subcarriers, obj.Nr.MCSymbols);
            for i_l= 1:obj.Nr.Subcarriers;
                x(i_l)=1;
                TXMatrixTemp(:,i_l) = obj.Modulation(x);
                x(i_l)=0;
            end
            for i_k=1:obj.Nr.MCSymbols
                TXMatrix(:,(1:obj.Nr.Subcarriers)+(i_k-1)*obj.Nr.Subcarriers)=circshift(TXMatrixTemp,[(i_k-1)*obj.Implementation.TimeSpacing,0]);
            end
        end
        
        function RXMatrix = GetRXMatrix(obj)
            % returns a matrix Q so that y=Q*r is the same as 
            % y=reshape(obj.Demodulation(r),[],1)     
            
            if obj.PHY.TransmitRealSignal
                obj.PHY.TransmitRealSignal = false;
                RXMatrix = sqrt(2)*obj.GetTXMatrix'*(obj.Nr.Subcarriers*obj.PHY.SubcarrierSpacing/(obj.PHY.SamplingRate));
                obj.PHY.TransmitRealSignal = true;
            else
                RXMatrix = obj.GetTXMatrix'*(obj.Nr.Subcarriers*obj.PHY.SubcarrierSpacing/(obj.PHY.SamplingRate));
            end
            IndexTemp = obj.Implementation.ZeroGuardSamples+bsxfun(@plus, (1:obj.Implementation.CyclicPrefix)', (0:obj.Nr.MCSymbols-1)*obj.Implementation.TimeSpacing);
            RXMatrix(:,IndexTemp(:))=0;
        end
        
        % Plot
        function [TransmitPower,Time] = PlotTransmitPower(obj, Rx)
            % plot the expected transmit power over time. The input 
            % argument represents the correlation of the data symbols. 
            % If no input argument is specified, an identity matrix is 
            % assumed (uncorrelated data) 
            
            if exist('Rx','var')
                [V,D] = eig(Rx);
            else 
                % assume that Rx is an identity matrix, that is,
                % uncorrelated symbols
                V = eye(obj.Nr.Subcarriers*obj.Nr.MCSymbols);
                D = V;               
            end
            D=sqrt(D);
            TransmitPower = zeros(obj.Nr.SamplesTotal,1);
            for i_lk = 1:obj.Nr.Subcarriers*obj.Nr.MCSymbols
                % For OFDM.PHY.TransmitRealValuedSignal == 0, the modulation is not linear => we have to consider 1 and +j.
                TransmitPower = TransmitPower+(abs(obj.Modulation(reshape(V(:,i_lk),obj.Nr.Subcarriers,obj.Nr.MCSymbols))*D(i_lk,i_lk)).^2+abs(obj.Modulation(reshape(1j*V(:,i_lk),obj.Nr.Subcarriers,obj.Nr.MCSymbols))*D(i_lk,i_lk)).^2)/2;
                if mod(i_lk,1000)==0
                    disp([int2str(i_lk/(obj.Nr.Subcarriers*obj.Nr.MCSymbols )*100) '%']);
                end
            end
            Time = (0:length(TransmitPower)-1)*obj.PHY.dt;
            if nargout==0
                plot(Time,TransmitPower);
                ylabel('Transmit Power');
                xlabel('Time(s)');
            end
        end
        
        function [PowerSpectralDensity,Frequency] = PlotPowerSpectralDensity(obj,Rx)
            % plot the power spectral density. The input argument 
            % represents the correlation of the data symbols. If no input
            % argument is specified, an identity matrix is assumed 
            % (uncorrelated data) 
            
            if exist('Rx','var')
                [V,D] = eig(Rx);
            else 
                V = eye(obj.Nr.Subcarriers*obj.Nr.MCSymbols);
                D = V;               
            end
            D=sqrt(D);
            PowerSpectralDensity = zeros(obj.Nr.SamplesTotal,1);
            for i_lk = 1:obj.Nr.Subcarriers*obj.Nr.MCSymbols
                PowerSpectralDensity = PowerSpectralDensity+abs(fft(obj.Modulation(reshape(V(:,i_lk),obj.Nr.Subcarriers,obj.Nr.MCSymbols))*D(i_lk,i_lk))).^2;
                if mod(i_lk,1000)==0
                    disp([int2str(i_lk/(obj.Nr.Subcarriers*obj.Nr.MCSymbols )*100) '%']);
                end
            end
            Frequency = (0:length(PowerSpectralDensity)-1)*1/(length(PowerSpectralDensity)*obj.PHY.dt);
            PowerSpectralDensity=PowerSpectralDensity/length(PowerSpectralDensity)^2/Frequency(2)^2;
            if nargout==0
                plot(Frequency,10*log10(PowerSpectralDensity));
                ylabel('Power Spectral Density (dB/Hz)');
                xlabel('Frequency (Hz)');
            end
        end
        
        function [PowerSpectralDensity,Frequency] = PlotPowerSpectralDensityUncorrelatedData(obj)
            % plot the power spectral density in case of uncorrelated data
            % symbols. Faster than FBMC.PlotPowerSpectralDensity
            
            PowerSpectralDensity = zeros(obj.Nr.SamplesTotal,1);
            for i_lk = 1:obj.Nr.Subcarriers
                V = zeros(obj.Nr.Subcarriers,obj.Nr.MCSymbols);
                V(i_lk,round(obj.Nr.MCSymbols/2))=1;
                PowerSpectralDensity = PowerSpectralDensity+abs(fft(obj.Modulation(V))).^2;
            end
            Frequency = (0:length(PowerSpectralDensity)-1)*1/(length(PowerSpectralDensity)*obj.PHY.dt);
            PowerSpectralDensity=obj.Nr.MCSymbols*PowerSpectralDensity/length(PowerSpectralDensity)^2/Frequency(2)^2;
            if nargout==0
                plot(Frequency,10*log10(PowerSpectralDensity));
                ylabel('Power Spectral Density (dB)');
                xlabel('Frequency (Hz)');
            end
        end
        
        % SIR and Noise power
        function Pn = GetSymbolNoisePower(obj, Pn_time)
            % returns the symbol noise power, that is, the noise power
            % after demodulation. The input argument is the noise power 
            % in the time domain. 
            
            Pn = (Pn_time*(obj.Nr.Subcarriers*obj.PHY.SubcarrierSpacing/(obj.PHY.SamplingRate)));
        end
        

        function [SignalPower,InterferencePower] = GetSignalAndInterferencePowerQAM(...
                obj,...                                 % OFDM object
                VectorizedChannelCorrelationMatrix,...  % Let the received signal be r=H*s with H representing a time-variant convolution matrix. Then "VectorizedChannelCorrelationMatrix" represents the expectation E{{H(:)*H(:)'}. We can obtain such matrix by ChannelModel.GetCorrelationMatrix
                DataSymbolCorrelationMatrix,...         % Correlation matrix of the vectorized data symbols
                TimeSamplesOffset,...                   % Time offset in samples (to improve the SIR)
                SubcarrierPosition,...                  % Subcarrier position for which the SIR is calculated.
                OFDMSymbolPosition...                   % OFDM symbol position in time for which the SIR is calculated.
                )
            % returns the signal and interference power for a
            % doubly-selective channel in case of QAM transmissions

            TXMatrix = obj.GetTXMatrix;
            RXMatrix = obj.GetRXMatrix;
            RXMatrix = [zeros(size(RXMatrix,1),TimeSamplesOffset),RXMatrix(:,1:end-TimeSamplesOffset)]; % Time offset compensation
              
            Index = SubcarrierPosition+(OFDMSymbolPosition-1)*obj.Nr.Subcarriers;

            % TempOld = kron(TXMatrix.',RXMatrix(Index,:))*VectorizedChannelCorrelationMatrix*kron(TXMatrix.',RXMatrix(Index,:))';
            % Much more efficient
            RXVectorRep = kron(sparse(eye(length(RXMatrix(Index,:)))),RXMatrix(Index,:)');
            Temp = RXVectorRep'*VectorizedChannelCorrelationMatrix*RXVectorRep;
            CorrMatrixNoData = TXMatrix.'*Temp*conj(TXMatrix);

            SignalDataSymbolCorrelationMatrix = zeros(size(DataSymbolCorrelationMatrix));
            SignalDataSymbolCorrelationMatrix(Index,Index) = DataSymbolCorrelationMatrix(Index,Index);
            InterferenceDataSymbolCorrelationMatrix = DataSymbolCorrelationMatrix;
            InterferenceDataSymbolCorrelationMatrix(Index,:)=0;
            InterferenceDataSymbolCorrelationMatrix(:,Index)=0;

            SignalPower = abs(sum(sum(CorrMatrixNoData.*SignalDataSymbolCorrelationMatrix)));
            InterferencePower = abs(sum(sum(CorrMatrixNoData.*InterferenceDataSymbolCorrelationMatrix)));
        end
       
        % Additional Functions
        function TimePos = GetTimeIndexMidPos(obj)
            % returns a vector which represents the discete time position
            % of the corresponding FBMC symbol (middle position)
            
            TimePos = obj.Implementation.ZeroGuardSamples+obj.Implementation.CyclicPrefix+round(obj.Implementation.FFTSize/2)+1+ (0:obj.Nr.MCSymbols-1)*obj.Implementation.TimeSpacing;
        end    
    end
end