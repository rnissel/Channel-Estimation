classdef SignalConstellation < handle 
    % =====================================================================        
    % This MATLAB class represents a QAM or PAM signal constellation. The 
    % parameters are initialized by the class contructor. 
    % Afterwards we can transform a bit stream into symbols by 
    % ".Bit2Symbol()" and back from symbol to bit by ".Symbol2Bit".
    % Furthermore LLR calculation for an AWGN channal and MIMO is also
    % included. 
    % =====================================================================    
    % Ronald Nissel, rnissel@nt.tuwien.ac.at
    % (c) 2017 by Institute of Telecommunications, TU Wien
    % www.tc.tuwien.ac.at
    % =====================================================================    

  properties (SetAccess = private)
       Method
       ModulationOrder
       BitMapping
       SymbolMapping
       Implementation
  end
  
  methods
  	function obj = SignalConstellation(ModulationOrder,Method)
        % Generates a QAM or PAM signal constellation object with 
        % corresponding Gray-coded bit mapping. The first argument 
        % represents the modulation order and the second argument the
        % method, either 'QAM' or 'PAM'. For example (4,'QAM'), (256,'QAM')
        % or (2,'PAM'), (16,'PAM')
        
        obj.ModulationOrder = ModulationOrder;
        obj.Method          = Method;
        
        % Gray coded bitmapping
        if strcmp( obj.Method,'QAM' )
            BitMappingAtom = [ones(sqrt(obj.ModulationOrder)/2,1);zeros(sqrt(obj.ModulationOrder)/2,1)];
            for i_temp = 2:log2(sqrt(obj.ModulationOrder))
                BinaryTemp = BitMappingAtom(1:2:end,i_temp-1);
                BitMappingAtom(:,i_temp) = [BinaryTemp;BinaryTemp(end:-1:1)];
            end
            IQ = 2*(1:sqrt(obj.ModulationOrder))-sqrt(obj.ModulationOrder)-1;
            [I_rep,Q_rep]=meshgrid(IQ,IQ);
            obj.SymbolMapping = I_rep(:)+1i*Q_rep(:);
            obj.SymbolMapping = obj.SymbolMapping/sqrt(mean(abs(obj.SymbolMapping).^2));
            obj.BitMapping = false(obj.ModulationOrder,log2(obj.ModulationOrder));
            for x_IQ = IQ
                obj.BitMapping(I_rep(:)==x_IQ,2:2:end) = BitMappingAtom;
                obj.BitMapping(Q_rep(:)==x_IQ,1:2:end) = BitMappingAtom;
            end
        elseif strcmp(obj.Method,'PAM')
            obj.BitMapping = [ones(obj.ModulationOrder/2,1);zeros(obj.ModulationOrder/2,1)];
            for i_temp = 2:log2(obj.ModulationOrder)
                BinaryTemp = obj.BitMapping(1:2:end,i_temp-1);
                obj.BitMapping(:,i_temp) = [BinaryTemp;BinaryTemp(end:-1:1)];
            end
            obj.SymbolMapping = (2*(1:obj.ModulationOrder)-obj.ModulationOrder-1).';
            obj.SymbolMapping = obj.SymbolMapping/sqrt(mean(abs(obj.SymbolMapping).^2));
        else
           error('Signal constellation method must be QAM or PAM!');
        end
        
        % Determine the underlying symbol alphabet and the corresponding
        % bit mapping
        [~,SortOrder] = sort(bi2de(obj.BitMapping),'ascend');
        obj.SymbolMapping = obj.SymbolMapping(SortOrder);
        obj.BitMapping = obj.BitMapping(SortOrder,:);
        
        % For the LLR detection, we determine all data symbols which have a
        % bit value of one (zero) at a certain position
        obj.Implementation.DataSymbolsBitvalueOne   = repmat(obj.SymbolMapping,1,log2(obj.ModulationOrder));
        obj.Implementation.DataSymbolsBitvalueOne   = reshape(obj.Implementation.DataSymbolsBitvalueOne(logical(obj.BitMapping)),obj.ModulationOrder/2,[]);
        obj.Implementation.DataSymbolsBitvalueZero  = repmat(obj.SymbolMapping,1,log2(obj.ModulationOrder));
        obj.Implementation.DataSymbolsBitvalueZero  = reshape(obj.Implementation.DataSymbolsBitvalueZero(not(logical(obj.BitMapping))),obj.ModulationOrder/2,[]);       
    end   
    
    function DataSymbols = Bit2Symbol(obj,BinaryStream)
        % Maps a bit stream to the correpsonding symbol alphabet
        tmpSize = size(BinaryStream);
        DataSymbols = obj.SymbolMapping( bi2de(reshape(BinaryStream(:),log2(obj.ModulationOrder),[])')+1 );
        DataSymbols = reshape( DataSymbols, tmpSize(1)/log2(obj.ModulationOrder), tmpSize(2) );
    end

    function EstimatedBitStream = Symbol2Bit(obj,EstimatedDataSymbols)
        % Maps symbols (nearest neighbor detection) to the corresponding
        % bit stream
        EstimatedDataSymbols = EstimatedDataSymbols(:);
        
        [~,b] = min(abs((repmat(EstimatedDataSymbols,1,obj.ModulationOrder)-repmat((obj.SymbolMapping).',size(EstimatedDataSymbols,1),1)).'));
        EstimatedBitStream = obj.BitMapping(b(:),:).';
        EstimatedBitStream = EstimatedBitStream(:);         
    end
    
    function QuantizedDataSymbols = SymbolQuantization(obj,EstimatedDataSymbols)
        % Performs quantization of the received symbols, that is nearest
        % neighbor detection
        EstimatedDataSymbols = EstimatedDataSymbols(:);

        [~,b] =min(abs((repmat(EstimatedDataSymbols,1,obj.ModulationOrder)-repmat((obj.SymbolMapping).',size(EstimatedDataSymbols,1),1)).'));
        QuantizedDataSymbols = obj.SymbolMapping(b(:),:).';
        QuantizedDataSymbols = QuantizedDataSymbols(:);
    end    
    
    function LLR = LLR_AWGN(obj,y,Pn)
        % Calculates the LLR for an AWGN channel, that is, y=x+n with y
        % denoting the received data symbol, x the transmitted data symbol
        % and n the Gaussian distributed noise with power Pn
    
        if numel(Pn)>1
            PnRepeated = reshape(repmat(Pn.',log2(obj.ModulationOrder)*obj.ModulationOrder/2,1),obj.ModulationOrder/2,[]);
        else 
            PnRepeated = Pn;
        end
  
        ReceivedDataSymbolsRepeated     = reshape(repmat(y.',log2(obj.ModulationOrder)*obj.ModulationOrder/2,1),obj.ModulationOrder/2,[]);
        DataSymbolsBitvalueOneRepeated  = repmat(obj.Implementation.DataSymbolsBitvalueOne,1,length(y));
        DataSymbolsBitvalueZeroRepeated = repmat(obj.Implementation.DataSymbolsBitvalueZero,1,length(y));
        
        LLR =   log(sum(exp(-abs(ReceivedDataSymbolsRepeated-DataSymbolsBitvalueOneRepeated).^2./PnRepeated),1)./...
                    sum(exp(-abs(ReceivedDataSymbolsRepeated-DataSymbolsBitvalueZeroRepeated).^2./PnRepeated),1)).';
        LLR(LLR==Inf)=10^10;
        LLR(LLR==-Inf)=-10^10;     
    end
    
    function LLR = LLR_MIMO_ML(obj,y,H,Rn,Precoder)
        % This method is a straighhfoward (high compuational complexity)
        % implementation of a MIMO maximum likelihood LLR calculation (per
        % bit). The model assumes y=Hx+n, whereas "y" represents the received
        % symbols of size y(NrRX,NrPos), "H" the channel matrix of size 
        % H(NrRX,NrTx,NrPos), "Rn" the correlation matrix of the Gaussian 
        % distributed noise with size Rn(NrRX,NrRX,NrPos), and 
        % "Precoder" the precoding matrix of size Precoder(NrTx,NrStreams)

        if not(exist('Precoder','var'))
            Precoder = eye(size(H,2));
        end
   
        NrTxStreams = size(Precoder,2);
        NrTimeFrequencyPositions = size(H,3);
                     
        IndexRef = repmat(1:obj.ModulationOrder,obj.ModulationOrder^(NrTxStreams-1),1);
        IndexRef = IndexRef(:);
        IndexDataSymbolAntennaCombination = zeros(length(IndexRef),NrTxStreams);
        for i=0:NrTxStreams-1
            IndexDataSymbolAntennaCombination(:,i+1) = repmat(IndexRef(1:obj.ModulationOrder^i:end),obj.ModulationOrder^i,1);
        end
        IndexDataSymbolAntennaCombination = IndexDataSymbolAntennaCombination(1:end/2,:)';
                          
        for i_BitPos = 1:size(obj.Implementation.DataSymbolsBitvalueOne,2)
            % Matlab changes dimensions!!!!!
            if NrTxStreams==2
                x_One_Atom = [obj.Implementation.DataSymbolsBitvalueOne(IndexDataSymbolAntennaCombination(1,:),i_BitPos).';obj.SymbolMapping(IndexDataSymbolAntennaCombination(2,:)).'];
                x_Zero_Atom = [obj.Implementation.DataSymbolsBitvalueZero(IndexDataSymbolAntennaCombination(1,:),i_BitPos).';obj.SymbolMapping(IndexDataSymbolAntennaCombination(2,:)).'];
            else
                x_One_Atom = [obj.Implementation.DataSymbolsBitvalueOne(IndexDataSymbolAntennaCombination(1,:),i_BitPos).';obj.SymbolMapping(IndexDataSymbolAntennaCombination(2:end,:))];
                x_Zero_Atom = [obj.Implementation.DataSymbolsBitvalueZero(IndexDataSymbolAntennaCombination(1,:),i_BitPos).';obj.SymbolMapping(IndexDataSymbolAntennaCombination(2:end,:))];            
            end 
            for i_TX_Antenna = 0:NrTxStreams-1                
                x_One(:,:,i_TX_Antenna+1,i_BitPos) = circshift(x_One_Atom,[i_TX_Antenna,0]);
                x_Zero(:,:,i_TX_Antenna+1,i_BitPos) = circshift(x_Zero_Atom,[i_TX_Antenna,0]);  
            end            
        end
        
        [~,b,c,d] = size(x_One);
        a = size(y,1);
        LLR = zeros(c,d,NrTimeFrequencyPositions);
        for i_TimeFrequencyPos = 1:NrTimeFrequencyPositions
            
            Rn_temp = Rn(:,:,i_TimeFrequencyPos);
            InvRn_temp = Rn_temp^-1;
            [S,V] = svd(InvRn_temp);            
            C = (S*sqrt(V))';

            y_temp = C*y(:,i_TimeFrequencyPos);
            H_temp = C*H(:,:,i_TimeFrequencyPos)*Precoder;
                        
            y_temp_repeated = repmat(y_temp,1,b,c,d);            
            Hx_One_temp = reshape(H_temp*x_One(:,:),a,b,c,d);
            Hx_Zero_temp = reshape(H_temp*x_Zero(:,:),a,b,c,d);
               
            LLR(:,:,i_TimeFrequencyPos)=squeeze(log(sum(exp(sum(-abs(y_temp_repeated-Hx_One_temp).^2,1)),2)./sum(exp(sum(-abs(y_temp_repeated-Hx_Zero_temp).^2,1)),2)));  
        end
        LLR = LLR(:,:).';
        LLR(LLR==Inf)=10000;
        LLR(LLR==-Inf)=-10000;   
    end  
    
    function [LLR,x_est,NoiseScaling] = LLR_MIMO_ZF(obj,y,H,Pn,Precoder)
        % This method calculates the LLR for a MIMO system after a ZF
        % equalization. To keep the compuational complexity low, the cross
        % correlation after the ZF equalizer is neglected. 
        % The model assumes y=Hx+n, whereas "y" represents the received
        % symbols of size y(NrRX,NrPos), "H" the channel matrix of size 
        % H(NrRX,NrTx,NrPos), "Pn" the power of the Gaussian noise n 
        % and "Precoder" the precoding matrix of size Precoder(NrTx,NrStreams) 
        
        if not(exist('Precoder','var'))
            Precoder = eye(size(H,2));
        end  
        
        N = size(H,3);
        NrTxStreams = size(Precoder,2);
        NrRxAntennas = size(y,1);
        
        
        x_est = zeros(NrTxStreams,N);
        NoiseScaling = zeros(NrTxStreams,N);        
        for i_n = 1:N
           H_temp = H(:,:,i_n)*Precoder;           
           if NrRxAntennas<=NrTxStreams
                ZF_Equalizer = H_temp'/(H_temp*H_temp');
           else
                ZF_Equalizer = (H_temp'*H_temp)^-1*H_temp';
           end  
%            ZF_Equalizer = pinv(H_temp);
           x_est(:,i_n) = (ZF_Equalizer*y(:,i_n));           
           NoiseScaling(:,i_n) = Pn*sum(ZF_Equalizer.*conj(ZF_Equalizer),2);
        end
        x_est = x_est.';
        NoiseScaling = NoiseScaling.';
        LLR = reshape(obj.LLR_AWGN(reshape(x_est,[],1), reshape(NoiseScaling,[],1)),[],NrTxStreams); 
    end 
    
    
    function [LLR,x_est,NoiseScaling,UnbiasedScaling] = LLR_MIMO_MMSE(obj,y,H,Pn,Precoder)
        % This method calculates the LLR for a MIMO system after an 
        % unbiased MMSE equalization. To keep the compuational complexity
        % low, the cross correlation after the MMSE equalizer is neglected. 
        % The model assumes y=Hx+n, whereas "y" represents the received
        % symbols of size y(NrRX,NrPos), "H" the channel matrix of size 
        % H(NrRX,NrTx,NrPos), "Pn" the power of the Gaussian noise n 
        % and "Precoder" the precoding matrix of size Precoder(NrTx,NrStreams) 
        
        if not(exist('Precoder','var'))
            Precoder = eye(size(H,2));
        end
        
        N = size(H,3);
        NrTxStreams = size(Precoder,2);
        NrRXAntennas = size(H,1); 
        
        x_est = zeros(NrTxStreams,N);
        NoiseScaling = zeros(NrTxStreams,N); 
        UnbiasedScaling = zeros(NrTxStreams,N); 
        for i_n = 1:N
           H_temp = H(:,:,i_n)*Precoder;
           MMSE_Equalizer = H_temp'/(H_temp*H_temp'+Pn*eye(NrRXAntennas));  
           x_est(:,i_n) = (MMSE_Equalizer*y(:,i_n));
           
           Temp = MMSE_Equalizer*H_temp;
           NoiseScaling(:,i_n) = Pn*sum(MMSE_Equalizer.*conj(MMSE_Equalizer),2)+sum(abs(Temp-diag(diag(Temp))).^2,2);
           UnbiasedScaling(:,i_n) = abs(sum(MMSE_Equalizer.*H_temp.',2));
        end
        x_est = x_est.';
        NoiseScaling = NoiseScaling.';
        UnbiasedScaling = UnbiasedScaling.';
        LLR = reshape(obj.LLR_AWGN(reshape(x_est./UnbiasedScaling,[],1), reshape(NoiseScaling./UnbiasedScaling.^2,[],1)),[],NrTxStreams); 
    end   
    
    function LLR = LLR_MIMO_Sphere(obj,y,H,Pn,Precoder)
        % This method calculates the LLR for a MIMO system by employing
        % a SphereDecoder. The MATLAB Communications Toolbox is REQUIRED!
        % The model assumes y=Hx+n, whereas "y" represents the received
        % symbols of size y(NrRX,NrPos), "H" the channel matrix of size 
        % H(NrRX,NrTx,NrPos), "Pn" the power of the Gaussian noise n 
        % and "Precoder" the precoding matrix of size Precoder(NrTx,NrStreams) 
        
        if exist('Precoder','var')
            Heff = nan(size(H,1),size(Precoder,2),size(H,3));
            for i_H = 1:size(H,3)
                Heff(:,:,i_H)=H(:,:,i_H)*Precoder;
            end
        else
            Heff = H;
        end
        
        Heff_permute = permute(Heff,[3 2 1]);
        CommToolboxSphereDecoder = comm.SphereDecoder('Constellation', obj.SymbolMapping,'BitTable', double(obj.BitMapping), 'DecisionType', 'Soft');
        LLR = step(CommToolboxSphereDecoder,y.',Heff_permute)*(1/Pn);       
    end    
  end
end
      