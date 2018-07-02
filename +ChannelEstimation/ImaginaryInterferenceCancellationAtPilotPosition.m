classdef ImaginaryInterferenceCancellationAtPilotPosition < handle 
   % =====================================================================        
   % This MATLAB class represents an implementation of two different
   % imaginary interference cancellation methods for pilot aided channel
   % estimation in FBMC.
   %        1) The classical auxiliary symbol. 
   %        2) The coding (data spreading) approach. 
   % The implementation is based on the paper 
   % R. Nissel,et.al. "On Pilot-Symbol Aided Channel Estimation in FBMC-OQAM"
   % and returns the precoding matrix "obj.PrecodingMatrix", which cancels
   % the imaginary interference at the pilot positions. 
   % =====================================================================    
   % Ronald Nissel, rnissel@nt.tuwien.ac.at
   % (c) 2017 by Institute of Telecommunications, TU Wien
   % www.nt.tuwien.ac.at
   % =====================================================================   

   properties (SetAccess = private)
       Method
       PilotMatrix
       PrecodingMatrix
       PostCodingChannelMatrix
       NrDataSymbols
       NrPilotSymbols
       NrAuxiliarySymbols
       NrTransmittedSymbols
       PilotToDataPowerOffset
       AuxiliaryToDataPowerOffset
       DataPowerReduction
       SIR_dB
       ConsideredInterferenceMatrix
   end
   
   
   methods
      % Class constructor, define default values. 
      function obj = ImaginaryInterferenceCancellationAtPilotPosition(varargin)
          % Initialize parameters, set default values
          Method                        = varargin{1};                  % Cancellation method, either 'Coding' or 'Auxiliary'
          PilotMatrix                   = varargin{2};                  % PilotMatrix, 0 = Data, 1 = Pilot, -1 = Auxiliary symbol                
          FBMCMatrix                    = varargin{3};                  % FBMC transmission matrix D, i.e., y = D*x with x transmitted data symbols and y received data symbols (before equalization)
          NrCanceledInterferersPerPilot = varargin{4};                  % Number of neighboring time-frequency positions which are canceled. The higher the number, the lower the SIR but the higher the complexity. For the coding approach, the canceled symbols must not overlap
          PilotToDataPowerOffset        = varargin{5};                  % Pilot to data power offset. 2 guarantees that the SNR is the same at pilot position and at data position => fair comparision.

          
          % Abs Interference Matrix, same as InterferenceMatrix = FBMC.GetInterferenceMatrix
          InterferenceMatrix_11         = abs(reshape(FBMCMatrix(:,1),size(PilotMatrix)));
          InterferenceMatrix_End1       = abs(reshape(FBMCMatrix(:,size(PilotMatrix,1)),size(PilotMatrix)));
          InterferenceMatrix_1End       = abs(reshape(FBMCMatrix(:,numel(PilotMatrix)-size(PilotMatrix,1)+1),size(PilotMatrix)));
          InterferenceMatrix_EndEnd     = abs(reshape(FBMCMatrix(:,numel(PilotMatrix)),size(PilotMatrix)));
          InterferenceMatrix            = [[InterferenceMatrix_EndEnd;InterferenceMatrix_1End(2:end,:)],[InterferenceMatrix_End1(:,2:end);InterferenceMatrix_11(2:end,2:end)]];         
          
          
          switch Method
              case 'Auxiliary'                  
                  NrPilotSymbols        = sum(PilotMatrix(:)==1);
                  NrDataSymbols         =  sum(PilotMatrix(:)==0);
                  NrAuxiliarySymbols    = sum(PilotMatrix(:)==-1);

                  PseudoInvers      = pinv(FBMCMatrix(PilotMatrix(:)==1,PilotMatrix(:)==-1));
                  AuxMatrixPilots   = PseudoInvers*(eye(NrPilotSymbols)-FBMCMatrix(PilotMatrix(:)==1,PilotMatrix(:)==1));
                  AuxMatrixData     = -PseudoInvers*FBMCMatrix(PilotMatrix(:)==1,PilotMatrix(:)==0);


                  AuxiliaryMatrix = zeros(numel(PilotMatrix),numel(PilotMatrix)-NrAuxiliarySymbols);
                  AuxiliaryMatrix(PilotMatrix==-1,1:NrPilotSymbols)     = AuxMatrixPilots;
                  AuxiliaryMatrix(PilotMatrix==-1,NrPilotSymbols+1:end) = AuxMatrixData;
                  AuxiliaryMatrix(PilotMatrix==1,1:NrPilotSymbols)      = eye(NrPilotSymbols)*sqrt(PilotToDataPowerOffset);
                  AuxiliaryMatrix(PilotMatrix==0,NrPilotSymbols+1:end)  = eye(NrDataSymbols);

                  if NrCanceledInterferersPerPilot>0            
                      [SortedInterferenceValues]    = sort(abs(InterferenceMatrix(:)),'descend');
                      ConsideredInterference_temp   = abs(FBMCMatrix(PilotMatrix(:)==1,:))>=SortedInterferenceValues(NrCanceledInterferersPerPilot+1);
                      Temp_PilotNumber(1,1,:)       = -(1:NrPilotSymbols);
                      ConsideredInterference        = sum(bsxfun(@times,reshape(ConsideredInterference_temp',size(PilotMatrix,1),size(PilotMatrix,2),NrPilotSymbols),Temp_PilotNumber),3);
                      ConsideredInterference(PilotMatrix(:)==1) = 1:NrPilotSymbols; 
                      Index_Pilots  = ConsideredInterference;
                      Index_Data    =  ConsideredInterference;
                      Index_Pilots(PilotMatrix(:)<1)=[];
                      Index_Data(PilotMatrix(:)~=0)=[];

                      AuxiliaryMatrix(PilotMatrix(:)==-1,[Index_Pilots,Index_Data]==0)=0;            
                  else
                      ConsideredInterference = 'All';
                  end

                  % Normalize
                  DataPowerReduction    = (numel(PilotMatrix)/sum(sum(AuxiliaryMatrix.*conj(AuxiliaryMatrix),2),1));
                  AuxiliaryMatrix       = AuxiliaryMatrix*sqrt(DataPowerReduction);


                  FBMCMatrix_Temp   = FBMCMatrix(PilotMatrix==1,:)*AuxiliaryMatrix;
                  SIR_dB= nan(NrPilotSymbols,1);
                  for i_pilot = 1: NrPilotSymbols
                      SIR_dB(i_pilot) = 10*log10(abs(FBMCMatrix_Temp(i_pilot,i_pilot)).^2/(sum(abs(FBMCMatrix_Temp(i_pilot,:)).^2)-abs(FBMCMatrix_Temp(i_pilot,i_pilot)).^2));
                  end

                  Power = diag(AuxiliaryMatrix*AuxiliaryMatrix');

                  AuxiliaryToDataPowerOffset = mean(Power(PilotMatrix(:)==-1))./mean(Power(PilotMatrix(:)==0));
                
                  PrecodingMatrix = AuxiliaryMatrix; 
                  obj.PostCodingChannelMatrix = nan;            
              
              
              case 'Coding'                  
                  NrPilotSymbols                = sum(PilotMatrix(:)==1);
                  NrDataSymbols                 =  numel(PilotMatrix) - 2*NrPilotSymbols;
                  NrAuxiliarySymbols            = 0;
                  AuxiliaryToDataPowerOffset    = 0;
                 

                  [SortedInterferenceValues]    = sort(abs(InterferenceMatrix(:)),'descend');
                  ConsideredInterference_temp   = abs(FBMCMatrix(PilotMatrix(:)==1,:))>=SortedInterferenceValues(NrCanceledInterferersPerPilot+1);

                  if sum(sum(ConsideredInterference_temp,1)>1)
                    error('Coding symbols must not overlap: The pilot-spacing is too small!');
                  end

                  Temp_PilotNumber(1,1,:)   = -(1:NrPilotSymbols);
                  ConsideredInterference    = sum(bsxfun(@times,reshape(ConsideredInterference_temp',size(PilotMatrix,1),size(PilotMatrix,2),NrPilotSymbols),Temp_PilotNumber),3);
                  ConsideredInterference(PilotMatrix(:)==1) = 1:NrPilotSymbols;  
                  NrUncodedDataSymbols      = sum(ConsideredInterference(:)==0);

                  CodingMatrix  = zeros(numel(PilotMatrix),numel(PilotMatrix)-NrPilotSymbols);
                  CodingMatrix(PilotMatrix==1,1:NrPilotSymbols) = eye(NrPilotSymbols)*sqrt(PilotToDataPowerOffset);
                  CodingMatrix(ConsideredInterference(:)==0,NrPilotSymbols+(1:NrUncodedDataSymbols)) = eye(NrUncodedDataSymbols);

                  ColumnIndex_outerLoop = NrPilotSymbols+NrUncodedDataSymbols;
                  for i_pilot = 1:NrPilotSymbols
                      Interference = FBMCMatrix(ConsideredInterference(:)==i_pilot,ConsideredInterference(:)==-i_pilot);
                      % round to compensate for numerical inaccuracies
                      Interference = round(imag(Interference)*10^10)/10^10; % 

                      NrCanceledInterferersPerPilot_temp = length(Interference);

                      [AbsInterferenceSorted,InterferenceSorted_index] = sort(abs(Interference),'descend');
                      InterferenceSorted = Interference(InterferenceSorted_index);          

                      [NumberUniqueInterference,UniqueInterference]=hist(abs(InterferenceSorted),unique(abs(InterferenceSorted)));

                      CodingMatrixOnePilot = zeros(NrCanceledInterferersPerPilot_temp,NrCanceledInterferersPerPilot_temp-1);
                      ColumnIndex = 0;
                      for i_uniqueInterference = 1:length(UniqueInterference)
                          NrElementsCluster = NumberUniqueInterference(i_uniqueInterference);

                          IndexInterference = (AbsInterferenceSorted==UniqueInterference(i_uniqueInterference));
                          InterferenceTemp = InterferenceSorted(IndexInterference).';
                          if mod(log2(NrElementsCluster),1)==0 
                              % check if power of two. If so, use hadamard!
                              C_temp = hadamard(NrElementsCluster)./repmat(InterferenceTemp,1,NrElementsCluster);
                              C_temp(:,1)=[];

                              CodingMatrixOnePilot(IndexInterference,ColumnIndex+(1:size(C_temp,2))) = C_temp;
                              ColumnIndex = ColumnIndex+size(C_temp,2);
                          elseif NrElementsCluster>1
                              % not very efficient! However, this case should not happen! A mor efficient implementation would also use the same cluster method as done later.                    
                              C_temp1 = eye(NrElementsCluster,NrElementsCluster-1)./repmat(InterferenceTemp,[1 NrElementsCluster-1]);
                              C_temp2 = circshift(eye(NrElementsCluster,NrElementsCluster-1),[1 0])./repmat(InterferenceTemp,[1 NrElementsCluster-1]);

                              C_temp = C_temp1-C_temp2;

                              CodingMatrixOnePilot(IndexInterference,ColumnIndex+(1:size(C_temp,2))) = C_temp;
                              ColumnIndex = ColumnIndex+size(C_temp,2);                                                
                          end                              
                      end

                      % combine clusters with the the smallest number of elements => computationaly more efficient
                      Clusters = repmat(abs(InterferenceSorted)',[1 length(NumberUniqueInterference)])==repmat(UniqueInterference,[NrCanceledInterferersPerPilot_temp 1]);       
                      for i_Clusters = 1:size(Clusters,2)-1
                          [~,IndexCluster1]=min(sum(Clusters,1));
                          Cluster1 = Clusters(:,IndexCluster1);
                          Clusters(:,IndexCluster1) = [];
                          [~,IndexCluster2]=min(sum(Clusters,1));
                          Cluster2 = Clusters(:,IndexCluster2);
                          Clusters(:,IndexCluster2) = [];
                          CombineClusterIndex = [find(Cluster1,1) find(Cluster2,1)];
                          Clusters = [Clusters Cluster1+Cluster2];

                          ColumnIndex = ColumnIndex+1;
                          CodingMatrixOnePilot(CombineClusterIndex,ColumnIndex) = [1 -1]./InterferenceSorted(CombineClusterIndex);
                      end

                      % Gram-Schmidt
                      CStart = CodingMatrixOnePilot;
                      CGram = CStart(:,1)/sqrt(CStart(:,1)'*CStart(:,1));
                      for i_gram=2:NrCanceledInterferersPerPilot_temp-1    
                        v= CStart(:,i_gram);    
                        CGram(:,i_gram)=v-sum(repmat(v'*CGram,NrCanceledInterferersPerPilot_temp,1).*CGram,2);
                        CGram(:,i_gram)=CGram(:,i_gram)/sqrt(CGram(:,i_gram)'*CGram(:,i_gram));
                      end  
                      CodingMatrixOnePilot_resorted = zeros(size(CGram));
                      CodingMatrixOnePilot_resorted(InterferenceSorted_index,:) = CGram;

                      % Map Code matrix to correct position
                      CodingMatrix(ConsideredInterference(:)==-i_pilot,ColumnIndex_outerLoop+(1:(NrCanceledInterferersPerPilot_temp-1))) = CodingMatrixOnePilot_resorted;       
                      ColumnIndex_outerLoop =  ColumnIndex_outerLoop+NrCanceledInterferersPerPilot_temp-1;
                  end

                  DataPowerReduction    = (numel(PilotMatrix)/sum(sum(CodingMatrix.*conj(CodingMatrix),2),1));
                  CodingMatrix          = CodingMatrix*sqrt(DataPowerReduction);

                  FBMCMatrix_Temp = FBMCMatrix(PilotMatrix==1,:)*CodingMatrix;
                  SIR_dB= nan(NrPilotSymbols,1);
                  for i_pilot = 1: NrPilotSymbols
                    SIR_dB(i_pilot) = 10*log10(abs(FBMCMatrix_Temp(i_pilot,i_pilot)).^2/(sum(abs(FBMCMatrix_Temp(i_pilot,:)).^2)-abs(FBMCMatrix_Temp(i_pilot,i_pilot)).^2));
                  end   
                  
                  PrecodingMatrix               = CodingMatrix; 
                  obj.PostCodingChannelMatrix   = abs(PrecodingMatrix').^2;    
              otherwise
                  error(['Method must be  ''Auxiliary'' or ''Coding''!']);
          end  
          
          % Set Properties
          obj.Method                        = Method;
          obj.PilotMatrix                   = PilotMatrix;
          obj.PrecodingMatrix               = PrecodingMatrix;
          obj.NrDataSymbols                 = NrDataSymbols;
          obj.NrPilotSymbols                = NrPilotSymbols; 
          obj.NrAuxiliarySymbols            = NrAuxiliarySymbols;
          obj.NrTransmittedSymbols          = size(PrecodingMatrix,1);
          obj.PilotToDataPowerOffset        = PilotToDataPowerOffset;
          obj.AuxiliaryToDataPowerOffset    = AuxiliaryToDataPowerOffset;
          obj.DataPowerReduction            = DataPowerReduction;
          obj.SIR_dB                        = SIR_dB;
          obj.ConsideredInterferenceMatrix  = ConsideredInterference;     
          
      end
           
   end 
end

