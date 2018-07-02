classdef PilotSymbolAidedChannelEstimation < handle 
   % =====================================================================        
   % This MATLAB class represents an implementation of pilot-aided channel
   % estimation, that is, it allows to estimate the channel by
   % interpolation/extrapolation.
   % The following interpolation methods are supported:
   %    1) 'linear','nearest','natural' based on the MATLAB built-in
   %        function "scatteredInterpolant"
   %    2) 'FullAverage' averages over all pilots => assumes a double flat
   %        channel
   %    3) 'MovingBlockAverage' averages over a few close pilots
   % The pilot pattern can be 'Diamond', 'Rectangular', 'Custom'
   % =====================================================================    
   % Ronald Nissel, rnissel@nt.tuwien.ac.at
   % (c) 2017 by Institute of Telecommunications, TU Wien
   % www.nt.tuwien.ac.at
   % =====================================================================   


   properties (SetAccess = private)
       NrPilotSymbols
       PilotPattern
       PilotSpacingFrequency
       PilotSpacingTime
       InterpolationMethod
       Implementation
       InterpolationProperties
       PilotMatrix
   end
   
   methods
      % Class constructor, define default values. 
      function obj = PilotSymbolAidedChannelEstimation(varargin)
          % Initialize parameters, set default values  
          obj.PilotPattern = varargin{1};                                   
          obj.InterpolationMethod = varargin{3};
 
          % Generate pilot matrix according to the specified pilot pattern.
          % A zero corresponse to a data symbol, a one to a pilot symbol
          switch obj.PilotPattern
              case 'Rectangular'
                  NrSubcarriers = varargin{2}(1,1);
                  obj.PilotSpacingFrequency = varargin{2}(1,2);
                  NrMCSymbols = varargin{2}(2,1);
                  obj.PilotSpacingTime = varargin{2}(2,2);
                  
                  obj.PilotMatrix = zeros(NrSubcarriers,NrMCSymbols);     
                  obj.PilotMatrix(round(mod(NrSubcarriers-1,obj.PilotSpacingFrequency)/2)+1:obj.PilotSpacingFrequency:NrSubcarriers,round(round(mod(NrMCSymbols-1,obj.PilotSpacingTime)/2)+1:obj.PilotSpacingTime:NrMCSymbols)) = true;             
              case 'Diamond'
                  NrSubcarriers = varargin{2}(1,1);
                  obj.PilotSpacingFrequency = varargin{2}(1,2);
                  NrMCSymbols = varargin{2}(2,1);
                  obj.PilotSpacingTime = varargin{2}(2,2);
                  
                  obj.PilotMatrix = zeros(NrSubcarriers,NrMCSymbols);    
                  % There should be a much smarter way of doing this ... but too lazy
                  FrequencyPositionShift = floor((NrSubcarriers-max([(1:2*obj.PilotSpacingFrequency:NrSubcarriers),(1+1/2*obj.PilotSpacingFrequency:2*obj.PilotSpacingFrequency:NrSubcarriers),(1+obj.PilotSpacingFrequency:2*obj.PilotSpacingFrequency:NrSubcarriers),(1+3/2*obj.PilotSpacingFrequency:2*obj.PilotSpacingFrequency:NrSubcarriers)]))/2)+1;
                  TimePositionShift = floor((NrMCSymbols-max([(1:2*obj.PilotSpacingTime:NrMCSymbols),(1+obj.PilotSpacingTime):2*obj.PilotSpacingTime:NrMCSymbols]))/2)+1;
                  obj.PilotMatrix(FrequencyPositionShift:2*obj.PilotSpacingFrequency:NrSubcarriers,TimePositionShift:2*obj.PilotSpacingTime:NrMCSymbols) = 1;
                  obj.PilotMatrix(FrequencyPositionShift+round(1/2*obj.PilotSpacingFrequency):2*obj.PilotSpacingFrequency:NrSubcarriers,round(TimePositionShift+obj.PilotSpacingTime):2*obj.PilotSpacingTime:NrMCSymbols) = 1;
                  obj.PilotMatrix(FrequencyPositionShift+round(obj.PilotSpacingFrequency):2*obj.PilotSpacingFrequency:NrSubcarriers,TimePositionShift:2*obj.PilotSpacingTime:NrMCSymbols) = 1;
                  obj.PilotMatrix(FrequencyPositionShift+round(3/2*obj.PilotSpacingFrequency):2*obj.PilotSpacingFrequency:NrSubcarriers,round(TimePositionShift+obj.PilotSpacingTime):2*obj.PilotSpacingTime:NrMCSymbols) = 1;
              case 'Custom'
                  obj.PilotSpacingFrequency = nan;
                  obj.PilotSpacingTime =nan;
                  
                  obj.PilotMatrix = varargin{2};
              otherwise
                  error('Pilot pattern is not supported! Chose Rectangular Diamond or Custom');
          end          
          obj.NrPilotSymbols = sum(obj.PilotMatrix(:));
          
          % preinitialize interpolation method
          switch obj.InterpolationMethod  
              case {'linear','nearest','natural'}
                  [x_pilot_pos,y_pilot_pos] = find(obj.PilotMatrix);
                  obj.InterpolationProperties = scatteredInterpolant(x_pilot_pos,y_pilot_pos,zeros(obj.NrPilotSymbols,1),obj.InterpolationMethod);                  
              case 'MovingBlockAverage'                             
                  PilotIndices = find(obj.PilotMatrix);
                  PilotMatrixPosNumbered = zeros(size(obj.PilotMatrix));
                  PilotMatrixPosNumbered(PilotIndices)=1:numel(PilotIndices);
                  
                  BlockLengthFrequency = varargin{4}(1);
                  BlockLengthTime = varargin{4}(2);
                 
                  maxF = size(obj.PilotMatrix,1);
                  maxT = size(obj.PilotMatrix,2);
                
                  IndexF = -BlockLengthFrequency:BlockLengthFrequency;
                  IndexT = -BlockLengthTime:BlockLengthTime;
    
                  obj.InterpolationProperties.InterpolationMatrix = zeros(numel(obj.PilotMatrix),obj.NrPilotSymbols);
                  for i_pos = 1:numel(obj.PilotMatrix)
                      Impulse = zeros(size(obj.PilotMatrix));
                      Impulse(i_pos)=1;
                      [posF,posT]=find(Impulse);
                      
                      IndexPosT = posT+IndexT;
                      IndexPosT(IndexPosT<1)=[];
                      IndexPosT(IndexPosT>maxT)=[];
                      IndexPosF = posF+IndexF;
                      IndexPosF(IndexPosF<1)=[];
                      IndexPosF(IndexPosF>maxF)=[];                     
                      Impulse(IndexPosF,IndexPosT)=1;
                      
                      InterpolationMatrixPosPilots = PilotMatrixPosNumbered(logical(Impulse) & logical(obj.PilotMatrix));
                      
                      obj.InterpolationProperties.InterpolationMatrix(i_pos,InterpolationMatrixPosPilots) = 1/numel(InterpolationMatrixPosPilots);
                  end                 
              case 'MMSE'
                  error('Needs to be implemented');
          end              
      end 
      
      function InterpolatedChannel = ChannelInterpolation(varargin)
          obj = varargin{1};
          LSChannelEstimatesAtPilotPosition = varargin{2};
              
          switch obj.InterpolationMethod  
              case {'linear','nearest','natural'}
                  obj.InterpolationProperties.Values = LSChannelEstimatesAtPilotPosition;
                  [yq,xq] = meshgrid(1:size(obj.PilotMatrix,2),1:size(obj.PilotMatrix,1));
                  InterpolatedChannel = obj.InterpolationProperties(xq,yq);
              case 'FullAverage'
                  InterpolatedChannel = ones(size(obj.PilotMatrix))*mean(LSChannelEstimatesAtPilotPosition);
              case 'MovingBlockAverage'
                  InterpolatedChannel = obj.InterpolationProperties.InterpolationMatrix*LSChannelEstimatesAtPilotPosition;       
              case 'MMSE'
                  error('Needs to be done');
              otherwise
                  error('Interpolation method not implemented');
          end
      end
      
      
      
      function AuxiliaryMatrix = GetAuxiliaryMatrix(varargin)
          obj = varargin{1};
          NrAxuiliarySymbols = varargin{2};
          
          AuxiliaryMatrix = obj.PilotMatrix;
          [index_l,index_k]=find(obj.PilotMatrix);  
          if (min(index_l)<2) || (max(index_l)>=size(obj.PilotMatrix,1))
              warning('Pilots should not be close to the border! There might be a problem!'); 
          elseif (min(index_k)<2) || (max(index_k)>=size(obj.PilotMatrix,2))
              warning('Pilots should not be close to the border! There might be a problem!'); 
          end
              
          for i_lk = 1:size(index_l,1)
              switch NrAxuiliarySymbols
                  case 1                 
                     AuxiliaryMatrix(index_l(i_lk),index_k(i_lk)+1) = -1;           
                  case 2
                     AuxiliaryMatrix(index_l(i_lk),index_k(i_lk)+1) = -1;     
                     AuxiliaryMatrix(index_l(i_lk),index_k(i_lk)-1) = -1;     
                  case 3
                     AuxiliaryMatrix(index_l(i_lk),index_k(i_lk)+1) = -1;
                     AuxiliaryMatrix(index_l(i_lk),index_k(i_lk)-1) = -1;    
                     AuxiliaryMatrix(index_l(i_lk)+1,index_k(i_lk)) = -1;    
                  case 4
                     AuxiliaryMatrix(index_l(i_lk),index_k(i_lk)+1) = -1;
                     AuxiliaryMatrix(index_l(i_lk),index_k(i_lk)-1) = -1;    
                     AuxiliaryMatrix(index_l(i_lk)+1,index_k(i_lk)) = -1;                          
                     AuxiliaryMatrix(index_l(i_lk)-1,index_k(i_lk)) = -1;                            
                  otherwise
                     error('Only 1,2,3,4 auxiliary symbols per pilot are supported');
              end
          end   
      end
      
      function InterpolationMatrix = GetInterpolationMatrix(varargin)
          obj = varargin{1};
          
          [x_pilot_pos,y_pilot_pos] = find(obj.PilotMatrix);
          InterpolationMatrix = zeros(numel(obj.PilotMatrix),numel(x_pilot_pos));
          for i_pos =1:length(x_pilot_pos)
            TestDirac = zeros(size(x_pilot_pos));
            TestDirac(i_pos)=1;
            ImpulseResponse = obj.ChannelInterpolation(TestDirac);
            
            InterpolationMatrix(:,i_pos)=ImpulseResponse(:);
          end
          
      end
      
      function PlotPilotPattern(varargin)
          if numel(varargin)==2
            PilotMatrixTemp = varargin{2};
          else
            PilotMatrixTemp = varargin{1}.PilotMatrix;
          end
          
          PilotMatrixRGB(:,:,1) = PilotMatrixTemp==0;
          PilotMatrixRGB(:,:,2) = not(PilotMatrixTemp==1);    
          PilotMatrixRGB(:,:,3) = not(PilotMatrixTemp==-1);          
          
          imagesc(PilotMatrixRGB); 
          hold on;
          for i_row = 1:size(PilotMatrixRGB,1)+1
             plot([.5,size(PilotMatrixRGB,2)+0.5],[i_row-.5,i_row-.5],'k-');
          end
          for i_column = 1:size(PilotMatrixRGB,2)+1
             plot([i_column-.5,i_column-.5],[.5,size(PilotMatrixRGB,1)+0.5],'k-');
          end
          xlabel('Time index');
          ylabel('Frequency index');
          title('Pilot pattern');
      end
  
   end


end