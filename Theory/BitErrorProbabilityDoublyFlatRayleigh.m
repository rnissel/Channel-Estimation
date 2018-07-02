% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2017 by Institute of Telecommunications, TU Wien
% www.tc.tuwien.ac.at    

% This function calculates the bit error probability for an arbitrary
% signal constellation in a doubly flat rayleigh channel.
% It is based on "OFDM and FBMC-OQAM in doubly-selective channels:
% Calculating the bit error probability", R. Nissel and M. Rupp, IEEE
% Communications Letters, 2017
function BitErrorProbability = BitErrorProbabilityDoublyFlatRayleigh(...
    SNR_dB, ...       % The SNR in the complex domain => SNR_FBMC = SNR_OFDM-3dB.
    SymbolMapping, ...% The symbol mapping with mean(SymbolMapping.*conj(SymbolMapping))==1. For example in 4QAM we have: SymbolMapping=[0.7071 + 0.7071i;-0.7071 + 0.7071i;0.7071 - 0.7071i;-0.7071 - 0.7071i];
    BitMapping)       % The bitmapping corresponding to the symbol mapping. For example for 4QAM we have: BitMapping = [0 0;1 0;0 1;1 1];



    % For the decision regions we assume a rectangular regular grid! The rest
    % of the function could also be used for an irregular grid but
    % the decision regions would have to be rewritten!        
    HalfDecisionInterval = min(abs(real(SymbolMapping)));
    DecisionRegions = [...
        real(SymbolMapping)- HalfDecisionInterval ...
        real(SymbolMapping)+ HalfDecisionInterval ...
        imag(SymbolMapping)- HalfDecisionInterval ...
        imag(SymbolMapping)+ HalfDecisionInterval  ];
    DecisionRegions(min(real(SymbolMapping))==real(SymbolMapping),1) = -inf;
    DecisionRegions(max(real(SymbolMapping))==real(SymbolMapping),2) = +inf;
    DecisionRegions(min(imag(SymbolMapping))==imag(SymbolMapping),3) = -inf;
    DecisionRegions(max(imag(SymbolMapping))==imag(SymbolMapping),4) = +inf;
    
    BitErrorProbability = nan(length(SNR_dB),1);
    for i_SNR = 1:length(SNR_dB)
        Pn = 10^(-SNR_dB(i_SNR)/10);
        ProbabilityMatrix = nan(size(SymbolMapping,1),size(SymbolMapping,1));
        for i_symbol = 1:size(SymbolMapping,1)
            x = SymbolMapping(i_symbol);
            Ey2 = abs(x).^2+Pn;
            Eh2 = 1;
            Eyh = x;
            ProbabilityMatrix(:,i_symbol)=GaussianRatioProbabilityRectangularRegion(Ey2,Eh2,Eyh,DecisionRegions(:,1),DecisionRegions(:,2),DecisionRegions(:,3),DecisionRegions(:,4));
        end 
        ErrorProbability = nan(2,size(BitMapping,2));
        for i_bit= 1:size(BitMapping,2)
            for i_zero_one = [0 1]
                index_x = (BitMapping(:,i_bit)==i_zero_one);
                ErrorProbability(i_zero_one+1,i_bit) = mean(sum(ProbabilityMatrix(not(index_x),index_x)));
            end   
        end
        BitErrorProbability(i_SNR) = mean(mean(ErrorProbability));
    end
end


% This function calculates the Probability that the complex Gaussian ratio
% y/h is within the rectrangular region "(zRlower zIlower] x (zRupper
% zIupper]". It requires the function "GaussianRatioCDF"
function Probability = GaussianRatioProbabilityRectangularRegion(...
    Ey2,...             % Expectation{|y|^2}
    Eh2,...             % Expectation{|h|^2}
    Eyh,...             % Expectation{y*conj(h)}
    zRlower,...         % Determines the rectangular region
    zRupper,...
    zIlower,...
    zIupper)


    CDF_RegionA = GaussianRatioCDF(Ey2,Eh2,Eyh,zRupper,zIupper);
    CDF_RegionB = GaussianRatioCDF(Ey2,Eh2,Eyh,zRlower,zIlower);
    CDF_RegionC = GaussianRatioCDF(Ey2,Eh2,Eyh,zRlower,zIupper);
    CDF_RegionD = GaussianRatioCDF(Ey2,Eh2,Eyh,zRupper,zIlower);
    
    Probability = CDF_RegionA+CDF_RegionB-CDF_RegionC-CDF_RegionD;
    
end



% This function calculates the CDF of the complex Gaussian ratio y/h, that
% is Pr(real(y/h)<zR & imag(y/h)<zR), whereas y and h are complex Gaussian random variables
function CDF = GaussianRatioCDF(...
    Ey2,...             % Expectation{|y|^2}
    Eh2,...             % Expectation{|h|^2}
    Eyh,...             % Expectation{y*conj(h)}
    zR,...              % Real part of the CDF, i.e, Pr(real(y/h)<zR &...)
    zI)                 % Imaginary part of the CDF, i.e, Pr(... & imag(y/h)<zR)

a = Eyh/Eh2; % alpha
b = Ey2/Eh2; % beta

Index0      = (zR == -inf) | (zI == -inf);
Index1      = (zR == inf) & (zI == inf);
IndexReal   = (zI == inf) & isfinite(zR);
IndexImag   = (zR == inf) & isfinite(zI);
IndexNormal = isfinite(zR) & isfinite(zI);

% See "Bit error probability for pilot-symbol aided channel estimation in
% FBMC-OQAM", R. Nissel and M. Rupp, 2016
CDF_Real = 1/2-...
    (real(a)-zR(IndexReal))./...
    (2*sqrt((real(a)-zR(IndexReal)).^2+b-abs(a).^2));

CDF_Imag = 1/2-...
    (imag(a)-zI(IndexImag))./...
    (2*sqrt((imag(a)-zI(IndexImag)).^2+b-abs(a).^2));


% General Case, see "OFDM and FBMC-OQAM in Doubly-Selective Channels:
% Calculating the Bit Error Probability" R. Nissel and M. Rupp, 2017
CDF_Normal = 1/4+...
(zR(IndexNormal)-real(a)).*...
 (2*atan(...
            (zI(IndexNormal)-imag(a))./sqrt((zR(IndexNormal)-real(a)).^2+b-abs(a).^2)...
         )+pi)./...
 (4*pi*sqrt((zR(IndexNormal)-real(a)).^2+b-abs(a).^2))+...
(zI(IndexNormal)-imag(a)).*...
 (2*atan(...
            (zR(IndexNormal)-real(a))./sqrt((zI(IndexNormal)-imag(a)).^2+b-abs(a).^2)...
         )+pi)./...
 (4*pi*sqrt((zI(IndexNormal)-imag(a)).^2+b-abs(a).^2));


% Map CDF to correct one
CDF              = nan(size(zR));
CDF(Index0)      = 0;
CDF(Index1)      = 1;
CDF(IndexReal)   = CDF_Real;
CDF(IndexImag)   = CDF_Imag;
CDF(IndexNormal) = CDF_Normal;

end