function [J, J_Index, Feature] = Jvalue(X,Y)
    
    Xfeature = [my_reshape(X.signal) X.var' my_reshape(X.Hist) X.skew' X.formFactor' X.modeFreq' ...
                X.meanFreq' X.medFreq' my_reshape(X.DST) my_reshape(X.DCT)...
                X.alphaEnergy' X.betaEnergy' X.thetaEnergy' X.deltaEnergy' my_reshape(X.STFT)...
                X.CSP1' X.CSP2' X.CSP1_alpha' X.CSP2_alpha' X.CSP1_beta' X.CSP2_beta'...
                X.CSP1_meanFreq' X.CSP2_meanFreq' X.CSP1_medFreq' X.CSP2_medFreq'];

    Yfeature = [my_reshape(Y.signal) Y.var' my_reshape(Y.Hist) Y.skew' Y.formFactor' Y.modeFreq' ...
                Y.meanFreq' Y.medFreq' my_reshape(Y.DST) my_reshape(Y.DCT)...
                Y.alphaEnergy' Y.betaEnergy' Y.thetaEnergy' Y.deltaEnergy' my_reshape(Y.STFT)...
                Y.CSP1' Y.CSP2' Y.CSP1_alpha' Y.CSP2_alpha' Y.CSP1_beta' Y.CSP2_beta'...
                Y.CSP1_meanFreq' Y.CSP2_meanFreq'  Y.CSP1_medFreq' Y.CSP2_medFreq'];        
%     zeros(20,288*64)...%
	Feature = [Xfeature;Yfeature];
            
    mu0 = mean([Xfeature;Yfeature],1);
    mu1 = mean(Xfeature,1);
    mu2 = mean(Yfeature,1);
    sigma1 = var(Xfeature,0,1);
    sigma2 = var(Yfeature,0, 1);
    
    J = (abs(mu0 - mu1).^2 + abs(mu0 - mu2).^2)./(sigma1 + sigma2+eps);
   

    
    
    %J_Index = struct(1,length(Xfeature));
    clear J_Index
    J_Index(1:360*64, :) = repmat('signal       ',360*64,1);
    J_Index(360*64+1: 361*64, :) = repmat('var          ',64,1);
    J_Index(361*64+1: 386*64, :) = repmat('Hist         ',25*64,1);
    J_Index(386*64+1: 387*64, :) = repmat('skew         ',64,1);
    J_Index(387*64+1: 388*64, :) = repmat('formFactor   ',64,1);
    J_Index(388*64+1: 389*64, :) = repmat('modeFreq     ',64,1);
    J_Index(389*64+1: 390*64, :) = repmat('meanFreq     ',64,1);
    J_Index(390*64+1: 391*64, :) = repmat('medFreq      ',64,1);
    J_Index(391*64+1: 751*64, :) = repmat('DST          ',360*64,1);
    J_Index(751*64+1: 1111*64, :) = repmat('DCT          ',360*64,1);
    J_Index(1111*64+1: 1112*64,:) = repmat('alphaEnergy  ',64,1);
    J_Index(1112*64+1: 1113*64,:) = repmat('betaEnergy   ',64,1);
    J_Index(1113*64+1: 1114*64,:) = repmat('thetaEnergy  ',64,1);
    J_Index(1114*64+1: 1115*64,:) = repmat('deltaEnergy  ',64,1);
    J_Index(1115*64+1: 1403*64,:) = repmat('STFT         ',288*64,1);
    J_Index(1403*64+1: 1403*64+360,:) = repmat('CSP1         ',360,1);
    J_Index(1403*64+360+1: 1403*64+720,:) = repmat('CSP2         ',360,1);
    J_Index(1403*64+721,:) = 'CSP1_alpha   ';
    J_Index(1403*64+722,:) = 'CSP2_alpha   ';
    J_Index(1403*64+723,:) = 'CSP1_beta    ';
    J_Index(1403*64+724,:) = 'CSP2_beta    ';
    J_Index(1403*64+725,:) = 'CSP1_meanfreq';
    J_Index(1403*64+726,:) = 'CSP2_meanfreq';
    J_Index(1403*64+727,:) = 'CSP2_medfreq ';
    J_Index(1403*64+728,:) = 'CSP2_medfreq ';
    
   
    
end
