Segmentation;

%&&&&&&& Initialize Weight Values (1st Network)
% Network Parameters
I  = 65;  % Input  Neurons + 1
ON = 4;   % Output Neurons
H  = 7;   % Hidden Neurons + 1
HN = H-1;  % Hidden Neurons
cc = 90;  % No of cycles
pp = 120;  % No of Patterns
l = 0.05;  % Learning Constant

% Random Initial Weights 
v = random('Normal', 0, 0.1, HN, I); % Hidden
w = random('Normal', 0, 0.1, ON, H); % Real

minErV = 0; % Initial Minimum Validation Cycle Error

% Target Output Vectors
d1  = [1 -1 -1 -1] ;
d2  = [-1 1 -1 -1] ;
d3  = [-1 -1 1 -1] ;
d4  = [-1 -1 -1 1] ;

h1 = 0;
h2 = 30;

for sm=1:4:120
% Recycling output matrices for each pattern

    d{sm}   = d1; 
    d{sm+1} = d2; 
    d{sm+2} = d3;
    d{sm+3} = d4;
    
% Recycling patterns 1-50 for training

    h1 = h1 + 1;
    
    x{sm}   = cube(:, h1);
    x{sm+1} = letter(:, h1);
    x{sm+2} = math(:, h1);
    x{sm+3} = nothing(:, h1);
    
end
 
% Recycling patterns 51-100 for validation
for sm=1:4:120
 
    h2 = h2 + 1;
      
    xx{sm}     = cube(:, h2);
    xx{sm+1}   = letter(:, h2);
    xx{sm+2}   = math(:, h2);
    xx{sm+3}   = nothing(:, h2);
end    

% Cycle testing
for c=1:cc 
    ec  = 0 ;
    ecv = 0 ;
    
% p : No of patterns
    for p=1:pp    
    x{p}(I)  = -1; % Augmentation of Input
    xx{p}(I) = -1;
    
        % Calculate Hidden Outputs
        for j=1:HN
            net   = abs(v(j,:) * x{p});
            y(j)  = (1-exp(-net))/(1+exp(-net));
        
            net   = abs(v(j,:) * xx{p});
            yy(j) = (1-exp(-net))/(1+exp(-net));    
        end
    
    % Augmentation of Hidden Outputs
    y(H)  = -1 ; % 65th row
    yy(H) = -1 ;
    
    % Initial Pattern Error
    EP  = 0 ;
    EPV = 0 ;
        
    for k=1:ON
        
        net = w(k,:) * y' ;
        % Output vectors
        z(k) = (1-exp(-net))/(1+exp(-net));
        
        % Delta learning rule
        outputDelta(k) = 0.5*((d{p}(k)-z(k))*(1-(z(k)^2)));

        % Pattern error calculation
        EP = EP + 0.5 * (d{p}(k) - z(k))^2;
        
        net   = w(k,:) * yy' ;
        zz(k) = (1-exp(-net))/(1+exp(-net));
        EPV   = EPV + 0.5 * (d{p}(k) - zz(k))^2;
        
        % Update weights based on error
        for cm=1:H
            w2(k,cm) = w(k,cm) + l * outputDelta(k)* y(cm);
        end      
    end
    
% Hidden Delta Learning Rule
    for j=1:HN
        
        f = 0.5 * (1-y(j)^2);
        net = outputDelta * w(:,j) ;
        hiddenDelta(j)= net * f;
        
        for cm=1:I
        v2(j,cm) = v(j,cm) + l * hiddenDelta(j) * x{p}(cm);
        end   
    end
    if sum(isnan(v2(:))) == 0
        v = v2 ;
    end
    if sum(isnan(w2(:))) == 0
        w = w2 ;
    end
    
    % Calculate Cycle Error
    ec  = ec  + EP  ;        
    ecv = ecv + EPV ;
    
    end
    
    EC(c)  = ec  ;
    ECV(c) = ecv ;
    
    if c == 1
        minErV = ecv ;
        finalCycleV = 1 ;
        vv = v;
        ww = w;
    % Update outputs based on minimum error
    elseif ecv < minErV 
        minErV = ecv ;
        finalCycleV = c ;
        vv = v;
        ww = w;
    end
end 
    
% Plotting the cycle error curves

minErT = min(EC) ;
finalCycleT = find(EC == min(EC)) ;
         
disp('Final Training Cycle');
disp(finalCycleT);
disp('Minimum Training Cycle Error');
disp(minErT);
disp('Final Validation Cycle');
disp(finalCycleV);
disp('Minimum Validation Cycle Error');
disp(minErV);  
 
c=1:cc;
 
f1 = figure ;
set(f1,'name','Cycle Error','numbertitle','off');
plot(c,EC,'-',c,ECV,'--') ;
xlabel('Cycles');
ylabel('Cycle Error');
legend('Training','Validation');

% Training Accuracy
 
SmS = 1;     % First Sample (Testing)
SmE = 30;    % Last  Sample (Testing)
    
% Cube
acCtTrain = 0;
        
for sm=SmS:SmE          
    xx    = cube(:,sm); 
    xx(I) = -1;   
    for j=1:HN
        net = vv(j,:) * xx;
        y(j) = (1-exp(-net))/(1+exp(-net));
    end
    
    y(H) = -1;    
    for k=1:ON
        net = ww(k,:) * y';
        z(k) = (1-exp(-net))/(1+exp(-net));
    end

    % Find which element is the largest
    index = find(z == max(z));
    
    if index == 1
        acCtTrain = acCtTrain + 1;  
    end
end
    
% Letter
dcCtTrain = 0;
        
for sm=SmS:SmE         
    x    = letter(:,sm);
    x(I) = -1;   
    for j=1:HN
        net = vv(j,:) * xx;
        y(j) = (1-exp(-net))/(1+exp(-net));
    end
    
    y(H) = -1;        
    for k=1:ON
        net = ww(k,:) * y';
        z(k) = (1-exp(-net))/(1+exp(-net));
    end
    
    index = find(z == max(z));
    
    if index == 2
        dcCtTrain = dcCtTrain + 1;  
    end
end
    
% Math
gCtTrain = 0 ;
        
for sm=SmS:SmE          
    xx    = math(:,sm) ;
    xx(I) = -1 ;   
    for j=1:HN
        net = vv(j,:) * xx ;
        y(j) = (1-exp(-net))/(1+exp(-net));
    end
    
    y(H) = -1;
    for k=1:ON
        net = ww(k,:) * y' ;
        z(k) = (1-exp(-net))/(1+exp(-net));
    end
    
    index = find(z == max(z)) ;
    
    if index == 3
      gCtTrain = gCtTrain + 1 ;  
    end
end
    
% Nothing
battCtTrain = 0 ;
        
for sm=SmS:SmE   
    xx    = nothing(:,sm) ;
    xx(I) = -1 ;
    for j=1:HN
        net = vv(j,:) * xx ;
        y(j) = (1-exp(-net))/(1+exp(-net));
    end
    
    y(H) = -1;  
    
    for k=1:ON
        net = ww(k,:) * y' ;
        z(k) = (1-exp(-net))/(1+exp(-net));
    end
    
    index = find(z == max(z)) ;
    
    if index == 4
      battCtTrain = battCtTrain + 1 ;  
    end
end
     
% Accuracy as percentages

acAyTrain   = (acCtTrain   / 30) * 100;
dcAyTrain   = (dcCtTrain   / 30) * 100;
gAyTrain    = (gCtTrain    / 30) * 100;
battAyTrain = (battCtTrain / 30) * 100;
    
accuracyTrain = (acAyTrain+dcAyTrain+gAyTrain+battAyTrain)/4 ;

disp('Training Accuracy');
disp(accuracyTrain);
 
% Validation Accuracy
 
SmS = 31;     % First Sample (Testing)
SmE = 60;     % Last  Sample (Testing)

acCtValidation = 0 ; % Cube
        
for sm=SmS:SmE 
    xx    = cube(:,sm); 
    xx(I) = -1 ;
    for j=1:HN
        net = vv(j,:) * xx ;
        y(j) = (1-exp(-net))/(1+exp(-net));
    end
    
    y(H) = -1;
        
    for k=1:ON
        net = ww(k,:) * y' ;
        z(k) = (1-exp(-net))/(1+exp(-net));
    end
    
    index = find(z == max(z)) ;
    
    if index == 1
      acCtValidation = acCtValidation + 1 ;  
    end
end
    
dcCtValidation = 0; % Letter
        
for sm=SmS:SmE 
    xx    = letter(:,sm) ;
    xx(I) = -1 ;
    for j=1:HN
        net = vv(j,:) * xx ;
        y(j) = (1-exp(-net))/(1+exp(-net));
    end
    
    y(H) = -1;
        
    for k=1:ON
        net = ww(k,:) * y' ;
        z(k) = (1-exp(-net))/(1+exp(-net));
    end
    
    index = find(z == max(z)) ;
    
    if index == 2
      dcCtValidation = dcCtValidation + 1 ;  
    end
end
    
gCtValidation = 0 ; % Math
        
for sm=SmS:SmE   
    xx    = math(:,sm) ;
    xx(I) = -1 ;
    for j=1:HN
        net = vv(j,:) * xx ;
        y(j) = (1-exp(-net))/(1+exp(-net));
    end
    
    y(H) = -1;
            
    for k=1:ON
        net = ww(k,:) * y' ;
        z(k) = (1-exp(-net))/(1+exp(-net));
    end
    
    index = find(z == max(z)) ;
    
    if index == 3
      gCtValidation = gCtValidation + 1 ;  
    end
    end

battCtValidation = 0 ; % Nothing

for sm=SmS:SmE 
    xx    = nothing(:,sm) ;
    xx(I) = -1 ;
    for j=1:HN
        net = vv(j,:) * xx ;
        y(j) = (1-exp(-net))/(1+exp(-net));
    end
    
    y(H) = -1;
        
    for k=1:ON
        net = ww(k,:) * y' ;
        z(k) = (1-exp(-net))/(1+exp(-net));
    end
    
    index = find(z == max(z)) ;
    
    if index == 4
      battCtValidation = battCtValidation + 1 ;  
    end
end
    
acAyValidation   = (acCtValidation   / 50) * 100;
dcAyValidation   = (dcCtValidation   / 50) * 100;
gAyValidation    = (gCtValidation    / 50) * 100;
battAyValidation = (battCtValidation / 50) * 100;
    
accuracyValidation = (acAyValidation+dcAyValidation+gAyValidation+battAyValidation)/4 ;
    
disp('Validation Accuracy');
disp(accuracyValidation);

% Offline Testing
    
SmS = 61; % First Sample (Testing)
SmE = 90; % Last  Sample (Testing)

acCt = 0;  % Cube
        
for sm=SmS:SmE 
    xx    = cube(:,sm); 
    xx(I) = -1;
    for j=1:HN
        net = vv(j,:) * xx;
        y(j) = (1-exp(-net))/(1+exp(-net));
    end
    
    y(H) = -1;
        
    for k=1:ON
        net = ww(k,:) * y';
        z(k) = (1-exp(-net))/(1+exp(-net));
    end
    
    index = find(z == max(z));
    
    if index == 1
      acCt = acCt + 1;  
    end
end
    
dcCt = 0 ; % Letter
        
for sm=SmS:SmE 
    xx    = letter(:,sm);
    xx(I) = -1;
    for j=1:HN
        net = vv(j,:) * xx;
        y(j) = (1-exp(-net))/(1+exp(-net));
    end
    
    y(H) = -1;
        
    for k=1:ON
        net = ww(k,:) * y';
        z(k) = (1-exp(-net))/(1+exp(-net));
    end
    
    index = find(z == max(z));
    
    if index == 2
      dcCt = dcCt + 1;  
    end
end
    
gCt = 0; % Math
        
for sm=SmS:SmE 
    xx    = math(:,sm);
    xx(I) = -1;
    for j=1:HN
        net = vv(j,:) * xx;
        y(j) = (1-exp(-net))/(1+exp(-net));
    end
    
    y(H) = -1;
            
    for k=1:ON
        net = ww(k,:) * y';
        z(k) = (1-exp(-net))/(1+exp(-net));
    end
    
    index = find(z == max(z));
    
    if index == 3
      gCt = gCt + 1;  
    end
end
    
battCt = 0; % Nothing
        
for sm=SmS:SmE 
    xx    = nothing(:,sm);
    xx(I) = -1;
    for j=1:HN
        net = vv(j,:) * xx;
        y(j) = (1-exp(-net))/(1+exp(-net));
    end
    
    y(H) = -1;
        
    for k=1:ON
        net = ww(k,:) * y';
        z(k) = (1-exp(-net))/(1+exp(-net));
    end
    
    index = find(z == max(z));
    
    if index == 4
      battCt = battCt + 1;  
    end
end
    
acAy   = (acCt   / 30) * 100;
dcAy   = (dcCt   / 30) * 100;
gAy    = (gCt    / 30) * 100;
battAy = (battCt / 30) * 100;
    
accuracy = (acAy+dcAy+gAy+battAy)/4;     
disp('Testing Accuracy');
disp(accuracy);

plot(1:64, FFT), xlim([0 64])

ChBcol=[32 32 64 64];
ChAcol=[0 0 32 32];
totalcol=[0 1200 1200 0];
patch(ChAcol,totalcol,[0 1 0],'FaceAlpha',0.2,'EdgeColor','none')
patch(ChBcol,totalcol,[1 0 0],'FaceAlpha',0.2,'EdgeColor','none')