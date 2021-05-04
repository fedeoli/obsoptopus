close all
%clear all

global ObserverTest
ObserverTest.EulerAngleNoiseOnMag = 0;
ObserverTest.AttitudeZeroErrors=1;
ObserverTest.GyroGaussianCovariance = [1 1 1]*1E-4;

%[MagTemp,GyrosTemp] = AttitudeObserver_GetMeasures_v2_1(Agent(k).iner_ECI(1:3,i+1)',Agent(k).attitude(1:4,i+1)',...
%Agent(k).attitude(5:7,i+1)',Agent(k).magnetoBias(:,i+1)',Agent(k).gyroBias(:,i+1)',0); %????
inerEciAgent = [0, 0,  6.997005529470000e+03 ];                     
eul = [0, pi/2, 0]; 
attitudeAgent = eul2quat(eul);%[    6.505062931568687e-01,     7.465385161643023e-01,     -1.053399846902226e-01,     9.178940065157556e-02];
omegaAgent = [0,     0,     1.059147735401656e-03   ];
magnetoBias  = [0,0,0];     
gyroBias  = [0,0,0];


[MagTemp1,GyrosTemp1] = AttitudeObserver_GetMeasures_v2_1(inerEciAgent,attitudeAgent,...
   omegaAgent,magnetoBias,gyroBias,0);

 myutc = [2019 12 15 10 20 36]; %CHANGE THIS...??                     
 LatLongAlt = eci2lla(inerEciAgent*1E3,myutc); %converto from ECI to latitude, longitude,  altitude
 LatLongAlt(3) = LatLongAlt(3) - 300000;
 inerEciAgent = lla2eci(LatLongAlt,myutc)*1E-3;
[MagTemp,GyrosTemp] = AttitudeObserver_GetMeasures_v2_1(inerEciAgent,attitudeAgent,...
   omegaAgent,magnetoBias,gyroBias,0);

%diffNorm = MagTemp/norm(MagTemp)-MagTemp1/norm(MagTemp1) %already normalized
diff = MagTemp-MagTemp1
 %eul = [0, -pi/2, 0]; 
%attitudeAgent = eul2quat(eul);%[    6.505062931568687e-01,     7.465385161643023e-01,     -1.053399846902226e-01,     9.178940065157556e-02];


%%


%Namplitude = [1E-2 1E-1 1 1E1];
Namplitude = [1E-2 1E-1 1 1E1]*pi/180;
Ntrials = 100;
resultsMag = zeros(3,Ntrials);
resultsGyros = zeros(3,Ntrials);
meanMag = zeros(length(Namplitude),3);
meanGyros = zeros(length(Namplitude),3);
stdMag = zeros(length(Namplitude),3);
stdGyros = zeros(length(Namplitude),3);

[r,p,y] = quat2angle(attitudeAgent);

for j=1:length(Namplitude)
    for i=1:Ntrials
       % [MagTemp,GyrosTemp] = AttitudeObserver_GetMeasures_v2_1(inerEciAgent + randn(1,3)*Namplitude(j),...
        %    attitudeAgent,omegaAgent,magnetoBias,gyroBias,0);
         
        temp = [r,p,y]'+randn(3,1)*Namplitude(j);
        attitudeAgent = quatnormalize(angle2quat(temp(1),temp(2),temp(3)));
        [MagTemp,GyrosTemp] = AttitudeObserver_GetMeasures_v2_1(inerEciAgent ,...
            attitudeAgent,omegaAgent,magnetoBias,gyroBias,0);
        resultsMag(:,i) = MagTemp-MagTemp1;
        resultsGyros(:,i) = norm(GyrosTemp-GyrosTemp1);
    end
    meanMag(j,:) = mean(resultsMag');
    meanGyros(j,:) = mean(resultsGyros');
    stdMag(j,:) = std(resultsMag');
    stdGyros(j,:) = std(resultsGyros');
end
%%
figure(1)
subplot(3,1,1)
errorbar(meanMag(:,1),stdMag(:,1));
subplot(3,1,2)
errorbar(meanMag(:,2),stdMag(:,2));
subplot(3,1,3)
errorbar(meanMag(:,3),stdMag(:,3));

figure(2)
subplot(3,1,1)
errorbar(meanGyros(:,1),stdGyros(:,1));
subplot(3,1,2)
errorbar(meanGyros(:,2),stdGyros(:,2));
subplot(3,1,3)
errorbar(meanGyros(:,3),stdGyros(:,3));
