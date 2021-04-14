clear all; close all
% get library
if not(libisloaded('observers'))
    loadlibrary('observers.so','observers.h')
end
% methods 
%libfunctions observers -full

% read file
data = csvread('link4.csv');
tm   = data(:,1);
q    = data(:,2:7);
qd   = data(:,8:13);
cur  = data(:,14:19);

K = [10.0,10.6956,8.4566,9.0029,9.48,10.1232];
res = libpointer('doublePtr',zeros(1,6));

% use 
% calllib('observers','reset',id)
% to reset observer state

% ========== Momentum Observer ================

% configuration
Kmo = [90,50,50,90,90,40];  

% get observer ID
id_mo = calllib('observers','configMomentumObserver',-1,Kmo);  %    -1 means create and configure

ext_mo = zeros(size(cur));
for i = 2:size(cur,1)
    calllib('observers','getExternalTorque',id_mo,res,q(i,:),qd(i,:),K.*cur(i,:),tm(i)-tm(i-1));
    ext_mo(i,:) = res.Value(:);
end

plot(tm,ext_mo);
figure;

% ======= Disturbance Observer ===============

% configuration
sigma = 21;
xeta = 18;
beta = 50;

% get ID
id_dis = calllib('observers','configDisturbanceObserver',-1,sigma,xeta,beta);

ext_dis = zeros(size(cur));
for i = 2:size(cur,1)
    calllib('observers','getExternalTorque',id_dis,res,q(i,:),qd(i,:),K.*cur(i,:),tm(i)-tm(i-1));
    ext_dis(i,:) = res.Value(:);
end

plot(tm,ext_dis);
title("Disturbance observer");
figure;

% ======= Sliding Mode Observer ============

% configuration
S1 = [20,30,20,30,20,30];
T1 = 2*sqrt(S1);
S2 = [10,10,10,10,10,10];
T2 = 2*sqrt(S2);

% get ID
id_sm = calllib('observers','configSlidingModeObserver',-1,T1,S1,T2,S2);

ext_sm = zeros(size(cur));
for i = 2:size(cur,1)
    calllib('observers','getExternalTorque',id_sm,res,q(i,:),qd(i,:),K.*cur(i,:),tm(i)-tm(i-1));
    ext_sm(i,:) = res.Value(:);
end

plot(tm,ext_sm);
title("Sliding Mode Observer");
figure;

% ========= Kalman Filter Observer =======

% configuration 
S = zeros(6,6);
H = eye(6,6);
Q = blkdiag(0.002*eye(6,6), 0.3*eye(6,6));
R = 0.05*eye(6,6);

% get ID
id_kf = calllib('observers','configDistKalmanObserver',-1,S,H,Q,R);

ext_kf = zeros(size(cur));
for i = 2:size(cur,1)
    calllib('observers','getExternalTorque',id_kf,res,q(i,:),qd(i,:),K.*cur(i,:),tm(i)-tm(i-1));
    ext_kf(i,:) = res.Value(:);
end

plot(tm,ext_kf);
title("Kalman Filter Observer");
figure;

% ========= Kalman Filter Continous System Observer ===== 

% configuration 
%S = zeros(6,6);
%H = eye(6,6);
%Q = blkdiag(0.2*eye(6,6), 30*eye(6,6));
%R = 0.0005*eye(6,6);

% get ID
%id_kfe = calllib('observers','configDistKalmanObserverExp',-1,S,H,Q,R);

%ext_kfe = zeros(size(cur));

%for i = 2:size(cur,1)
%    calllib('observers','getExternalTorque',id_kfe,res,q(i,:),qd(i,:),K.*cur(i,:),tm(i)-tm(i-1));
%    ext_kfe(i,:) = res.Value(:); 
%end

%plot(tm,ext_kfe);
%title("Continous Kalman Filter Observer");
%figure;

% ====== Filtered Dynamics ==============

% cinfiguration 
cutOff = 8;  % rad/s
timeStep = 0.01;  % s

% get ID
id_df = calllib('observers','configFilterDynObserver',-1,cutOff,timeStep);

ext_df = zeros(size(cur));
for i = 2:size(cur,1)
    calllib('observers','getExternalTorque',id_df,res,q(i,:),qd(i,:),K.*cur(i,:),tm(i)-tm(i-1));
    ext_df(i,:) = res.Value(:);
end

plot(tm,ext_df);
title("Filtered dynamics");
figure;

% ===== Filtered Range ===================

% cinfiguration 
%cutOff = 8;  % rad/s
%timeStep = 0.01;  % s
delta = 0.1;

% get ID
id_dr1 = calllib('observers','configFilterRangeObserver',-1,cutOff,timeStep, delta);
id_dr2 = calllib('observers','configFilterRangeObserver',-1,cutOff,timeStep,-delta);

ext_up = zeros(size(cur));
ext_low = zeros(size(cur));
for i = 2:size(cur,1)
    calllib('observers','getExternalTorque',id_dr1,res,q(i,:),qd(i,:),K.*cur(i,:),tm(i)-tm(i-1));
    ext_up(i,:) = res.Value(:);
    calllib('observers','getExternalTorque',id_dr2,res,q(i,:),qd(i,:),K.*cur(i,:),tm(i)-tm(i-1));
    ext_low(i,:) = res.Value(:);
end

plot(tm,[ext_up(:,3),ext_low(:,3),ext_df(:,3)]); 
title("Range");

% =============== exit ===================

calllib('observers','freeAll'); % clear memory

unloadlibrary('observers');