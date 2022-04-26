clear ; clc; close all
%% create time series
fre = 400;
t = 1/fre:1/fre:500;
t = berow(t);
N = length(t);
% generate theta rhythm
theta = 8.5;
period = 1/theta;
gamma = [80, 100];
amp1 = 0.2;
phi1 = 0;
amp2 = 0.35;
phi2 = pi/3;
a1 = amp1*sin(2*pi*theta.*t+phi1);
a2 = amp2*sin(2*pi*theta.*t+phi2);

% get the trough of thetas  start of gamma oscillation
thetatrghs1 = (2.3*pi/2-phi1)./(2*pi)*period:period:t(end);
thetatrghs2 = (2.3*pi/2-phi2)./(2*pi)*period:period:t(end);

% get the trough of thetas  end of gamma oscillation
thetatrghe1 = (3.7*pi/2-phi1)./(2*pi)*period:period:t(end);
thetatrghe2 = (3.7*pi/2-phi2)./(2*pi)*period:period:t(end);

% make gamma band amplitude modulated to trough of theta
ampgb1 = zeros(size(a1));
ampgb2 = zeros(size(a2));
sd = 0.012; %std of gaussian bump
for i = 1:length(thetatrghs1)
    % the start and end of gamma would not be precise, add some noise
    thetatrghs1(i) = thetatrghs1(i) + (rand-0.5)*period/20;
    thetatrghe1(i) = thetatrghe1(i) + (rand-0.5)*period/20;
    ttghtemp = (thetatrghs1(i) + thetatrghe1(i))/2;
    ampgb1 = ampgb1 + amp1*0.3*normpdf(t, ttghtemp, sd);
end
for i = 1:length(thetatrghs2)
    % the start and end of gamma would not be precise, add some noise
    thetatrghs2(i) = thetatrghs2(i) + (rand-0.5)*period/20;
    thetatrghe2(i) = thetatrghe2(i) + (rand-0.5)*period/20;
    ttghtemp = (thetatrghs2(i) + thetatrghe2(i))/2;
    ampgb2 = ampgb2 + amp2*0.2*normpdf(t, ttghtemp, sd);
end

% generate gamma with phase aligned to theta trough
gamma1 = zeros(size(a1));
gamma2 = zeros(size(a2));
for i = 1:length(thetatrghs1)
    ttind = t >= thetatrghs1(i) &  t <= thetatrghe1(i);
    % for each theta trough, generate a gamma rhythm, with random frequency
    % in [80,120], random amplitude, and 0 phase with respect to the start
    gaa = 80;
    gab = 100;
    ftemp = (gab-gaa).*rand + gaa;
    ampa = 0.02;
    ampb = 0.05;
    amptemp = (ampb-ampa).*rand + ampa;
    tst = thetatrghs1(i);
    phitemp = -2*pi*ftemp.*tst;
    gatemp = amptemp*sin(2*pi*ftemp.*t+phitemp);
    % add a little white noise
    gatemp = gatemp+(rand(1,length(t))-0.5)*amptemp;
    mgb1temp = berow(gatemp).*berow(ttind);
    gamma1 = gamma1 + mgb1temp;
end

for i = 1:length(thetatrghs2)
    ttind = t >= thetatrghs2(i) &  t <= thetatrghe2(i);
    % for each theta trough, generate a gamma rhythm, with random frequency
    % in [80,120], random amplitude, and 0 phase with respect to the start
    gaa = 80;
    gab = 100;
    ftemp = (gab-gaa).*rand + gaa;
    ampa = 0.03;
    ampb = 0.08;
    amptemp = (ampb-ampa).*rand + ampa;
    tst = thetatrghs2(i);
    phitemp = -2*pi*ftemp.*tst;
    gatemp = amptemp*sin(2*pi*ftemp.*t+phitemp);
    % add a little white noise
    gatemp = gatemp+(rand(1,length(t))-0.5)*amptemp;
    mgb2temp = berow(gatemp).*berow(ttind);
    gamma2 = gamma2 + mgb2temp;
end

% generate gamma rhythm
gamma1 = berow(gamma1).* berow(ampgb1);
gamma2 = berow(gamma2).* berow(ampgb2);

% generate the time series, add weak pink noise
s1 = berow(a1) + gamma1 + berow(pinknoise_my(N,fre,amp1*0.5));
s2 = berow(a2) + gamma2 + berow(pinknoise_my(N,fre,amp2*0.75));

% plot to check the amplitude modulate
tind = t<=1;

figure()
subplot(3,2,1)
hold on
plot(t(tind),a1(tind))
plot(t(tind),gamma1(tind))

subplot(3,2,2)
hold on
plot(t(tind),s1(tind))

subplot(3,2,3)
hold on
plot(t(tind),a2(tind))
plot(t(tind),gamma2(tind))

subplot(3,2,4)
hold on
plot(t(tind),s2(tind))

%% compute coherence between s1 and s2

[s11,ff] = cpsd(s1,s1,hamming(fre),0,fre,fre);
[s22,~] = cpsd(s2,s2,hamming(fre),0,fre,fre);
[s12,~] = cpsd(s1,s2,hamming(fre),0,fre,fre);
cspec = s12./(sqrt(s11).*sqrt(s22));
subplot(3,2,5)
plot(ff,abs(cspec))





