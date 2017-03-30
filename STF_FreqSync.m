function [FreqOffset,finish,nstatus] = STF_FreqSync(CoarseDet_En,PKT,dt,i) 
% input = STF symbols
% fsMHz = sampling frequency in MHz
fsMHz=1/dt*10^-6;
if CoarseDet_En ==0 
    FreqOffset = 0;
    finish = 0;
    nstatus=0;
    return
end
input = PKT(1:i);
%input = PKT;
DelayBuffer = zeros(1,0.8*fsMHz); % 0.8*fsMHz = the number of samples in one repetition = N
op = zeros(size(input));

for n = 1:length(input)
   op(n) = conj(input(n))*DelayBuffer(end); %r*(n) x r(n-N)
   
   % shifting samples in the delay buffer
   DelayBuffer(2:end) = DelayBuffer(1:end-1);
   DelayBuffer(1) = input(n);
end

FreqOffset = -1*angle(op)/(2*pi*0.8*1e-6);
finish = 1;
nstatus=1;
FreqOffset=mean(FreqOffset(FreqOffset>0));