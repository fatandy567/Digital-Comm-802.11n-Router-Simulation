clear;clc;close all;
t=linspace(0,10*10^-6,10^4);
dt=10*10^-6/10^4;

itx=1;
MCS=4;
HTLENGTH=0.25*10^3;
NumHTLTF=4;
STBC=0;
state=0;
[LSTF,state]=FieldGnerator('L-STF',t,itx,MCS,NumHTLTF,STBC,HTLENGTH,state);

pulse=zeros(1,length(t));
for i=1:length(t)
    [pulse(i),~]=LSTFDetection(1,real([zeros(1,10^3) LSTF]),dt,i);
end
plot(t,pulse,t,LSTF)