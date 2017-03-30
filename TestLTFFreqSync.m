clear;clc;close all;
t=linspace(8*10^-6,16*10^-6,8*10^2);
dt=8*10^-6/10^2/8;

itx=1;
MCS=4;
HTLENGTH=24;
NumHTLTF=4;
STBC=0;
state=0;
NSR=28;
TCS=0;
tHT_LTF1=8*10^-6;
[HTLTF1,state]=FieldGnerator('HT-LTF1',t-tHT_LTF1,itx,MCS,NumHTLTF,STBC,HTLENGTH,state);

x=HTLTF1;

[FreqOffset, ~]=LTF_FreqSync(1,x,dt)
[y,~]=LTF_ChanEstimation(1,x,FreqOffset,dt);
plot(y)
