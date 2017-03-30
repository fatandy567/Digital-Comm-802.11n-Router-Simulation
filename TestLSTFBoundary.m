clear;clc;close all;
t=linspace(6*10^-6,18*10^-6,2*10^3);
dt=12*10^-6/10^3/2;

itx=1;
MCS=4;
HTLENGTH=24;
NumHTLTF=4;
STBC=0;
state=0;
tHT_LTF1=8*10^-6;
tHT_SIG=tHT_LTF1+8*10^-6;
tHT_LTFs=tHT_SIG+8*10^-6;
tHT_DATA=tHT_LTFs+(NumHTLTF-1)*4*10^-6;

[LSTF,state]=FieldGenerator('L-STF',t,itx,MCS,NumHTLTF,STBC,HTLENGTH,state);
[HTLTF1,state]=FieldGenerator('HT-LTF1',t-tHT_LTF1,itx,MCS,NumHTLTF,STBC,HTLENGTH,state);
[HTSIG,state]=FieldGenerator('HT-SIG',t-tHT_SIG,itx,MCS,NumHTLTF,STBC,HTLENGTH,state);

x=LSTF+HTLTF1+HTSIG;
x1=x+[0 0 0 0 0 0 x(1:end-6)];% test with delay signal
pulse1=zeros(1,length(t));
for i=1:length(t)
    [pulse1(i),~]=STF_LTFBoundary(1,x,dt,i);
end
pulse2=zeros(1,length(t));
for i=1:length(t)
    [pulse2(i),~]=STF_LTFBoundary(1,x,dt,i);
end

[V I]=max(pulse1);
pulse1(I-floor(8*10^-6/dt)+1)=1;
[V I]=max(pulse2);
pulse2(I-floor(8*10^-6/dt)+1)=1;
%plot(t,pulse1,t,pulse2,t,real(x1),t,real(x))
plot(t,pulse2,t,real(x))
