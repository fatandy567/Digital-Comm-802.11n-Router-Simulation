clc;close all;
Starttime=7.2*10^-6; %Packet Start time;
Endtime=40*10^-6;%Packet End time;
TimeReso=10^-8; % Time Resolution;
CarrierFreq=2.4*10^9;

itx=1;
MCS=4;
HTLENGTH=250;
dt=TimeReso;
t=linspace(Starttime,Endtime,(Endtime-Starttime)/TimeReso);

%Generate the packet
[PKT,NSYM,Mod]=PacketGenerator(t,itx,MCS,HTLENGTH,CarrierFreq);
save('Mcs=4_PKT.mat')

%Down Convert to baseband
PKT=PKT.*exp(-j*2*pi*CarrierFreq*1.00001*t);
STF_LTFBoundary_En=1;
DATADecoder_En=0;
LTF_ChanEstimation_En=0;
for i=1:length(t)
    
    [STF_LTFBoundary_Done,~]=STF_LTFBoundary(STF_LTFBoundary_En,PKT,dt,i);
    
    %Channel Estimation---------------------
    if STF_LTFBoundary_Done==true
        LTF_ChanEstimation_En=1;
        LTF_FreqSync_En=0;
    end
    [H,LTF_ChanEstimation_Done,~]=LTF_ChanEstimation(LTF_ChanEstimation_En,PKT,dt,i);
  
    %Decode Data------------------------------------
    if LTF_ChanEstimation_Done==true
        DATADecoder_En=1;
    end
    [y,DATADecoder_Done,nstatus]=DATADecoder(DATADecoder_En,Mod,H,PKT,dt,NSYM,i);
    
    if DATADecoder_Done==true
        break;
    end
    
end