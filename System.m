clear all;clc;close all;
Starttime=0; %Packet Start time;
Endtime=6*10^-5;%Packet End time;
TimeReso=10^-8; % Time Resolution;
CarrierFreq=2.4*10^9; %Carrier Frequency

itx=1;
MCS=4;
HTLENGTH=250;%0.25KB data payload
dt=TimeReso;
t=linspace(Starttime,Endtime,(Endtime-Starttime)/TimeReso);

%Generate the packet
[PKT,NSYM,Mod]=PacketGenerator(t,itx,MCS,HTLENGTH,CarrierFreq);
figure(1)
plot(t,abs(PKT));grid on
title('Packet Time domain')
xlabel('Time')
ylabel('Abs(.)')

%Channel
PKT = Channel(PKT,1,dt); % 1 for 2-ray model & 2 for 4-ray model
PKT = PKT(1:length(t));
figure(3);
plot(abs(PKT));grid on
title('Packet Time domain After Channel')
xlabel('Time')
ylabel('Abs(.)')


%Receiver
AGC_En                = 1; %AGC is on
FreqCorrectionLoop_En = 1;%CarrierDownConvertor Enable
LSTFDet_En            = 1;%LSTFDetection Enable
STF_LTFBoundary_En    = 0;%LSTFDetection disable
LTF_FreqSync_En       = 0;%Frequency offset estimation disable
LTF_ChanEstimation_En = 0;%Channel Estimation disable
STF_FreqSync_En       = 0;%coarse frequeny detection
DATADecoder_En        = 0;%DATADecoder disable

STF_FreqsyncWaittime=-2;
STF_LTFBoundaryWaittime=-2;

CoarseFrequencyOffset=0;
FineFrequencyOffset=0;

%Down Convert to baseband
PKT=PKT.*exp(-j*2*pi*CarrierFreq*t);
pulse=zeros(1,length(t));

for i=1:length(t)
    %AGC
    [PKT,G] = AGC(AGC_En,PKT,0.0027,5);
    AGC_En=0;% Turn off AGC

    %Frequency Correction
    [PKT,~]=FreqCorrectionLoop(FreqCorrectionLoop_En,PKT,t,CoarseFrequencyOffset,FineFrequencyOffset,i);
    
    %Packet Detection-------------
    [iscoming,~]=LSTFDetection(LSTFDet_En,PKT,dt,i);
    %Boundary Detection for HT LTF1------------------------------
    if iscoming==true;%LSTF is detected
        LSTFDet_En=0;
        STF_FreqsyncWaittime=floor(4*10^-6/dt);%After LSTF is detected, we have to wait for 4us to activate Coarse Frequency Detection(STF_Freqsync)
        STF_LTFBoundaryWaittime=floor(5.6*10^-6/dt);%After LSTF is detected, we have to wait for 5.6us to activate boundary detection
    end
    %After LSTF is detected, we have to wait for 4us to activate boundary detection
    if STF_FreqsyncWaittime>=0;
        STF_FreqsyncWaittime=STF_FreqsyncWaittime-1;
    elseif STF_FreqsyncWaittime==-1;
        STF_FreqSync_En=1;%4us is up, activate boundary detection
    else
        STF_FreqSync_En=0;
    end
    
    %Coarse Frequency Offset Detection via LSTF------------------------------
    [CoarseFrequencyOffset,STF_FreqSync_Done,~] = STF_FreqSync(STF_FreqSync_En,PKT,dt,i);
    if STF_FreqSync_Done==true
        fprintf('Coarse Frequency Offset= %d\n',CoarseFrequencyOffset);
        STF_FreqSync_En=0;
        STF_FreqsyncWaittime=-2;
    end
    
    %After LSTF is detected, we have to wait for 5.6us to activate boundary detection
    if STF_LTFBoundaryWaittime>=0;
        STF_LTFBoundaryWaittime=STF_LTFBoundaryWaittime-1;
    elseif STF_LTFBoundaryWaittime==-1;
        STF_LTFBoundary_En=1;%5.6us is up, activate boundary detection
    else
        STF_LTFBoundary_En=0;
    end
    [STF_LTFBoundary_Done,~]=STF_LTFBoundary(STF_LTFBoundary_En,PKT,dt,i);
    pulse(i)=STF_LTFBoundary_Done;
    
    %Frequency Offset Estimation------------------
    if STF_LTFBoundary_Done==true
        LTF_FreqSync_En=1;
        STF_LTFBoundary_En=0;STF_LTFBoundaryWaittime=-2;%disable boundary detection
    end
    [FineFrequencyOffset,LTF_FreqSync_Done,~]=LTF_FreqSync(LTF_FreqSync_En,PKT,dt,i);
    
    %Channel Estimation---------------------
    if LTF_FreqSync_Done==true
        fprintf('Fine Frequency Offset = %d\n',FineFrequencyOffset);
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
        break; %After decoding the data, stop the system
    end
    %{
    Note:
    This is a time-squential system, but for simplicity,once DataDecoder is turned on, it will decode all the symbol, instead
    of waiting for the right time.
    %}
    
    
end
