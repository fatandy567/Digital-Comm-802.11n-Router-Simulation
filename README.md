## Introduction

This is a virtual router implemented by MATLAB. For simplicity, several simplifications are made :

1. Only "Greenfield" mode is considered.
2. HT-SIG is ignored. (You still could generate the field with the code, but it is not decoded)
3. Data is randomly generated, coding is ignored.
4. Only SISO is considered.

## Files

1. **System.m**

   The main simulation system, starting with packet generation. You can customize your simulation by adjusting *Starttime*, *Endtime* ,*TimeReso* ...parameters. 

   For example :

   `Starttime=0; %Packet Start time;`
   `Endtime=6*10^-5;%Packet End time;`
   `TimeReso=10^-8; % Time Resolution;`
   `CarrierFreq=2.4*10^9; %Carrier Frequency`

   Most of the blocks in the router is simulated according to time sequence, which means they would not use the data until its time arrives even though all data are sent as input. Exceptions are *AGC* and *DATADecoder* block.

2. **PacketGenerator.m**

   Generator 802.11n packet in time domain. You could uncomment corresponding line in the code to generate *HT-SIG, HT-LTFs,...etc*.

   ![Digital-Comm-802.11n-Router-Simulation](/Img1.jpg)

3. **Channel.m**

   2-ray model and 4-ray model are implemented. Choosing model by adjusting input parameter.

   ![Digital-Comm-802.11n-Router-Simulation](/Img2.jpg)

4. **AGC.m**

   Automatic gain control.

5. **FreqCorrectionLoop.m**

   Correct carrier frequency offset.

6. **LSTFDetection.m**

   L-STF field arrival detector. If L-STF is detected, output is set to 1, which would turn on other blocks.

7. **STF_FreqSync.m**

   Doing coarse frequency estimation with the repetition of L-STF symbol.

8. **LTF_FreqSync.m**

   Doing fine frequency estimation with the repetition of HT-LTF symbol.

9. **STF_LTFBoundary.m**

   Detect the end of the HT-LTF and turn on DATADecoder block. (HT-SIG is ignored)

10. **DATADecoder.m**

    Decode the data. Real data bits, the decode result and error bits are displayed.

    ![Digital-Comm-802.11n-Router-Simulation](/Img3.JPG)

11. **Other Test.m**

    Testbench for blocks.

## More Detail

[802.11 Standard](http://standards.ieee.org/getieee802/download/802.11-2012.pdf)
