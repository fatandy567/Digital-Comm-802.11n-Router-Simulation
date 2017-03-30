function [y,nstatus]=FreqCorrectionLoop(En,x,t,CoarseFreqOffset,FineFreqOffset,i)
    if En==0
       y=x;
       nstatus=0;
       return 
    end
    FreqOffset=0;
    if CoarseFreqOffset~=0
       FreqOffset=CoarseFreqOffset;
    end
    
    if FineFreqOffset~=0
        FreqOffset=FineFreqOffset;
    end
    y=[x(1:i) x(i+1:end).*exp(j*2*pi*(-FreqOffset)*t(i+1:end))];
    nstatus=1;
end