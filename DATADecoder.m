function [y,Done,nstatus]=DATADecoder(status,Mod,H,x,dt,NSYM,i)
y=[];
if status==0
   y=0;
   nstatus=0;
   Done=0;
   return;
end
T=4*10^-6;
pts=T/dt;

wrong=0;
for a=0:min(floor(length(x(i:end))/pts)-1,NSYM-1)
    fprintf('Data Symbol %d',a);
    fprintf('\n')
    
    data=x(i+1+a*pts:i+(a+1)*pts);
    TGI=0.8*10^-6;
    TGIpts=TGI/dt;
    data=data(TGIpts+1:end);
    
    X=fftshift(fft(data));
    lx=length(X);
    r=[X(lx/2-27:lx/2) X(lx/2+1:lx/2+29)];
    
    FreqAmp=r./H;
    FreqAmp=PilotCancelor(FreqAmp);
    y=Demodulator(Mod,FreqAmp);
    
    %---Create ground truth for comparison-----
    rand('state',0);
    QAM16=[-3 -1 1 3];
    DATA=(QAM16(randi(4,[1,52]))+j*QAM16(randi(4,[1,52])))/sqrt(10);
    realy=Demodulator(Mod,DATA);
    
    for q=1:length(y)
       if strcmp(y(1,q),realy(1,q))==false
           wrong=wrong+1;
       end
       if strcmp(y(2,q),realy(2,q))==false
           wrong=wrong+1;
       end
    end
    fprintf('Correct bits: \n')
    fprintf(realy);
    fprintf('\n');
    
    fprintf('Result bits: \n')
    fprintf(y);
    fprintf('\n\n');
end
fprintf('Total error bits: %d',wrong)
fprintf('\n')
Done=1;
nstatus=1;

function y=Demodulator(Mod,FreqAmp)
Realpart=real(FreqAmp);
Imagpart=imag(FreqAmp);
I=[];
Q=[];
N=length(FreqAmp);
if strcmp(Mod,'BPSK')
    Realpart=Realpart*1;
    Imagpart=Imagpart*1;
    for i=1:N
        if Realpart(i)>0
            I=[I '1'];
        else
            I=[I '0'];
        end
    end
elseif strcmp(Mod,'QPSK')
    Realpart=Realpart*sqrt(2);
    Imagpart=Imagpart*sqrt(2);
    for i=1:N
        if Realpart(i)>0
            I=[I '1'];
        else
            I=[I '0'];
        end
        
        if Imagpart(i)>0
            Q=[Q '1'];
        else
            Q=[Q '0'];
        end
        
    end
elseif strcmp(Mod,'16QAM')
    Realpart=Realpart*sqrt(10);
    Imagpart=Imagpart*sqrt(10);
    for i=1:N
        if Realpart(i)<-2
            I=[I '00'];
        elseif Realpart(i)<0 && Realpart(i)>-2
            I=[I '01'];
        elseif Realpart(i)<2 && Realpart(i)>0
            I=[I '11'];
        else
            I=[I '10'];
        end
        
        if Imagpart(i)<-2
            Q=[Q '00'];
        elseif Imagpart(i)<0 && Imagpart(i)>-2
            Q=[Q '01'];
        elseif Imagpart(i)<2 && Imagpart(i)>0
            Q=[Q '11'];
        else
            Q=[Q '10'];
        end
        
    end
    
elseif strcmp(Mod,'64QAM')
    Realpart=Realpart*sqrt(42);
    Imagpart=Imagpart*sqrt(42);
     for i=1:N
        if Realpart(i)<-6
            I=[I '000'];
        elseif Realpart(i)<-4 && Realpart(i)>-6
            I=[I '001'];
        elseif Realpart(i)<-2 && Realpart(i)>-4
            I=[I '011'];
        elseif Realpart(i)<0 && Realpart(i)>-2
            I=[I '010'];
        elseif Realpart(i)<2 && Realpart(i)>0
            I=[I '110'];
        elseif Realpart(i)<4 && Realpart(i)>2
            I=[I '111'];
        elseif Realpart(i)<6 && Realpart(i)>4
            I=[I '101'];
        elseif Realpart(i)>6
            I=[I '100'];
        end
        
         if Imagpart(i)<-6
            Q=[Q '000'];
        elseif Imagpart(i)<-4 && Imagpart(i)>-6
            Q=[Q '001'];
        elseif Imagpart(i)<-2 && Imagpart(i)>-4
            Q=[Q '011'];
        elseif Imagpart(i)<0 && Imagpart(i)>-2
            Q=[Q '010'];
        elseif Imagpart(i)<2 && Imagpart(i)>0
            Q=[Q '110'];
        elseif Imagpart(i)<4 && Imagpart(i)>2
            Q=[Q '111'];
        elseif Imagpart(i)<6 && Imagpart(i)>4
            Q=[Q '101'];
        elseif Imagpart(i)>6
            Q=[Q '100'];
        end
        
     end
end
y=[I;Q];
function y=PilotCancelor(x)
y=x([1:7,9:21,23:28,30:35,37:49,51:57]);

function y=DATA_M(k)
    y=zeros(1,length(k));
    for i=1:length(k)
        if 0<=k(i) && k(i)<=6
            y(i)=k(i)-28;
        elseif 7<=k(i) && k(i)<=19
            y(i)=k(i)-27;
        elseif 20<=k(i) && k(i)<=25
            y(i)=k(i)-26;
        elseif 26<=k(i) && k(i)<=31
            y(i)=k(i)-25;
        elseif 32<=k(i) && k(i)<=44
            y(i)=k(i)-24;
        elseif 45<=k(i) && k(i)<=51
            y(i)=k(i)-23;
        end
    end


