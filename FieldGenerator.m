%===============================================
%Green Field Mode:
%   LSTF + HT-LTF + HT-SIG + HT-LTFs + Data
%
%
%===============================================

% LSTF(Legacy Short Training Field)
function [y,NSYM,Mod,nstate]=FieldGenerator(field,t,itx,MCS,NumHTLTF,STBC,HTLENGTH,inistate)
    global state
    state=inistate;
    DATA=0;
    NSYM=0;% the number of data symbol,initiaized
    
    iss=itx-1;
    [NCBPS,Mod,NBPSC,Nss,NES]=MCStoParameters(MCS);
    y=zeros(1,length(t));
    df=20*10^6/64;
    TGI=0.8*10^-6;%guard interval
    gamma=1;
    TR=100*10^-9;%windows transition time
    PHT_LTF=[1,-1,1,1;1,1,-1,1;1,1,1,-1;-1,1,1,1];
    
    TCS=[0,0,0,0;0,-400,0,0;0,-400,-200,0;0 -400,-200,-600];
    TCS=TCS(Nss,itx);
    
    if strcmp(field,'L-STF') % LSTF(Legacy Short Training Field)
        fprintf('Generating L-STF....\n');
        S_LSTF=1/sqrt(2)*[0,0,1+j,zeros(1,3),-1-j,zeros(1,3),1+j,...
        zeros(1,3),-1-j,zeros(1,3),-1-j,zeros(1,3),1+j,...
        zeros(1,7),-1-j,zeros(1,3),-1-j,zeros(1,3),1+j,...
        zeros(1,3),1+j,zeros(1,3),1+j,zeros(1,3),1+j,0,0];
    
        NL_STF=12;
        T=8*10^-6;
   
    
        for i=1:length(t)
            y(i)=1/sqrt(NL_STF)*Wt(TR,T,t(i))*(sum(S_LSTF([1:27]).*exp(j*2*pi*[-26:1:0]*df*(t(i)-TCS)))+...
                gamma*sum(S_LSTF([28:53]).*exp(j*2*pi*[1:1:26]*df*(t(i)-TCS))));
            
        end
    elseif strcmp(field,'L-LTF') % LLTF(Legacy Long Training Field)
        fprintf('Generating L-LTF....\n');
        S_LLTF=[1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,0,...
            1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1];
        
        NL_LTF=52;
        TGI=TGI*2;
       
    
        for i=1:length(t)
            y(i)=1/sqrt(NL_LTF)*Wt(TR,T,t(i))*(sum(S_LLTF([1:27]).*exp(j*2*pi*[-26:1:0]*df*(t(i)-TGI2-TCS)))+...
                gamma*sum(S_LLTF([28:53]).*exp(j*2*pi*[1:1:26]*df*(t(i)-TGI2-TCS))));
            
        end
        
        %HT_LTF High Throughput Long Training Field
    elseif (strcmp(field,'HT-LTF1')||strcmp(field,'HT-LTF2')||strcmp(field,'HT-LTF3')||strcmp(field,'HT-LTF4'))
    
        if strcmp(field,'HT-LTF1')
            n=1;
            fprintf('Generating HT-LTF1....\n');
        elseif strcmp(field,'HT-LTF2')
            n=2;
            fprintf('Generating HT-LTF2....\n');
        elseif strcmp(field,'HT-LTF3')
            n=3;
            fprintf('Generating HT-LTF3....\n');
        elseif strcmp(field,'HT-LTF4')
            n=4;
            fprintf('Generating HT-LTF4....\n');
        end
        
        HT_LTF1=[1,1,1,1,-1, -1,1,1,-1,1, -1,1,1,1,1, 1,1,-1,-1,1 ,1,-1,1,-1,1, 1,1,1,0,...
            1,-1,-1,1,1, -1,1,-1,1,-1, -1,-1,-1,-1,1, 1,-1,-1,1,-1, 1,-1,1,1,1, 1,-1,-1];
        %HT_LTF1=zeros(1,57);
        %HT_LTF1(end)=1;
        
        NHT_LTF=56;
        T=8*10^-6;
        NSR=28;
        for i=1:length(t)
            tp1=sum(PHT_LTF(iss+1,n)*HT_LTF1([1:29]).*exp(j*2*pi*[-28:1:0]*df*(t(i)-2*TGI-TCS)));
            tp2=gamma*sum(PHT_LTF(iss+1,n)*HT_LTF1([30:57]).*exp(j*2*pi*[1:1:28]*df*(t(i)-2*TGI-TCS)));
            y(i)=1/sqrt(NHT_LTF)*Wt(TR,T,t(i))*(tp1+tp2);
        end 
        
    elseif strcmp(field,'HT-SIG') %HT-SIG (High Throughput Signal Field)
        fprintf('Generating HT-SIG....\n');
        MCS=de2bi(MCS,7);
        BW2040=de2bi(0,1);
        HTLENGTH=de2bi(HTLENGTH,16);
        Aggre=de2bi(0,1);
        RESERVED=[1 1 1];
        STBC=de2bi(STBC,2);
        AdvCoding=de2bi(0,1);
        ShortGI=1;
        NumHTLTF=de2bi(NumHTLTF-1,2);
        CRC=calc_crc([MCS,BW2040,HTLENGTH,Aggre,RESERVED,STBC,ShortGI,NumHTLTF]);
        Tailbit=de2bi(0,6);
        NBPSC_here=1;%the number of coded bits per subcarrier
        %HT-SIG is transmitted using the most robust data rate (BPSK), 
        %and it encodes information about how the rest of the data will be transmitted
        
        NCBPS_here=96;% the number of bits in a single OFDM symbol
        d=[MCS,BW2040,HTLENGTH,RESERVED,Aggre,STBC,AdvCoding,ShortGI,NumHTLTF,CRC,Tailbit];
        d=ConvCoder(d);%convolutional code
        d=d(interleaver(NCBPS_here,NBPSC_here));%interleave data
        
        d=Mapper(d, NBPSC_here);
        
        Wix=PHT_LTF(itx,1);
        TSYM=4*10^-6;
        NHT_SIG=56;
        T=8*10^-6;
        
        
        Ppilot=[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,...
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0];
        
        for i=1:length(t)
            for n=0:1
                if n==0
                    td=d(1:1:48);
                else
                    td=d(49:1:96);
                end
                y(i)=y(i)+1/sqrt(NHT_SIG)*Wix*Wt(TR,T,t(i)-n*TSYM)*...
                    (j*(sum(td.*exp(j*2*pi*HTSIG_M([0:1:47])*df*(t(i)-n*TSYM-TGI-TCS))))+...
                    p(state+n)*sum(Ppilot.*exp(j*2*pi*[-26:1:26]*df*(t(i)-n*TSYM-TGI-TCS))));
            end
        end 
        StateChange(2)%change state for p();
        
        elseif strcmp(field,'DATA') %DATA
            fprintf('Generating DATA....\n');
            iss=iss;
            Nss=Nss;
            NHT_DATA=56;
            TSYM=4*10^-6;
           
            NSR=28;
            NSD=52;
            rand('state',0);
            if strcmp(Mod,'BPSK')
                NDBPS=1;
                BPSK=[-1,1];
                DATA=BPSK(randi(2,[1,NSD]));
                fprintf('Fix Data and Map to BPSK\n')
            elseif strcmp(Mod,'QPSK')
                NDBPS=2;
                QPSK=[-1,1];
                DATA=(QPSK(randi(2,[1,NSD]))+j*QPSK(randi(2,[1,NSD])))/sqrt(2);
                fprintf('Fix Data and Map to QPSK\n')
            elseif strcmp(Mod,'16QAM')
                NDBPS=4;
                QAM16=[-3 -1 1 3];
                DATA=(QAM16(randi(4,[1,NSD]))+j*QAM16(randi(4,[1,NSD])))/sqrt(10);
                fprintf('Fix Data and Map to 16QAM\n')
            elseif strcmp(Mod,'64QAM')
                NDBPS=6;
                QAM64=[-7 -5 -3 -1 1 3 5 7];
                DATA=(QAM64(randi(8,[1,NSD]))+j*QAM64(randi(8,[1,NSD])))/sqrt(42);
                fprintf('Fix Data and Map to 64QAM\n')
            end
            NSYM= ceil((8*HTLENGTH+16+6*NES)/NDBPS);%number of symbols
            
            fprintf('HT-Legnth = ');fprintf(num2str(HTLENGTH));fprintf('\n');
            fprintf('Total number of Symbols of Data = ');fprintf(num2str(NSYM));fprintf('\n');
            
            if(TSYM*NSYM>=max(t))
                fprintf('Warning !! Input time is too short to contain whole packet !! \n')
            end
            
            %------Wait Bar-----
            h = waitbar(0,'Generating DATA...');
            
            for i=1:length(t)
                waitbar(i/length(t))
                
                tp3=0;
                for n=0:NSYM-1
                    if (Wt(TR,TSYM,t(i)-n*TSYM)==0)
                       continue; 
                    end
                    
                    tp1=0;
                    for k=0:NSD-1
                        tp1=tp1+DATA(k+1)*exp(j*2*pi*DATA_M(k)*df*(t(i)-TGI-n*TSYM-TCS));
                    end
                    
                    tp2=0;
                    for k=-NSR:NSR
                        tp2=tp2+Pilot(Nss,iss,n,k)*exp(j*2*pi*df*(t(i)-TGI-n*TSYM-TCS));
                    end
                    tp3=tp3+Wt(TR,TSYM,t(i)-n*TSYM)*(tp1+tp2*p(state+n));
                    
                end
                y(i)=1/sqrt(NHT_DATA)*tp3;
            end 
            close(h);%close Wait bar
            StateChange(NSYM);
            
        else
            error('Field Generator Error !');
    end
        nstate=state;%output new state;
        
        
    
        
        
    function y=Wt(TR,T,t)
        y=zeros(1,length(t));
        for i=1:length(t)
            if -1*TR/2 < t(i) && t(i)<1*TR/2
                y(i)=(sin(pi/2*(0.5+t(i)/TR)))^2;
            elseif 1*TR/2 <= t(i)&& t(i)<=T-1*TR/2
                y(i)=1;          
            elseif T-TR/2<=t && t<=T+TR/2
                y(i)=(sin(pi/2*(0.5-(t(i)-T)/TR)))^2;
            else
                y(i)=0;
            end
        end
        
    function c=calc_crc(bits)
            c=ones(1,8);

        for i=1:length(bits)
            next_c=zeros(1,8);
            next_c(1) = bitxor(bits(i),c(8));
            next_c(2) = bitxor(bitxor(bits(i),c(8)),c(1));
            next_c(3) = bitxor(bitxor(bits(i),c(8)),c(2));
            next_c(4) = c(3);
            next_c(5) = c(4);
            next_c(6) = c(5);
            next_c(7) = c(6);
            next_c(8) = c(7);
            c = next_c;
        end
        for i=1:8;
           if c(i)==1;
               c(i)=0;
           else
               c(i)=1;
           end
        end
        c=c([end:-1:1]);
       
    function y=HTSIG_M(k)
        y=zeros(1,length(k));
        for i=1:length(k)
            if 0<=k(i) && k(i)<=4
                y(i)=k(i)-26;
            elseif 5<=k(i) && k(i)<=17
                y(i)=k(i)-25;
            elseif 18<=k(i) && k(i)<=23
                y(i)=k(i)-24;
            elseif 24<=k(i) && k(i)<=29
                y(i)=k(i)-23;
            elseif 30<=k(i) && k(i)<=42
                y(i)=k(i)-22;
            elseif 43<=k(i) && k(i)<=47
                y(i)=k(i)-21;
            end
        end
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
     function y= Pilot(Nss,iss,n,k)
        pilot=[1 1 1 -1;...
            1 1 -1 -1;...
            1 -1 -1 1;...
            1 1 -1 -1;...
            1 -1 1 -1;...
            -1 1 1 -1;...
            1 1 1 -1;...
            1 1 -1 1;...
            1 -1 1 1 ;...
            -1 1 1 1 ];
        pt=[];
        if Nss==1 && iss==0
            pt=pilot(1,:);
        elseif Nss==2 && iss==0
            pt=pilot(2,:);
        elseif Nss==2 && iss==1
            pt=pilot(3,:);
        elseif Nss==3 && iss==0
            pt=pilot(4,:);
        elseif Nss==3 && iss==1
            pt=pilot(5,:);
        elseif Nss==3 && iss==2
            pt=pilot(6,:);
        elseif Nss==4 && iss==0
            pt=pilot(7,:);
        elseif Nss==4 && iss==1
            pt=pilot(8,:);
        elseif Nss==4 && iss==2
            pt=pilot(9,:);
        elseif Nss==4 && iss==3
            pt=pilot(10,:);
        else
            error('Error Map to Pilot !!')
             
        end
        y=[zeros(1,7),pt(mod(n,4)+1),zeros(1,13),pt(mod(n+1,4)+1),zeros(1,13),pt(mod(n+2,4)+1),zeros(1,13),pt(mod(n+3,4)+1),zeros(1,7)];
        y=y(k+29);% k start from -28
     
     function StateChange(z)
         global state
         state=state+z;
         
     function y=p(idx)
         seq=[1,1,1,1 ,-1,-1,-1,1, -1,-1,-1,-1, 1,1,-1,1, -1,-1,1,1, -1,1,1,-1, 1,1,1,1, 1,1,-1,1,...
             1,1,-1,1, 1,-1,-1,1, 1,1,-1,1, -1,-1,-1,1, -1,1,-1,-1, 1,-1,-1,1, 1,1,1,1, -1,-1,1,1,...
             -1,-1,1,-1, 1,-1,1,1, -1,-1,-1,1, 1,-1,-1,-1, -1,1,-1,-1, 1,-1,1,1, 1,1,-1,1, -1,1,-1,1,...
             -1,-1,-1,-1, -1,1,-1,1, 1,-1,1,-1, 1,1,1,-1, -1,1,-1,-1, -1,1,1,1, -1,-1,-1,-1, -1,-1,-1];
         y=seq(mod(idx,126)+1);% Matlab no negative index
         
     function y=interleaver(NCBPS,NBPSC)
         s = max(NBPSC/2,1);
         I=zeros(1,length(NCBPS));
         for k = 0:NCBPS-1
            I(k+1) =(NCBPS/16)*(mod(k,16)) + floor(k/16);
         end
         
         J=zeros(1,length(NCBPS));
         for t= 0: NCBPS-1
             J(t+1) = s*floor(I(t+1)/s) + mod((I(t+1) + NCBPS-floor(16 * I(t+1)/NCBPS)),s);
         end
         y=J+1;% matlab no zero index
     
      function y=Mapper(bits, NBPSC)
          y=zeros(1,length(bits));
          %NBPSC 2:QPSK 4:16QAM,6:64QAM
          if NBPSC==1
              fprintf('Map to BPSK\n')
              Kmod=1;
          elseif NBPSC==2
              fprintf('Map to QPSK\n')
              Kmod=1/sqrt(2);
          elseif NBPSC ==4
              fprintf('Map to 16QAM\n')
              Kmod=1/sqrt(10);
          elseif NBPSC==6
              fprintf('Map to 64QAM\n')
              Kmod=1/sqrt(42);
          else 
              error('Mapper Error !!')
               
          end
          
          for i=1:length(bits)/NBPSC
             if NBPSC==1
                 if bits(i)==0
                    y(i)=-1;
                 else
                    y(i)=1;
                 end
             else
                idx1=[(i-1)*NBPSC+1:1:(i-1)*NBPSC+NBPSC/2];
                idx2=[(i-1)*NBPSC+NBPSC/2+1:1:i*NBPSC];
                idx=[idx1;idx2];
                for t=1:2
                    ty=zeros(1,length(bits));
               
                    if NBPSC==2 %QPSK
                        if bits(idx(t,:))==[0];
                            ty(idx(t,:))=-1;
                        elseif bits(idx(t,:))==[1];
                            ty(idx(t,:))=11;
                        else
                            error('Error Map to QPSK')
                        end
                    end
                    if NBPSC==4 %16QAM 
                        if bits(idx(t,:))==[0 0];
                            ty(idx(t,:))=-3;
                        elseif bits(idx(t,:))==[0 1];
                            ty(idx(t,:))=-1;
                        elseif bits(idx(t,:))==[1 1];
                            ty(idx(t,:))=1;
                        elseif bits(idx(t,:))==[1 0];
                            ty(idx(t,:))=3;
                        else
                            error('Error Map to 16QAM')
                        end
                    end
                    if NBPSC==6 %64QAM 
                        if bits(idx(t,:))==[0 0 0];
                        ty(idx(t,:))=-7;
                        elseif bits(idx(t,:))==[0 0 1];
                        ty(idx(t,:))=-5;
                        elseif bits(idx(t,:))==[0 1 1];
                        ty(idx(t,:))=-3;
                        elseif bits(idx(t,:))==[0 1 0];
                        ty(idx(t,:))=-1;
                        elseif bits(idx(t,:))==[1 1 0];
                        ty(idx(t,:))=1;
                        elseif bits(idx(t,:))==[1 1 1];
                        ty(idx(t,:))=3;
                        elseif bits(idx(t,:))==[1 0 1];
                        ty(idx(t,:))=5;
                        elseif bits(idx(t,:))==[1 0 0];
                            ty(idx(t,:))=7;
                        else
                            error('Error Map to 64QAM')
                             
                        end
                    end
                
                    if t==1
                        y=ty;
                    elseif t==2
                        y=y+ty*j;
                    end
                end
             end
          end
          y=y*Kmod;
          
         function y=ConvCoder(data)
                  trellis=poly2trellis(7,[133,171]);
                  y=convenc(data,trellis);
                  
         function output_bits = scrambler_and_descrambler(input_bits, initial_state)
                  states = initial_state; 
                  %BE AWARE! example: [0 1 1 0 1 0 0] -> [MSB ... LSB] that means... [state7 ... state1]
                       for n=1:length(input_bits)
                           save_state_4 = states(4);
                           save_state_7 = states(7);
                           output_bits(n) = xor(input_bits(n),xor(states(4),states(7)));
                           states = circshift(states',1)';
                           states(1)= xor(save_state_4,save_state_7);
                        end
          
          function [NCBPS,Mod,NBPSC,Nss,NES]=MCStoParameters(MCS)
              if MCS==0;
                  NCBPS=52;Nss=1;
                  NBPSC=1;NES=1;
                  Mod='BPSK';
              elseif MCS==1;
                  NCBPS=104;Nss=1;
                  NBPSC=2;NES=1;
                  Mod='QPSK';
              elseif MCS==2;
                  NCBPS=104;Nss=1;
                  NBPSC=2;NES=1;
                  Mod='QPSK';
              elseif MCS==3;
                  NCBPS=208;Nss=1;
                  NBPSC=4;NES=1;
                  Mod='16QAM';
              elseif MCS==4;
                  NCBPS=208;Nss=1;
                  NBPSC=4;NES=1;
                  Mod='16QAM';
              elseif MCS==5;
                  NCBPS=312;Nss=1;
                  NBPSC=6;NES=1;
                  Mod='64QAM';
              elseif MCS==6;
                  NCBPS=312;Nss=1;
                  NBPSC=6;NES=1;
                  Mod='64QAM';
              elseif MCS==7;
                  NCBPS=312;Nss=1;
                  NBPSC=6;NES=1;
                  Mod='64QAM';
              elseif MCS==8;
                  NCBPS=104;Nss=2;
                  NBPSC=1;NES=1;
                  Mod='BPSK';
              elseif MCS==9;
                  NCBPS=208;Nss=2;
                  NBPSC=2;NES=1;
                  Mod='QPSK';
              elseif MCS==10;
                  NCBPS=208;Nss=2;
                  NBPSC=2;NES=1;
                  Mod='QPSK';
              elseif MCS==11;
                  NCBPS=416;Nss=2;
                  NBPSC=4;NES=1;
                  Mod='16QAM';
              elseif MCS==12;
                  NCBPS=416;Nss=2;
                  NBPSC=4;NES=1;
                  Mod='16QAM';
              elseif MCS==13;
                  NCBPS=624;Nss=2;
                  NBPSC=6;NES=1;
                  Mod='64QAM';
              elseif MCS==14;
                  NCBPS=624;Nss=2;
                  NBPSC=6;NES=1;
                  Mod='64QAM';
              elseif MCS==15;
                  NCBPS=624;Nss=2;
                  NBPSC=6;NES=1;
                  Mod='64QAM';
              elseif MCS==16;
                  NCBPS=156;Nss=3;
                  NBPSC=1;NES=2;
                  Mod='BPSK';
              elseif MCS==17;
                  NCBPS=312;Nss=3;
                  NBPSC=2;NES=2;
                  Mod='QPSK';
              elseif MCS==18;
                  NCBPS=312;Nss=3;
                  NBPSC=2;NES=2;
                  Mod='QPSK';
              elseif MCS==19;
                  NCBPS=624;Nss=3;
                  NBPSC=4;NES=2;
                  Mod='16QAM';
              elseif MCS==20;
                  NCBPS=624;Nss=3;
                  NBPSC=4;NES=2;
                  Mod='16QAM';
              elseif MCS==21;
                  NCBPS=936;Nss=3;
                  NBPSC=6;NES=2;
                  Mod='64QAM';
              elseif MCS==22;
                  NCBPS=936;Nss=3;
                  NBPSC=6;NES=2;
                  Mod='64QAM';
              elseif MCS==23;
                  NCBPS=936;Nss=3;
                  NBPSC=6;NES=2;
                  Mod='64QAM';
              elseif MCS==24;
                  NCBPS=208;Nss=4;
                  NBPSC=1;NES=2;
                  Mod='BPSK';
              elseif MCS==25;
                  NCBPS=416;Nss=4;
                  NBPSC=2;NES=2;
                  Mod='QPSK';
              elseif MCS==26;
                  NCBPS=416;Nss=4;
                  NBPSC=2;NES=2;
                  Mod='QPSK';
              elseif MCS==27;
                  NCBPS=832;Nss=4;
                  NBPSC=4;NES=2;
                  Mod='16QAM';
              elseif MCS==28;
                  NCBPS=832;Nss=4;
                  NBPSC=4;NES=2;
                  Mod='16QAM';
              elseif MCS==29;
                  NCBPS=1248;Nss=4;
                  NBPSC=6;NES=2;
                  Mod='64QAM';
              elseif MCS==30;
                  NCBPS=1248;Nss=4;
                  NBPSC=6;NES=2;
                  Mod='64QAM';
              elseif MCS==31;
                  NCBPS=1248;Nss=4;
                  NBPSC=6;NES=2;
                  Mod='64QAM';
              else
                  error('Error MCS to Parameters');
                   
              end
              
              
             
  
    