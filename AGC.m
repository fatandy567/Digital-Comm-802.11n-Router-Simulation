function [output,G] = AGC(status,input,K,A)  

if status==0;
    output=input;
    G=0;
    return;
end
output = zeros(1,length(input)); % output signal
z = 0;                 % detector output (running sum of the power of the input signal)
g = 0;                 % instantaneous loop gain
e = 0;                 % error
G = zeros(1,length(input));      % loop gain
for i = 1:length(input)
    output(i) = input(i)*exp(g);
    z = (z + input(i)*conj(input(i)))/i*exp(2*g);
    e = A - log(z);
    g = g + K*e;
    G(i) = g;
end
output(301:end) = input(301:end)*G(300);

%figure(4)
%plot(G);grid on;hold on;
%title('The Gain of the AGC Loop');
%xlabel('Time (10^-^8 s)');
%G(301:end) = G(300);
%plot(G,'r');
