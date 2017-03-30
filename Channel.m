function [output] = Channel(input,mode,dt)
t = 0:dt:3e-7;
if mode == 1    % 2-ray model
    h1 = zeros(1,length(t));
    h1(1) = 1;
    h1(6) = 1;
    %h1 = sqrt(1/length(t)/((norm(h1))^2/length(h1)))*h1;
    output = conv(input,h1);
    %figure(2)
    %stem(h1)
    %set(gca,'XTick',1:length(t)); % Change x-axis ticks
    %set(gca,'XTickLabel',t);      % Change x-axis ticks labels to desired values.
    %title('2-Ray Model Impulse Response');
    %grid on
else            % 4-ray model
    h2 = zeros(1,length(input));
    h2(1)  = 1;
    h2(8)  = 10^(-3/10);
    h2(16) = 10^(-6/10);
    h2(21) = 10^(-10/10);
    %h2 = sqrt(1/length(t)/((norm(h2))^2/length(h2)))*h2;
    output = conv(input,h2);
    %figure(2)
    %stem(h2)
    %set(gca,'XTick',1:length(input)); % Change x-axis ticks
    %set(gca,'XTickLabel',t);      % Change x-axis ticks labels to desired values.
    %title('4-Ray Model Impulse Response');
    %grid on
end    