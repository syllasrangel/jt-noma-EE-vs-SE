function [h,n] = wireless_channel(h_length, TX_RX_distance)
%Input: Impulse response length, and distance between TX and RX
%
%Output: Impulse response and AWGN

%power_gain = 1/(TX_RX_distance^path_loss_exponent);

PL_dB = 128.1 + 37.6*log10(TX_RX_distance) + sqrt(10)*randn(1); % Macro cell path loss from https://www.arib.or.jp/english/html/overview/doc/STD-T104v4_20/5_Appendix/Rel13/36/36931-d00.pdf
%PL_dB = 38 + 30*log10(TX_RX_distance) + sqrt(6)*randn(1); % Pico cell
power_gain = 10^(-PL_dB/10);

h = sqrt(power_gain)*(((randn(1,h_length)+1i.*randn(1,h_length)))/sqrt(2)); % Impulse response

n = 1/sqrt(2).*(randn(1)+1i*randn(1)); % AWGN
end