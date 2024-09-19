
%---------------------------------------------------------------------------
%  OFDM modulator
%  NFFT: FFT length
%  chnr: number of subcarrier
%  G: guard length
%---------------------------------------------------------------------------

function [y] = OFDM_Modulator(data,NFFT,G)

chnr = length(data);

x = [data,zeros(1,NFFT - chnr)]; %Zero padding

a = ifft(x); % should be: fft a = sqrt(N)*ifft(x); 
y = [a(NFFT-G+1:NFFT),a]; % insert the guard interval

	
