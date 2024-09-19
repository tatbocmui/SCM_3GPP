
%---------------------------------------------------------------------------
%  OFDM modulator
%  NFFT: FFT length
%  chnr: number of subcarrier
%  G: guard length
%  N_P: channel impulse response length
%---------------------------------------------------------------------------

function [y] = OFDM_Demodulator(data,chnr,NFFT,G)


x_remove_guard_interval = data(G+1:NFFT+G); % remove the guard interval

x = fft(x_remove_guard_interval); % should be: x = (1/sqrt(N))*fft(x_remove_guard_interval);

y = x(1:chnr); %Zero removing



	
