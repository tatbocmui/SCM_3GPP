%=================================================================
% Start the program
    tic;
    clc;
%     clear all;
%     close all;
%    
%=================================================================
% System's parameters
%
NFFT = 512;         % FFT length
G = 212;            % Guard interval length
M_ary = 4;         % QPSK

NofOFDMSymbol = 1;            %Number of data and pilot OFDM symbol
length_data =  NofOFDMSymbol * NFFT;  % the total data length

%------------------------------------------------
% Sampling duration of LTE 5M bandwidth: Watch Table on the Internet
    t_a = 130.2e-9;
% OFDM symbol duration    
    symbol_duration = NFFT*t_a;%TS = t_a*NFFT: formular in page 29
%

%-------------
% Source bites
%-------------

source_data1 = randi([0 1],length_data,2);
source_data2 = randi([0 1],length_data,2);
%--------------------
% bit to symbol coder
%--------------------

symbols1 = bi2de(source_data1);  
symbols2 = bi2de(source_data2);
%----------------------------
% 4 QAM modulator in base band
%----------------------------

QAM_Symbol1 = qammod(symbols1, M_ary);%x1
QAM_Symbol2 = qammod(symbols2, M_ary);%x2

%--------------------------------------
% STBC Encoder and Preparing data pattern
%--------------------------------------
% Alamouti Code: x = [x1 -x2*;x2 x1*]
Data_Pattern1 = []; % Transmitted Signal before IFFT
% m = 0;
for i=0:NofOFDMSymbol-1;
    
    QAM_tem = [];
    for n=1:NFFT;
            QAM_tem = [QAM_tem,QAM_Symbol1(i*NFFT+n)];
    end;
    Data_Pattern1 = [Data_Pattern1;QAM_tem];
    
    QAM_tem = [];
    for n=1:NFFT;
        QAM_tem = [QAM_tem,-conj(QAM_Symbol2(i*NFFT+n))];
    end;
    Data_Pattern1 = [Data_Pattern1;QAM_tem];    
    clear QAM_tem;
end;

Data_Pattern2 = []; % Transmitted Signal before IFFT
% m = 0;
for i=0:NofOFDMSymbol-1    
    QAM_tem = [];
    for n=1:NFFT;
            QAM_tem = [QAM_tem,QAM_Symbol2(i*NFFT+n)];
    end;
    Data_Pattern2 = [Data_Pattern2;QAM_tem];
    
    QAM_tem = [];
    for n=1:NFFT;
        QAM_tem = [QAM_tem,conj(QAM_Symbol1(i*NFFT+n))];
    end;
    Data_Pattern2 = [Data_Pattern2;QAM_tem];    
    clear QAM_tem;
end;

% main program

SER_R = [];

snr_min =0.0;
snr_max = 10.0;
step = 1.0;
SER=[];
% snr=10;
NT=10;   % change number of loop 

d_s = 10;
d_u = 1/2;

for snr = snr_min:step:snr_max; 
    % snr = 10;
    SER_t=0;    
    snr = snr - 10*log10((NFFT-G)/NFFT); %Miss matching effect  
    SNR = 10^(snr/10);
    load ParametersForSpatialChannelModel.mat;
    for t_i=1:NT %(1000*1920)
     %=====================================================================
     % Transmitter
     %===================================================================== 
        rs11_frame = [];
        rs12_frame = [];
        rs21_frame = [];
        rs22_frame = [];

        h11_frame = [];
        h12_frame = [];
        h21_frame = [];
        h22_frame = [];

        t = rand(1) * 1000;
        t_init = t;%Consider in time t

        for i=0:2*NofOFDMSymbol-1    
            %
            OFDM_signal_tem11 = OFDM_Modulator(Data_Pattern1(i+1,:),NFFT,G); 
            OFDM_signal_tem12 = OFDM_Modulator(Data_Pattern2(i+1,:),NFFT,G);

            h11 = SpatialChannelModel(1,1,d_u,d_s,t,N,M,Pn,sigma_SF,G_BS,G_MS,G_BS_theta_n_m_AoD,G_MS_theta_n_m_AoA,theta_BS,theta_MS,theta_n_m_AoD,theta_n_m_AoA,Phi_n_m,Phi_LOS,v,theta_v,K);

            h12 = SpatialChannelModel(1,2,d_u,d_s,t,N,M,Pn,sigma_SF,G_BS,G_MS,G_BS_theta_n_m_AoD,G_MS_theta_n_m_AoA,theta_BS,theta_MS,theta_n_m_AoD,theta_n_m_AoA,Phi_n_m,Phi_LOS,v,theta_v,K);

            h21 = SpatialChannelModel(2,1,d_u,d_s,t,N,M,Pn,sigma_SF,G_BS,G_MS,G_BS_theta_n_m_AoD,G_MS_theta_n_m_AoA,theta_BS,theta_MS,theta_n_m_AoD,theta_n_m_AoA,Phi_n_m,Phi_LOS,v,theta_v,K);

            h22 = SpatialChannelModel(2,2,d_u,d_s,t,N,M,Pn,sigma_SF,G_BS,G_MS,G_BS_theta_n_m_AoD,G_MS_theta_n_m_AoA,theta_BS,theta_MS,theta_n_m_AoD,theta_n_m_AoA,Phi_n_m,Phi_LOS,v,theta_v,K);


            h11_frame = [h11_frame; h11];             
            h12_frame = [h12_frame; h12];
            h21_frame = [h21_frame; h21];
            h22_frame = [h22_frame; h22];

            %Singal pass channel
            rs11 = conv(OFDM_signal_tem11, h11);            
            rs12 = conv(OFDM_signal_tem12, h12); 
            rs21 = conv(OFDM_signal_tem11, h21); 
            rs22 = conv(OFDM_signal_tem12, h22); 

            % The received signal over multhipath channel is created
            rs11 = awgn(rs11,snr,'measured','dB');
            rs12 = awgn(rs12,snr,'measured','dB');
            rs21 = awgn(rs21,snr,'measured','dB');
            rs22 = awgn(rs22,snr,'measured','dB');

            % frame receive 
            rs11_frame = [rs11_frame; rs11];  
            rs12_frame = [rs12_frame; rs12];
            rs21_frame = [rs21_frame; rs21];  
            rs22_frame = [rs22_frame; rs22];

            t = t + 1000*symbol_duration; % Van Duc modified on 29.08.05 to obtained a reliable result in a short time simulation
            clear OFDM_signal_tem11 OFDM_signal_tem12;     
        end;
          %=====================================================================
          % Receiver
          %=====================================================================
    %         SignalPostFFT11 = [];
    %         SignalPostFFT21 = [];        
    %         SignalPostFFT12 = [];
    %         SignalPostFFT22 = [];

        s_est_matrix=[];
        s1_est = [];
        s2_est = [];

        for n=1:NofOFDMSymbol            
            %-----------------------------------------------------------
            % received data in RX
            %-----------------------------------------------------------
            i = 2*n-1;
            % 1st time slot
            rs11i = rs11_frame(i,:) + rs12_frame(i,:); %antenna1 
            rs21i = rs21_frame(i,:) + rs22_frame(i,:); %antenna2                       
            % 2nd time slot
            rs12i = rs11_frame(i+1,:) + rs12_frame(i+1,:); %antenna1 
            rs22i = rs21_frame(i+1,:) + rs22_frame(i+1,:); %antenna2                         

            %-----------------------------------------------------------
            % OFDM demodulator
            %-----------------------------------------------------------              
            Demodulated_signal11i = OFDM_Demodulator(rs11i,NFFT,NFFT,G);  
    %             SignalPostFFT11 = [SignalPostFFT11; Demodulated_signal11i];

            Demodulated_signal21i = OFDM_Demodulator(rs21i,NFFT,NFFT,G);  
    %             SignalPostFFT21 = [SignalPostFFT21; Demodulated_signal21i];

            Demodulated_signal12i = OFDM_Demodulator(rs12i,NFFT,NFFT,G);  
    %             SignalPostFFT12 = [SignalPostFFT12; Demodulated_signal12i];
    %             
            Demodulated_signal22i = OFDM_Demodulator(rs22i,NFFT,NFFT,G);  
    %             SignalPostFFT22 = [SignalPostFFT22; Demodulated_signal22i];


            %-----------------------------------------------------------
            % Channel Gains - done
            %-----------------------------------------------------------   
            h11_1 = h11_frame(i,:);
            H11_1 = fft([h11_1,zeros(1,NFFT-length(h11_1))]);
            h12_1 = h12_frame(i,:);
            H12_1 = fft([h12_1,zeros(1,NFFT-length(h12_1))]);

            h21_1 = h21_frame(i,:);
            H21_1 = fft([h21_1,zeros(1,NFFT-length(h21_1))]);
            h22_1 = h22_frame(i,:);
            H22_1 = fft([h22_1,zeros(1,NFFT-length(h22_1))]);

            h11_2 = h11_frame(i+1,:);
            H11_2 = fft([h11_2,zeros(1,NFFT-length(h11_2))]);
            h12_2 = h12_frame(i+1,:);
            H12_2 = fft([h12_2,zeros(1,NFFT-length(h12_2))]);

            h21_2 = h21_frame(i+1,:);
            H21_2 = fft([h21_2,zeros(1,NFFT-length(h21_2))]);
            h22_2 = h22_frame(i+1,:);
            H22_2 = fft([h22_2,zeros(1,NFFT-length(h22_2))]);

            %--------------------------------------------------------------
            % Zero Forcing Equalizer
            %--------------------------------------------------------------

            s_est = [];
            W = [];

            r11=Demodulated_signal11i;
            r21=Demodulated_signal21i; 

            r12=Demodulated_signal12i;
            r22=Demodulated_signal22i; 
            % vong lap so OFDM symbol, va vong lap cho tung subcarrier
            for m =1:NFFT                                       
                H_tmp = [H11_1(m), H12_1(m); H21_1(m), H22_1(m);...
                         conj(H12_2(m)), -conj(H11_2(m)); ...
                         conj(H22_2(m)), -conj(H21_2(m))];

%                 W = (transpose(conj(H_tmp))*H_tmp+(1/SNR)*eye(2))\transpose(conj(H_tmp));
                W = inv(transpose(conj(H_tmp))*H_tmp)*transpose(conj(H_tmp));
                r = [r11(m); r21(m); conj(r12(m)); conj(r22(m))];
                s_est_tmp = W*r;
                s_est = [s_est, s_est_tmp];
            end; 
            s_est_matrix = [s_est_matrix,s_est];                      

        end  % end of the OFDM index

        % %-----------------------------------------------
        % % Decoder (Because Zero Forcing Equalizer decode STBC)
        % %-----------------------------------------------
        s1 = qamdemod(s_est_matrix(1,:),M_ary);
        s2 = qamdemod(s_est_matrix(2,:),M_ary);      
                
        [number1, ratio1] = symerr(symbols1,s1');
        [number2, ratio2] = symerr(symbols2,s2');

        SER_R = (ratio1+ratio2)/2;         
        SER_t=SER_t+SER_R;     
    end % end of averaging index t_i
    tE=SER_t/t_i;       
    SER = [SER, tE]; 
end  % end of snr
  
%save and plot 
% psnr = snr_min-3:step:snr_max-3; %for two antennas divide the total energy

psnr = snr_min:step:snr_max;

% old- old delta, new- new delta
% save ZF_STBC_Ser_itr_18.mat psnr  SER;  
% save STBC_SCM_ZF_10_0.5.mat psnr  SER;  

semilogy(psnr, SER,'gs-');

xlabel('SNR in dB');
ylabel('SER');
hold on