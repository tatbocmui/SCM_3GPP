%
%SUMMARY
% This script file will generate al the parameters used in function
% SpatialChannelModel.m
%
    clc;
    clear all;
%
%-------------------------------------------------------------------
% scenario 
%  scenario = 1 => Suburban macrocell   (approximately 3Km distance BS to BS) 
%  scenario = 2 => Urban macrocell      (approximately 3Km distance BS to BS)
%  scenario = 3 => Urban microcell      (less than 1Km distance BS to BS) 
%
    scenario = 3;
%
%-------------------------------------------------------------------
% various distance and orientation parameters
%
%     R = 10*10^3;
    % cell size in macro-cellular networks in meters
    
%     d = round(R*rand(1));
    % distance between the BS and MS in meters
    
    Omega_BS = 0;
    % angle between the absolute North and the BS antenna broadside
    
    Omega_MS = 0;% round(360*rand(1));
    % angle between the absolute North and the MS antenna broadside
    
    theta_BS = 60;%round(360*rand(1));
    % angle between the BS antenna broadside and the LOS connection
   
    theta_MS = 60;%round(abs(Omega_BS-Omega_MS+theta_BS+180));
    % angle between the MS antenna broadside and the LOS connection
    
    v = 30*1000/3600;%km/h->m/s
    % magnitude of the MS velocity in meters
    
    theta_v = 0;%round(360*rand(1));
    % angle between the direction of the MS and the LOS connection
%
%-------------------------------------------------------------------
%
    N = 20;         % number of paths
    M = 20;         % number of subpaths per path
%
%-------------------------------------------------------------------
% lognormal shadowing standard deviation
%
	if(scenario ~= 3)
		sigma_SF = 10^0.8;	% 8 dB in suburban and urban macro 
							% environments
	else
		sigma_SF = 10^1;		% 4dB (LOS) or 10dB (NLOS) in suburban 
							% and urban macro environment
	end
%
%-------------------------------------------------------------------
% random avarage powers for each of the N multipath components
%
    if(scenario == 1)
		r_DS = 1.4;
		sigma_DS = 0.17*10^-6;
	else
		r_DS = 1.7;
		sigma_DS = 0.65*10^-6;
	end
%     z_n = rand(1,N);
%     
%     for k=0:N-1
%         t_n(k+1)= -r_DS*sigma_DS*log(z_n(k+1));
%     end
%     t_n = sort(t_n);
%     t_n = t_n - t_n(1);
% 
%     s = 10^(0.1*3);
%     x_n = s*randn(1,N);
% 
%     for k=0:N-1
%         P(k+1)=exp((((1-r_DS)*t_n(k+1))/(r_DS*sigma_DS)))*10^(-x_n(k+1)/10);
%     end

    load extended_Pn.mat
%     P_dB = [0, -1, -9, -10, -15, -20];
%     P = 10.^(P_dB/10);
    P = ext_Pn;
    sumP = sum(P);
    Pn = P/sumP;
%
%-------------------------------------------------------------------
% DoDs and AoDs of each of the N multipath components
%
	% for suburban and urban macro environments
	if(scenario ~= 3)	
		if(scenario == 1)
			r_AS = 1.2;
			sigma_AS = 5;
		else
			r_AS = 1.3;
			sigma_AS = 8;  % 8 deg or 15 deg
		end
		sigma_AoD = r_AS*sigma_AS;
		delta_n_AoD = sigma_AoD.*randn(1,N); 
		for n = 1:N-1
			for nn = n+1:N
				if(abs(delta_n_AoD(n))>abs(delta_n_AoD(nn)))
					tmp = delta_n_AoD(n);
					delta_n_AoD(n) = delta_n_AoD(nn);
					delta_n_AoD(nn) = tmp;
				end
			end
		end
		
		sigma_n_AoA = 104.12*(1-exp(-0.2175*abs(10*log10(Pn))));
		delta_n_AoA = sigma_n_AoA.*randn(1,N); 
        
          %------------------------------------
%             % Uniform distribution as used in [2]
%             Phi_MS_0 = 180/(2*N);
%             k = 4;
%             Phi_MS_n = [];
%             for n = 1:N
%                 phi_n = k*180/(2*N)*(n-1/2)+Phi_MS_0;
%                 Phi_MS_n = [Phi_MS_n,phi_n]; 
%             end
%             delta_n_AoA = Phi_MS_n;

        %     %------------------------------------
%         	% Natural uniform distribution
%             delta_n_AoA = 360*rand(1,N);
        %

	else 	% for urban micro environment
		delta_n_AoD = 80*rand(1,N)-40;
		 
		sigma_n_AoA = 104.12*(1-exp(-0.265*abs(10*log10(Pn))));
		delta_n_AoA = sigma_n_AoA.*randn(1,N);
        
        % Natural uniform distribution
        delta_n_AoA = 360*rand(1,N);
        
        % Uniform distribution as used in [2]
            Phi_MS_0 = 180/(2*N);
            k = 4;
            Phi_MS_n = [];
            for n = 1:N
                phi_n = k*180/(2*N)*(n-1/2)+Phi_MS_0;
                Phi_MS_n = [Phi_MS_n,phi_n]; 
            end
%             delta_n_AoA = Phi_MS_n;
    end
		
%--------------------------------------------------------------------
% Phase of the mth subpath of the nth path
%
    Phi_n_m = 2*pi*rand(N,M);
    Phi_LOS = 2*pi*rand;
%
%--------------------------------------------------------------------
% antenna gains of the BS and MS subpaths
%
    load Delta_n_m_AoD_5deg
    load Delta_n_m_AoA_35deg
    theta_n_m_AoD = [];
    theta_n_m_AoA = [];
    theta_n_m_AoD_norm = [];
    theta_n_m_AoA_norm = [];
    for n = 1:N
        for m = 1:M
            theta_n_m_AoD(n,m) = theta_BS+delta_n_AoD(n)+Delta_n_m_AoD_5deg(m);
            theta_n_m_AoA(n,m) = theta_MS+delta_n_AoA(n)+Delta_n_m_AoA_35deg(m);
            
            
            A_m=20;         %3 Sector antenna maximum attenuation in dB
            theta_3dB=70;   %3 Sector 3 dB beamwidth in degrees
            
            if (theta_n_m_AoD(n,m)>180)
                theta_n_m_AoD_norm(n,m) = 360-theta_n_m_AoD(n,m);
            elseif (theta_n_m_AoD(n,m)<-180)
                theta_n_m_AoD_norm(n,m) = 360+theta_n_m_AoD(n,m);
            else
                theta_n_m_AoD_norm(n,m) = theta_n_m_AoD(n,m);
            end
            A_theta_n_m_AoD(n,m) = -min((12*(theta_n_m_AoD_norm(n,m)/theta_3dB)^2),A_m);
            G_BS_theta_n_m_AoD(n,m) = 10^(A_theta_n_m_AoD(n,m)/10);
            
            
%             if theta_n_m_AoA(n,m)>180
%                 theta_n_m_AoA_norm(n,m) = 360-theta_n_m_AoA(n,m);
%             elseif theta_n_m_AoA(n,m)<-180
%                 theta_n_m_AoA_norm(n,m) = 360+theta_n_m_AoA(n,m);
%             else
%                 theta_n_m_AoA_norm(n,m) = theta_n_m_AoA(n,m);
%             end
            A_theta_n_m_AoA(n,m) = 0;
            G_MS_theta_n_m_AoA(n,m) = 10^(A_theta_n_m_AoA(n,m)/10);           
        end
    end
    A_theta_BS = -min((12*(theta_BS/theta_3dB)^2),A_m);
    G_BS = 10^(A_theta_BS/10);
    A_theta_MS = 0;
    G_MS = 10^(A_theta_MS/10);
    
    d = 250;
    K = 10^((13-0.03*d)/10);
%
% save ParametersForSpatialChannelModel.mat theta_BS theta_MS N M Pn sigma_SF G_BS G_MS G_BS_theta_n_m_AoD G_MS_theta_n_m_AoA theta_n_m_AoD theta_n_m_AoA Phi_LOS Phi_n_m v theta_v delta_n_AoA delta_n_AoD K


save ParametersForSpatialChannelModel_uniform_AoA.mat theta_BS theta_MS N M Pn sigma_SF G_BS G_MS G_BS_theta_n_m_AoD G_MS_theta_n_m_AoA theta_n_m_AoD theta_n_m_AoA Phi_LOS Phi_n_m v theta_v delta_n_AoA delta_n_AoD K


















