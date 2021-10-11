function [z,rec_m,col_m,m_prime] = receiver_v2(fc, input_signal, M, rt,S,M_order,size_img,m_type)

     N=length(input_signal)/M;
     fl=floor(200);
     %%%BPF at the front-end of the receiver%%%
     %Used to improve SNR by filtering out-band portion of white noise
     %and used to eliminate narrowband interferences
     freqs=[0 0.64 0.68 0.72 0.76 1];    
     amps=[0 0 1 1 0 0]; 
     bpf=firpm(fl,freqs,amps);     
     
     figure(4);
     plotspec(bpf,1/100);
     title('BPF at the receiver');

     x=filter(bpf,1,input_signal); % do the filtering
     
     figure(5);
     plotspec(x,1/2*M);
     title('Filtered received signal');

    
    %%%%%%%%AUTOMATIC GAIN CONTROL%%%%%%%%%%%%%
    len = length(rt);
    average_length = 10;
    a = zeros(1,len);
    a(1) = 1;
    step_size = 0.001;
    sampled = zeros(1,len);
    avg = zeros(1,average_length);
    
    for i=1:len-1
        
        sampled(i) = x(i)*a(i);
        new_sample = sign(a(i))*(sampled(i)^2-S);
        avg = [new_sample avg(1:average_length-1)];
        a(i+1) = a(i)-step_size*mean(avg);
        
    end
     
     figure(6);
     plotspec(sampled,1/M);
     title('Output of AGC');
     figure(7);
     plot(a);
     title('AGC gain (a) convergence');
     
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%Preprocessing of the received signal%%%
    
    rp=sampled.*sampled;     
    freqs=[0 0.54 0.58 0.62 0.66 1]; %freq 60   
    amps=[0 0 1 1 0 0]; 
    bpf=firpm(fl,freqs,amps);          
    figure(8);
    plotspec(bpf,1/M);
    title('BPF to preprocess the received signal for phase recovery');
    rp=filter(bpf,1,rp); 
    figure(9);
    plotspec(rp,1/M);
    title('Preprocessed received signal rp');

%     %%%DIGITAL DOWNCONVERSION TO BASEBAND%%% PLL
%     %%%Note: It is assumed that, transmitter frequency is known(fo=fc)
%     %%%at the receiver and jitter has no effect, so only phase
%     %%%tracking will be implemented by carrier synchronization

    
    freqs=[0 .01 .02 1];
    amps=[1 1 0 0];
    h=firpm(fl,freqs,amps);               
    step_size=0.001;                      
    fo=fc;                                
    theta=zeros(1,len);
    theta(1)=0;                          
    w=zeros(1,fl+1);                     
    for k=1:len-1                        
      w=[w(2:fl+1), rp(k)*sin(4*pi*fo*rt(k)+2*theta(k))];
      new_sample=fliplr(h)*w';               
      theta(k+1)=theta(k)-step_size*new_sample;   
    end
    
    figure(10);
    plot(rt,theta)
    title('Phase Tracking via the Phase Locked Loop')
    xlabel('time'); ylabel('phase offset')


    baseband = sampled.*cos(2*pi*fo*rt+theta(len))*2;
    baseband = [baseband zeros(1,200)];%% son sembolü almak için
    
    figure(11);
    plotspec(baseband,1/M);
    title('correlator input');
        
    freqs=[0 0.05 0.1 1];    
    amps=[1 1 0 0];          
    lpf=firpm(fl,freqs,amps);          
    figure(12);
    plotspec(lpf,1/M);
    title('LPF after carrier synchronization');
    baseband=filter(lpf,1,baseband); % do the filtering
    figure(13);
    plotspec(baseband,1/M);
    title('Digital downconverted baseband signal');

     
    %%%%Correlator%%%
    %%%Enhances peakiness as correlation output%%%
    %%%Same logic with inserting marker sequence to the message bits
    %%%Maximum correlation occurs at the peak values of each symbol
 
    ps=hamming(M);%%Since message shaped with hamming, correlation is achieved with hamming
    % y=filter(fliplr(ps)/(pow(ps)*M),1,baseband);%%Correlation is the same thing with
    %%Filtering baseband signal with flipped hamming. Power is divided
    %%because of power of M samples are included
    y = conv(baseband,fliplr(ps)/(pow(ps)*M));
    figure(14);
    plotspec(y,1/M);
    title('Correlator output');

    %%%Downconversion%%%
    %%%We need peak samples of symbols as soft decisions,
    %%%Timing synchronization is required to find peak sample
    %%%0.5*fl is resulted from LPF delay since its a linear filter(firpm)
    
    tnow=(fl/2)+1; tau=0;
    xs=zeros(1,N);
    tausave=zeros(1,N); 
    tausave(1)=tau; 
    i=0; 
    mu=0.1; 
    delta =0.18;
    
    while tnow<length(y)-fl/2-1 % output power maximization method for time sync. 
        i=i+1;
        xs(i)=interpsinc(y,tnow+tau,fl/8);
        x_deltap=interpsinc(y,tnow+tau+delta, fl/8); 
        x_deltam=interpsinc(y,tnow+tau-delta, fl/8); 
        dx=x_deltap-x_deltam;
        tau = tau+mu*dx*xs(i);
        tnow=tnow+M; %m
        tausave(i)=tau; 
    end
    
%     while tnow<length(y)-fl/2-1    % decision directed method for time sync.     
%       i=i+1;
%       xs(i)=interpsinc(y,tnow+tau,fl/8);   % interp value at tnow+tau
%       x_deltap=interpsinc(y,tnow+tau+delta,fl/8); % value to right
%       x_deltam=interpsinc(y,tnow+tau-delta,fl/8); % value to left
%       dx=x_deltap-x_deltam;             % numerical derivative
%       qx=quantalph(xs(i),[-3,-1,1,3]);  % quantize to alphabet
%       tau=tau+mu*dx*(qx-xs(i));         % alg update: DD
%       tnow=tnow+M; tausave(i)=tau;      % save for plotting
%     end
    
    figure(18);
    subplot(2,1,1), plot(xs(1:i-2),'b.')        % plot constellation diagram
    title('Constellation Diagram');
    ylabel('Estimated Symbol Values')
    subplot(2,1,2), plot(tausave(1:i-2))        % plot trajectory of tau
    ylabel('Offset Estimates'), xlabel('iterations') 

    z=y(fl+abs(tau)*M+M:M:N*M+2*M-1);%set delay to first symbol-sample and increment by M

    %%%Quantization(bits to symbols again considering amplitudes)%%%
    
    alphabet=-(M_order-1):+2:M_order-1;
        
    m_prime=quantalph(z,alphabet)';         % quantize to +/-1 and +/-3 alphabet
    
    figure(19);
    subplot(4,2,1);
    plot([1:3000],baseband(1:3000),'.')
    title('Data at the Input of the Correlator');
    subplot(4,2,2);
    ul=floor((length(baseband)-149)/(4*M));
    plot(reshape(baseband(150:ul*4*M+149),4*M,ul))  % plot the eye diagram
    title('Eye Diagram at the Input of the Correlator')

    subplot(4,2,3);
    plot([1:3000],y(1:3000),'.')
    title('Data at the Output of the Correlator');
    subplot(4,2,4);
    ul=floor((length(y)-149)/(4*M));
    plot(reshape(y(150:ul*4*M+149),4*M,ul))  % plot the eye diagram
    title('Eye Diagram at the Output of the Correlator');
    
    subplot(4,2,5);
    plot([1:length(z)],z,'.');
    title('Constellation Daigram at the Output of the Downsampler')
    subplot(4,2,6);
    ul=floor((length(z))/(4));
    plot(reshape(z(1:ul*4),4,ul));  % plot the eye diagram
    title('Eye Diagram at the Output of the Downsampler');

    if m_type == 0
        
        demod_m=pamdemod(m_prime,M_order);
        bin_m = de2bi(demod_m);   
        L=numel(bin_m);
        col_m = reshape(bin_m,L,1);   % Reshape data into binary k-tuples, k = log2(M)

        matrix_m=reshape(col_m,length(col_m)/8,8); %8 for gray image
        dec_m=bi2de(matrix_m);

        if length(size_img)==3
            rec_m=reshape(dec_m,size_img(1),size_img(2),size_img(3));
        elseif length(size_img)==2
            rec_m=reshape(dec_m,size_img(1),size_img(2));
        end
        rec_m=cast(rec_m,'uint8');
    
    else
        reconstructed_message=pam2letters(m_prime)
        rec_m=1;
        col_m=1;
    end

    
    subplot(4,2,7);
    plot([1:length(m_prime)],m_prime,'.');
    title('Constellation Daigram at the Output of the Quantizer')
    subplot(4,2,8);
    ul=floor((length(m_prime))/(4));
    plot(reshape(m_prime(1:ul*4),4,ul));  % plot the eye diagram
    title('Eye Diagram at the Output of the Quantizer');

      
    %%%Reconstruction of message from symbols using decoding%%%
    %reconstructed_message=pam2letters(mprime)