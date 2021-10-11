clc;
clear all;
close all;

M_order=4; 
k=log2(M_order); 
m_type=0; % 0 for image, 1 for text

% message as text
if m_type ==1
    str='furkan ayberk yagdi';
    m=letters2pam(str);
% message as image
else
    image=imread('panda.jpeg');
    resized=imresize(image, 0.08);
    var_image=var(double(resized(:)))/length(resized(:));

    size_img=size(resized);
    dec_img=de2bi(resized);
    L=numel(dec_img);
    col_img=reshape(dec_img,L,1);

    dataInMatrix = reshape(col_img,length(col_img)/k,k);   % Reshape data into binary k-tuples, k = log2(M)
    dataSymbolsIn = bi2de(dataInMatrix);                   % Convert to integers
    y=pammod(dataSymbolsIn,M_order);
    m=y';
end

po=-1;
t_offset=0.3;

%%TRANSMITTER

N=length(m);
M=200; mup=zeros(1,N*M); mup(1:M:N*M)=m;  

% figure();
% plot(mup(1:10*200))    % baseband signal spectrum
% title('Zero-padded message signal');

% Hamming pulse filter with T/M-spaced impulse response
p=hamming(M);

delayed_p = [zeros(t_offset*M, size(p,2)); p];
% pulse filter output
x=filter(delayed_p,1,mup); % convolve pulse shape with data

% figure(1)
% plotspec(x,1/M)    % baseband signal spectrum
% title('Baseband signal shaped with hamming pulse');

% am modulation with suppressed carrier
t=1/M:1/M:length(x)/M;        
fc=70;                        
c=cos(2*pi*fc*t+po);          
r=c.*x;                       
rt=(1)/M:1/M:length(r)/M;     

% figure(2);
% plotspec(r,1/M);
% title('Modulated signal');

channel_out = channel(r,t,M);

% figure(3);
% plotspec(channel_out,1/M);
% title('Channel output after noise and flat fading')

S = pow(r);
EIRP = S; %%Effective isotropic radiated power (using isotropic radiator(hypothetic)), it radiates same power with source through all directions
lambda = (3*10^8/fc);
R = fc; %%Information rate (net bitrate) is assumed equal to the bitrate since there is no extra use of bits that reduce channel capacity
GR = (16*pi^2*fc^4)/(3*10^8)^2; %%Receiver antenna gain

PR = EIRP/((4*pi/lambda)^2*R^2)*GR; %Transmitted signal power with receiver antenna gain=1



if m_type ==0
    [z,rec_image,col_m,m_prime]=receiver_v2(fc,channel_out,M,rt,PR,M_order,size_img,m_type);
    figure(20);
    subplot(1,2,1);
    imshow(resized);
    title("Transmitted Image");
    subplot(1,2,2);
    imshow(rec_image);
    title("Received Image")

    %symbol recovery error
    e_symbol=y-z';
    symbol_recovery_error=sum(e_symbol.^2)/length(e_symbol) %average squarred error

    %hard decision error
    e_hard_d=y-m_prime';
    hard_error=sum(e_hard_d.^2)/length(e_hard_d) %average squarred error

    %cluster variance error
    e_cv=m_prime-z;
    cluster_variance=sum(e_cv.^2)/length(e_cv) 

    % bit error rate (BER)
    col_im=cast(col_img,'double');
    ber=100*sum(abs(col_m-col_im))/length(col_m)

    %symbol error rate
    ser= 100*sum(abs(sign(y-m_prime')))/length(y)
    
else
    size_img=1;
    [z,rec_image,col_m,m_prime]=receiver_v2(fc,channel_out,M,rt,PR,M_order,size_img,m_type);
    e_cv=m_prime-z;
    cluster_variance=sum(e_cv.^2)/length(e_cv) 
    pererr=100*sum(abs(sign(m_prime-m(1:length(m_prime)))))/length(m_prime)
    
end
    
