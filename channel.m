function channel_output = channel(channel_input,rt,M)

len = length(rt);

sigma_broad = 0.5;
mu = 0;

broad_noise =sigma_broad*randn(1,len)+mu;

A1 = 0;
A2 = 0;

narrow_fr1 = randi(M/2);
narrow_fr2 = randi(M/2);

narrow_noise1 = A1*cos(2*pi*narrow_fr1*rt);
narrow_noise2 = A2*cos(2*pi*narrow_fr2*rt);

range = 0.4;
h=zeros(1,len);
h = [ones(1,floor(range*len)) 0.8*ones(1,len-floor(range*len))];
%%2 Narrowband noise that caused from nearby working systems
%%Broadband(white) noise in the background

channel_output = h.*channel_input+broad_noise+narrow_noise1+narrow_noise2;
