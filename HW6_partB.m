clear
itn = 4e3; % # of data samples.... 4000 symbols

SNR_dB=20; 
%SNR_dB=15;
%SNR_dB=10; 
sigman2 = 10^(-SNR_dB/10); 
sigman=sqrt(sigman2);

ch = [1 0.8]; % channel coefficients
N = 8; %number of taps
mu=.01; % step size

runs = 1; % number of independent trials

xi=zeros(itn,1);
for k=1:runs
    inphase=(randi(2,itn,1)-1)*2-1; %%inphase modulation
    quad1 = (randi(2,itn,1)-1)*2-1; %%quad component
    quad= 1j*quad1; %%quad modulation
    inphaseBits = (inphase>0); %%inphase bits
    quadBits = (quad1>0);      %%quad bits
    x = (inphase + quad)/sqrt(2);  %%qpsk signal
    noise = sigman*(randn(itn,1)+1j*randn(itn,1))/sqrt(2);
    d = filter(ch,1,x) + noise; % output signal without equalizer
    d_inphase = real(d(2001:4000));   % real component W/O equalizer
    d_quad = imag(d(2001:4000));   % imaginary component W/O equalizer
    d_inphase_bits =(d_inphase>0); % real bits without equalizer
    d_quad_bits =(d_quad>0); % quad bits without equalizer
    d_inphase_error = sum(abs(d_inphase_bits-inphaseBits(2001:4000))); %inphase bit error 
    d_quad_error = sum(abs(d_quad_bits-quadBits(2001:4000))); %inphase bit error
    d_total_error = (d_quad_error + d_inphase_error)/itn; %  
    weights = zeros(N,1);

    for n=N:itn; %%for n = 8:4000????
    xtdl = d(n:-1:n-N+1); %%y(t)???
    xout(n) = weights' * xtdl; % output signal w/ equalizer  (w' * y_k)
    e = x(n) - xout(n);    %% e = x_k - w_k'y_k = error
    xi(n) = xi(n) + abs(e)^2;  %% MSE ???
   weights = weights + mu * conj(e) * xtdl; %% w(t+1) = w(t) + mu * e* * y(t)
   
    end
end

xoutT = transpose(xout(2001:4000));  %%output with eq transpose last 2000 symbols
eq_inphase = real(xoutT);   %inphase component
eq_quad = imag(xoutT);      %quad component
eq_inphase_bits = (eq_inphase>0);
eq_quad_bits = (eq_quad>0);
eq_inphase_error = sum(abs(eq_inphase_bits - inphaseBits(2001:4000))); %eq inphase bit error 
eq_quad_error = sum(abs(eq_quad_bits-quadBits(2001:4000))); %eq inphase bit error
eq_total_error =(eq_inphase_error+eq_quad_error)/(2*itn);

xi=xi/runs;
realweights = real(weights);
a = d(3901:4000); %%last 100 symbols of NO EQ
a_real= real(a);  %%last 100 symbols of inphase NO EQ
a1= xout(3901:4000); %last 100 symbols WITH EQ
a1_real = real(a1);  %last 100 symbols of inphase WITH EQ
samples = [1:4000];
last100samples = samples(3901:4000);
plot (last100samples, a_real,'b',last100samples, a1_real, 'r')
title('Equalizer vs No Equalizer')
xlabel('blue=no EQ   red=EQ'); ylabel('In Phase component')

%x2 = [1:8];
%figure
%plot(x2, realweights)
%title('Weights')

%%semilogy(xi,'k');
%%xlabel('No. of iterations'); ylabel('MSE')
