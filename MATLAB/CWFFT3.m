function [peak,freq]=CWFFT3(w,res)
% w waveform, res FFT 2^rest, xlim flag for 20KHz
% close all 
Fs = 3.5e5;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = size(w,1);                     % Length of signal
% t = (0:L-1)*T;                % Time vector
% shift=(min(w)+max(w))/2;
shift=mean(w);
w=w-shift;
% figure, plot(w,':x')
h1=hamming(L,'periodic');
b1=blackman(L,'periodic');
f1=flattopwin(L,'periodic');
wh1=h1(:).*w(:);
wb1=b1(:).*w(:);
wf1=f1(:).*w(:);
ws=sgolayfilt(w,3,7);

% Sum of a 50 Hz sinusoid and a 120 Hz sinusoid
% x = 0.1636*sin(2*pi*2e4*t); 
% y = x;     % Sinusoids plus noise
% figure, plot(Fs*t,y,':bo')
% title('Signal Corrupted with Zero-Mean Random Noise')
% xlabel('time (milliseconds)')
% 
% hold on, plot(w,'r*')

NFFT = 2^res; % Next power of 2 from length of y
% Y = fft(y,NFFT)/L;
Y1=fft(w,NFFT)/L;
Y2=fft(wh1,NFFT)/L;
Y3=fft(wb1,NFFT)/L;
Y4=fft(wf1,NFFT)/L;
Y5=fft(ws,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2);

% Plot single-sided amplitude spectrum.
% figure, plot(f,2*abs(Y3(1:NFFT/2)),':bo') 
% hold on, plot(f,2*abs(Y1(1:NFFT/2)),':r*') 
% hold on, plot(f,2*abs(Y2(1:NFFT/2)),':kx') 

% figure, plot(f,2*abs(Y3(1:NFFT/2)),':bo')
% hold on, plot(f,2*abs(Y1(1:NFFT/2)),':r*') 
% hold on, plot(f,2*abs(Y2(1:NFFT/2)),':kx') 
% hold on, plot(f,2*abs(Y4(1:NFFT/2)),':gs') 
% hold on, plot(f,2*abs(Y5(1:NFFT/2)),':y+') 
% figure, plot(f,abs(Y1(1:NFFT/2)),':bo')
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')
% axis([15e3 25e3 0 1.5e5])
% 
% % max(2*abs(Y(1:NFFT/2)))
% 
% 
% 
% 
% 
% 
powY3=sqrt(real(Y3(1:NFFT/2)).^2+imag(Y3(1:NFFT/2)).^2)/NFFT;
powY1=sqrt(real(Y1(1:NFFT/2)).^2+imag(Y1(1:NFFT/2)).^2)/NFFT;
powY2=sqrt(real(Y2(1:NFFT/2)).^2+imag(Y2(1:NFFT/2)).^2)/NFFT;
powY4=sqrt(real(Y4(1:NFFT/2)).^2+imag(Y4(1:NFFT/2)).^2)/NFFT;
powY5=sqrt(real(Y5(1:NFFT/2)).^2+imag(Y5(1:NFFT/2)).^2)/NFFT;
% figure, plot(f,powY3,':bo') 
% hold on,plot(f,powY1,':r*')
% hold on,plot(f,powY2,':kx')

% power3=(abs(Y3).^2)/NFFT;

% 
% ylabel('Power')
%  axis([0 5e4 1e4 7e4])

% max(powY1(1:NFFT/2))
% 

%  if(xlim >0)
% 
% 
% 
% figure, plot(f,20*log10((powY3(1:NFFT/2))),':bo')
% figure, plot(f,20*log10((powY1(1:NFFT/2))),':r*') 
% figure, plot(f,20*log10((powY2(1:NFFT/2))),':kx') 
% figure, plot(f,20*log10((powY4(1:NFFT/2))),':gs') 
% figure, plot(f,20*log10((powY5(1:NFFT/2))),':y+') 
%  xlim([0 5e4])
%  end
 


% display('FFT peak=')
peak=20*log10(max(powY3));
% display('Center Freq(Hz)')
freq=f(find(powY3==max(powY3)));
% figure, plot(f,20*log10((power3(1:NFFT/2))),':rx')
% 
% hold on, plot(f,20*log10((powY3(1:NFFT/2))),':bo')
% xlim([0 5e4])
% 
% axis([0 5e4 -160 -100])
% 
% max(20*log10((powY1(1:NFFT/2))))

end