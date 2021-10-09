% Please keep all the 3files in the same folder, the files are
%LMSRAF, fliterlengthopt and the main code file
%ANCfinalwithpath_propagation
%After running the code, 3 sounds can be heard.
%The first one is the original signal, 2nd one is the distorted signal and
%the third one is the corrected signal.


clear all;
clc;

t = 0:.0001:1;
xc = 25.*sin((2*pi*450).*t) + 15.*cos((2*pi*500).*t) + 20.*sin((2*pi*450).*t + pi/3);%undistorted signal.

disp('original sound');
soundsc(xc,10000);
disp('end');
pause(3);



df = randi([200,2000],1,10);
dfs = zeros(1,length(xc));
for i = 1:length(df)
    dfs = dfs + randi([3 10],1,1).*sin((2*pi*df(i)).*t); %external interference or unwanted signal

end
xrefc = zeros(1,length(xc));
for i = 7:40:3207
    xrefc = xrefc + 5*cos((2*pi*i).*t);
end
h_pathc = (sinc((2*pi*500).*t)).^2;
path_distc = conv(h_pathc,xc+xrefc);%signal being distorted by path.
path_dist = path_distc(1:length(t));
dc = path_dist + dfs ;%distorted signal

disp('distorted sound');
soundsc(dc,10000);
pause(3);

% Now,the signals would be sampled.

fs = 10000;
ts = 1/fs;
Ns = fs*max(t);
n = 0:Ns-1;
f = -fs/2:fs/Ns:fs/2-fs/Ns;%frequency axis variable in Hz
x = xc(1:((length(t)-1)/Ns):length(t)-1);%sampling of signals
d = dc(1:((length(t)-1)/Ns):length(t)-1);
xref = xrefc(1:((length(t)-1)/Ns):length(t)-1);
h_path = h_pathc(1:((length(t)-1)/Ns):length(t)-1);
%signals being analyzed in frequency domain.
X = fftshift(fft(x,length(f)));
Xref = fftshift(fft(xref,length(f)));
D = fftshift(fft(d,length(f)));
H_Path = fftshift(fft(h_path,length(f)));
M = length(x);

Dist_X = zeros(1,length(X));
Dist_X(find(abs(X)>(.01*max(abs(X))))) = D(find(abs(X)>(.01*max(abs(X)))));
dist_x = real(ifft(ifftshift(Dist_X)));
Dist_Xref = zeros(1,length(X));
Dist_Xref(find(abs(Xref)>(.5*max(abs(Xref))))) = D(find(abs(Xref)>(.5*max(abs(Xref)))));
dist_xref = real(ifft(ifftshift(Dist_Xref)));
Interf = D-Dist_X-Dist_Xref;
interf = real(ifft(ifftshift(Interf)));

P = sum((x+xref).^2)/M;
DBI = 10*log10(sum(x.^2)/sum((d-x).^2)); %SNR of distorted signal
% // No need to run this code but can be used to find the perfect 'N'
% E = filterlengthopt(x,dist_x)                       
% N = min(round(M/6)+find(E < ((max(d)/100)+min(E))))


N = 100;


del = 1/(10*N*P);
[y,h] = LMSRAF(dist_x + dist_xref,x + xref,N,del);       %LMS algorithm for adaptive filter
z = y(1:M);
H = fftshift(fft(h,length(f)));%adaptive filter frequency domain
SPK = -(1./H).*(Interf + Dist_X + Dist_Xref -X);
spk = ifft(ifftshift(SPK)); %anti_noise generated at the ANC speaker.
ANC = H_Path.*SPK;
ans = real(ifft(ifftshift(ANC)));; %anti_noise signal at interference point.
ansc = interp1((ts.*n),ans,t);


output = d + ans;                    %anti_noise added to distorted signal
OUTPUT = fftshift(fft(output));
output_real = interp1((ts.*n),output,t);

disp('corrected sound');
soundsc(output_real,10000);


DBO = 10*log10(sum(x.^2)/sum((x-output).^2)); % SNR of OUTPUT signal

figure(1)
subplot(4,1,1),plot(n(end-200:end),x(end-200:end));
xlabel('Samples'),ylabel('original signal[n]')
subplot(4,1,2),plot(n(end-200:end),d(end-200:end));
xlabel('Samples'),ylabel('distorted signal[n]')
subplot(4,1,3),plot(n(end-200:end),ans(end-200:end));
xlabel('Samples'),ylabel('anti_noise[n]')
subplot(4,1,4),plot(n(end-200:end),output(end-200:end));
xlabel('Samples'),ylabel('corrected_signal[n]')
figure(2)
subplot(2,1,1),plot(f,abs(H_Path));
xlabel('Frequency(Hz)'),ylabel('Frequency response of Propagation Path')
subplot(2,1,2),plot(f,abs(H));
xlabel('Frequency(Hz)'),ylabel('Frequency response of adaptive filter')
