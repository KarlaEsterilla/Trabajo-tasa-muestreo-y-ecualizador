%clc
%clear all
%close all

%% Leer audio
[a,fs]=audioread("rush.wav");
info = audioinfo("rush.wav");

%% Características del audio
d=info.Duration;
Fs =info.SampleRate;  % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(a);        % Length of signal
t = (0:L-1)*T;        % Time vector

%% canales de audio
a_m = a(:,1); % canal izquierdo
a_m=a_m.';

%% Dibujar sonido original en el tiempo
figure
plot(t,a_m);
title('sonido original en el tiempo');
xlabel('time');

%% Decimacion e Interpolacion
%Decimacion=1
%Interpolacion=2
figure
select = 1;
r = 10;
if select == 2
    xn = decimate(a_m,r);
    t = linspace(0,d,ceil(length(t)/r));
    stem(t,xn);
    title('decimacion');
else
    xn = interp(a_m,r);  
    t=linspace(0,d,(length(t)*r));
    stem(t,xn);
    title('interpolacion');
end

%% Transformada de fourier de la señal muestriada

figure

A_m =fftshift( fft(xn));
f = linspace(-Fs/2,Fs/2,length(A_m));
plot(f,abs(A_m)/max(abs(A_m)));

title('sonido muestriado en frecuencia');
xlabel('hz');

%% Buscar Frecuencia maxima 

[pks,locs]=findpeaks(abs(A_m),f);
fc = (max(locs));
hold on
filtro=1.*((abs(f)<=fc)&((abs(f)>=0)));
plot(f,filtro,'r');

%% Filtrado antialiasing
A_lpf=A_m.*filtro;

%% Bandas para la ecualizacion

%Ganancias
G0=db2mag(1.2);
G1=db2mag(1);
G2=db2mag(0.5);
G3=db2mag(0);
G4=db2mag(0);
G5=db2mag(0);

figure

filtro_Sub_Bass = G0.*((abs(f)<=60)&((abs(f)>=16)));
plot(f,filtro_Sub_Bass,'r','LineWidth',1.5);
hold on
filtro_Bass = G1.*((abs(f)<=250)&((abs(f)>63)));
plot(f,filtro_Bass,'g','LineWidth',1.5);
hold on
filtro_Low_Mids = G2.*((abs(f)<=2000)&((abs(f)>253)));
plot(f,filtro_Low_Mids,'b','LineWidth',1.5);
hold on
filtro_High_Mids = G3.*((abs(f)<=4000)&((abs(f)>2003)));
plot(f,filtro_High_Mids,'c','LineWidth',1.5);
hold on
filtro_Presence = G4.*((abs(f)<=6000)&((abs(f)>4003)));
plot(f,filtro_Presence,'y','LineWidth',1.5);
hold on
filtro_Brilliance = G5.*((abs(f)<=16000)&((abs(f)>6003)));
plot(f,filtro_Brilliance,'m','LineWidth',1.5);
hold on
plot(f,abs(A_lpf)/max(abs(A_lpf)),'k');


%% Ecualizacion

figure
A_lpf=A_lpf.*(filtro_Sub_Bass+filtro_Bass+filtro_Low_Mids+filtro_High_Mids+filtro_Presence+filtro_Brilliance);
plot(f,abs(A_lpf)/max(abs(A_lpf)));
title('Audio ecualizado');

%% Transformada inversa

figure
a_lpf = real(ifft(fftshift(A_lpf)));
plot(t,a_lpf)
title('Audio ecualizado dominio del tiempo');


















    






