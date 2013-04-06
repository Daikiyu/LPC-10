% Especificaciones del estandar LPC-10:
    %So called LPC10 because 10 LP coefficients are used
    %Bandwidth: 2.4kbps
    %Samples/frame : 180 samples
    %Bits/frame: 54 bits
    %Frame Size: 22.5ms = 44.44 frames/sec
    %Target stream : 8khz sampling rate, 16bit quantization
    


%% AN�LISIS DE LOS PAR�METROS DE LA SE�AL DE VOZ POR TRAMAS

%Lectura de la se�al de entrada
[audio_original,fs,nbits]=wavread('PLC.wav');
audio_original=audio_original(:,1);
N=round(fs*(22.5)*10^(-3)); %n�mero de muestras en una trama
L=length(audio_original); %longitud de la se�al de audio
num_tramas=floor(L/N); %redondeamos el n� de tramas

lpc_total=[]; %Inicializamos un vector que almacena los par�metros lpc calculados para transmitirlos
prediccion_total = []; %Inicializamos este vector que contendr� los valores de la prediccion que nos da el filtro 
energia_total = []; %Inicializamos el vector que contiene los valores calculados de energ�a de la se�al predecida
pitch = []; %Inicializamos un vector de salida en el que guardamos los distintos valores calculados en el algoritmo de c�lculo del pitch
sonoridad_total = []; %inicializamos un vector que almacena los distintos resultados del an�lisis de sonoridad de cada trama y nos dice si es sorda o sonora

j=0;

for i=1:num_tramas
    trama=audio_original(j+1:j+N);
    ventana=hamming(N); % Enventanamos la se�al
    trama_enventanada=ventana.*trama;
    
    %FILTRADO LPC10
    coef_A=lpc(trama_enventanada,10); %Calculamos los coeficientes del filtro
    lpc_total=[lpc_total;coef_A]; %Concatenamos los valores de los coeficientes calculados para cada trama
    prediccion=filter(coef_A,1,trama_enventanada); %A la salida del filtro obtenemos la diferencia entre la se�al original y la predecida
    prediccion_total=[prediccion_total prediccion];%Concatenamos los valores de los errores obtenidos a la salida del filtro para cada trama
    
    % ENERGIA DE LA SE�AL y sonoridad
    energia= sum((prediccion.^2)); %Calculamos la energ�a para cada trama y lo almacenamos en energia_total
    energia_total=[energia_total,energia];
    if energia>(10^-4)
        sonoridad=1;
    else
        sonoridad=0;
    end
    
    %ELIMINACI�N DE CONTINUA Y FILTRADO PASO BAJO DE LA SE�AL ENVENTANADA
    valormedio=mean(trama_enventanada);
    trama_noDC =trama_enventanada-valormedio;
    
    fpb=(2/fs)*900; %Aplicamos un filtro paso bajo
    filtro_pb = fir1(25,fpb,'low');
    filtro_pb2 = filtro_pb;
    trama_filtrada = filter(filtro_pb2,1,trama_noDC);
 
    %B�SQUEDA DEL PITCH 
    y = fft(trama_filtrada, N);
    % Cepstrum is IDFT (or DFT) of log spectrum
    c = ifft(log(abs(y)+eps));


    ms1=fs/1000; % 1ms. maximum speech Fx at 1000Hz
    ms20=fs/50;  % 20ms. minimum speech Fx at 50Hz

    % plot waveform
    t=(0:N-1)/fs;        % times of sampling instants     
    ms2=floor(fs*0.002); % 2ms
    ms20=floor(fs*0.02); % 20ms
    [maxi,idx]=max(abs(c(ms2:ms20)));
    f0 = fs/(ms2+idx-1)
    
    %HACE FALTA INCLUIR LOS FILTROS QUE FACILITAN EL C�LCULO DEL PITCH  
    
    
    pitch = [pitch f0]; %Almacenamos en este vector los distintos valores del pitch para cada trama
    sonoridad_total = [sonoridad_total sonoridad]; %Almacenamos los distintos valores de sonoridad para cada trama
    j=j+N;
end;

%% OBTENCION DE LA SE�AL A PARTIR DE LA PREDICCION

audio_recuperado=[];
[F1,C1]=size(prediccion_total);
for i=1:F1
    trama_prediccion=(:,i);
    coef_lpc=lpc_total(i,:);
    original = filter(1,param_lpc,residuo_d);
    audio_recuperado = [audio_recuperado;original]; 
end;

%% Obtenci�n de la se�al a partir de los par�metros

voz_sintetizada=[];
for i=1:length(pitch)
    if (sonoridad_total(1,i) == 0)% Si el sonido es sordo,excitamos con una fuente de ruido aleatorio
        ruido = rand(1,N);
%         prediccion_sintetizadaN=ruido./((sum(ruido).^2));% Normalizamos la energia de la predicci�n
%     prediccion_sintetizada=prediccion_sintetizadaN.*sqrt(energia_total(1,i));% Asignamos a la trama el valor de energ�a original 
%     voz_sintetizada = [voz_sintetizada prediccion_sintetizada];%generamos nuestra se�al total por concatenaci�n
        prediccion_sintetizada = filter(1,lpc_total(i,:),ruido);% Obtenemos la predicci�n mediante filtrado
    else
        tren = zeros(1,N);% Creamos un tren de deltas de la frecuencia del pitch
        salto=0.000125;%segundos
        k=1/(salto*pitch(i))
        for j=1:round(k):N 
            tren(j)=1; 
        end;
%         prediccion_sintetizadaN=tren./((sum(tren).^2));% Normalizamos la energia de la predicci�n
%     prediccion_sintetizada=prediccion_sintetizadaN.*sqrt(energia_total(1,i));% Asignamos a la trama el valor de energ�a original 
%     voz_sintetizada = [voz_sintetizada prediccion_sintetizada];%generamos nuestra se�al total por concatenaci�n
        prediccion_sintetizada = filter(1,lpc_total(i,:),tren);% Obtenemos la predicci�n mediante filtrado
    end
    prediccion_sintetizadaN=prediccion_sintetizada./((sum(prediccion_sintetizada).^2));% Normalizamos la energia de la predicci�n
    prediccion_sintetizada=prediccion_sintetizadaN.*sqrt(energia_total(1,i));% Asignamos a la trama el valor de energ�a original 
    voz_sintetizada = [voz_sintetizada prediccion_sintetizada];%generamos nuestra se�al total por concatenaci�n
end;