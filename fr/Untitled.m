%% Homework 3
%Q1
load('Hwk2.mat')
sampling_frequency=250; % Set Sampling frequency(line 1) 
dt=1/sampling_frequency; % Sampling interval (line 2) 
rec_time=length(data(1,:))./250; % total recording time (line 3) 
freq_dat_1= fft(data(1,:)); %Calculate the Fourier Transform (line 4) 
Pow_spec=(2*(dt^2)/rec_time)*abs(freq_dat_1); %Do the calculation for the power spectrum (line 5) 
Pow_spec2=Pow_spec(1:(length(data(1,:))/(2))+1); % Only use the beginning half of the spectrum vector (line 6) 
  
df=1/max(rec_time); % Frequency resolution (line 7) 
fnq= sampling_frequency/2; % Nyquist frequency= half the sampling frequency. (line 8) 
freq_axis=(0:df:fnq); %Gives you the frequency axis. (line 9) 
  
plot(freq_axis,Pow_spec2) % Plot power spectrum (line 10) 
xlim([0 80]) % Set x-axis limits (line 11) 
xlabel('Frequency Hz') % (line 12) 
ylabel ('Power') %(line 13)

%Q3
old_freq_dat_1=abs(freq_dat_1);
old_freq_dat_1(1:10)
freq_dat_1(1:10)

%Q4

freq_dat_22= fft(data(2,:)); 
Pow_spec3=(2*(dt^2)/rec_time)*abs(freq_dat_22);
Pow_spec4=Pow_spec3(1:(length(data(2,:))/(2))+1); 
figure,
hold on
plot(freq_axis,Pow_spec2) 
plot(freq_axis,Pow_spec4)
xlim([0 80]) 
xlabel('Frequency Hz')
ylabel ('Power') 

%Q5
dat_22_mean_sub=(data(2,:)-mean(data(2,:)));
freq_dat_22_mean= fft(dat_22_mean_sub);
Pow_spec5=(2*(dt^2)/rec_time)*abs(freq_dat_22_mean);
Pow_spec6=Pow_spec5(1:(length(data(2,:))/(2))+1);
figure,
hold on
plot(freq_axis,Pow_spec4) 
plot(freq_axis,Pow_spec6)
xlabel('Frequency Hz')
ylabel ('Power')

%Q6
[V,m,h,n,t]=hhrun_sv(20, 500, -65, 0, 0, 0, 0); 
hold on
plot(t,V)
xlabel('Time (ms)')
ylabel ('Voltage (mV)')
[V,m,h,n,t]=hhrun_sv(50, 500, -65, 0, 0, 0, 0); 
plot(t,V)
xlabel('Time (ms)')
ylabel ('Voltage (mV)')

%Q7
figure,
[V,m,h,n,t]=hhrun_sv(20, 500, -65, 0, 0, 0, 0); 
hold on
plot(t,V)
[V,m,h,n,t]=hhrun_sv(20, 500, -65, 0, 0, 0, 0); 
plot(t,V)
xlabel('Time (ms)')
ylabel ('Voltage (mV)')