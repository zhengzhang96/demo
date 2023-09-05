function [fftmax,energy,FFT_signal,fftFreqs] = fftmax_with_boundarise(input,fs,boundarise)

[a b] = size(input);
if a>b
    input = input';
end
% Add padding to increase the resolution
pad_factor = 10;
windowLength = length(input)*pad_factor;
FFT_signal = fft(input,windowLength);
% Zero padding, increase the frequency resolution
fftFreqs = (0:(windowLength-1))*fs/windowLength;
fftFreqs(fftFreqs >= fs/2) = fftFreqs(fftFreqs >= fs/2)-fs;
total_energy = sum(abs(FFT_signal).*abs(FFT_signal));
FFT_signal(fftFreqs<boundarise(1) | fftFreqs>boundarise(2)) = 0;

[max_value,freqMax] = max(abs(FFT_signal)); 
fftmax = fftFreqs(freqMax);

energy = sum(abs(FFT_signal).*abs(FFT_signal));

energy = energy/(total_energy - energy);
end

