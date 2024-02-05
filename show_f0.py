import matplotlib.pyplot as plt
import numpy as np
import librosa, librosa.display
from scipy.io import wavfile
from os.path import basename, join

F = 160
fs = 16000

out_dir = "./test_files/results/"
# wavpath = "audio_files/percepnet-clip-16k.wav"
wavpath = "test_files/p225_242_mic1_16k.wav"
fs1, wav_audio = wavfile.read(wavpath)

f0__ = np.fromfile(join(out_dir, basename(wavpath)[:-3]+"f0bin"), dtype=np.float32)
f0__ = np.array(f0__)
# rnnoise_f0_48k = 48000 / rnnoise_f0

fig, ax = plt.subplots(figsize=(15, 5))
times = librosa.times_like(f0__)
D = librosa.feature.melspectrogram(y=wav_audio / 32768, n_fft=2048, hop_length=F, sr=fs, n_mels=290, power=2)
D = librosa.power_to_db(D, ref=np.max)

img = librosa.display.specshow(D, x_axis='time', y_axis='mel', ax=ax)
fig.colorbar(img, ax=ax, format="%+2.f dB")

# f0 = 2595.0 * np.log10(1.0 + f0 / 700.0)
ax.plot(times, f0__, label='f0', color='cyan', linewidth=1)
ax.legend(loc='upper right')

plt.show()