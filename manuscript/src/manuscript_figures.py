import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

np.random.seed(123456)      # for reproducibility

# computational time complexity sample
n = np.arange(1,101,1)      # number of terms
O_2N2 = 2*n**2              # Discrete Fourier Transform
O_2NlgN = 2*n*np.log2(n)    # Fast-Fourier Transform

# decomposition of a signal
L = 10                          # period time length, s
Ts = 1/50                       # sampling frequency, 1/s
t = np.arange(0,L-Ts,Ts)        # time, s

a1 = 2                          # amplitude
w1 = 5                          # angular frequency in sig1
sig1 = a1*np.sin(2*np.pi*w1*t)  # signal 1

a2 = 1                          # amplitude
w2 = 20                         # angular frequency in sig2
sig2 = a2*np.sin(2*np.pi*w2*t)  # signal 2

# combine the two signals
signal = sig1 + sig2

# discrete Fourier transform via FFT
fft = np.fft.fft(signal)
# frequency vector for signal sampling in frequency space
fs = 1/Ts
frq = np.arange(0,len(fft),1)*fs/len(fft)

# center shift the fft calculated above
fftshftd = np.fft.fftshift(fft)
frqshftd = np.arange(-len(signal)/2,len(signal)/2,1)*fs/len(signal)
# power spectral density (psd)
psd = np.abs(fftshftd)**2
# normalize psd
psd_norm = psd/np.max(psd)

# create a noisy signal
sig_noisy = signal*np.random.rand(len(signal))

# discrete Fourier transform via FFT
fftnoisy = np.fft.fft(sig_noisy)
# center shift the fft
fftnoisyshftd = np.fft.fftshift(fftnoisy)
# power spectral density (psd)
psd_noisy = np.abs(fftnoisyshftd)**2
# normalize psd
psd_noisy_norm = psd_noisy/np.max(psd_noisy)

# create a 2D black and white image based on the signals above
synthetic_img = np.zeros([len(signal),len(signal)])
for i in np.arange(0,len(signal),1):
    for j in np.arange(0,len(signal),1):
        synthetic_img[i,j] = signal[i]*np.random.rand()*j/255

# 2D FFT of the image
fft2 =np.fft.fft2(synthetic_img)

# center shift the fft
fft2shftd = np.fft.fftshift(fft2)
# 2D power spectral density
psd2D = np.abs(fft2shftd)**2
# normalize the 2D psd
psd2D_norm = psd2D/np.max(psd2D)

#############################################################################################
# Surface rougheness map
sq = [0.0634, 0.0437, 0.0349]
Rsk = [1.901, 2.037, 2.094]
Rku = [4.69, 5.83, 6.602]

# in vivo samples
# https://doi.org/10.1111/2041-210X.12778
# [Red maple leaf (Acer rubrum), Back of hand (Homo sapiens), ]
Rsk_vivo = [0.42, -0.19, -0.04]
Rku_vivo = [4.3, 3.5, 2.4]

x = np.arange(-3, 3, 0.01)
y = x**2 + 1

#############################################################################################
# visualization
# format sizing constants
TICKSFONT = 13; TITLEFONT = 15; TEXTFONT = 15   # fonts
FIGWIDTH = 6.4; FIGHEIGHT = 4.8                 # figure size
LINEWIDTH = 3; ROLLWINDOW = 100                 # plot spec
CLRS = list(mcolors.TABLEAU_COLORS)             # list of tableau colors

# signal_decomposed
NROWS = 1; NCOLS = 1
f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))

axs.plot(t, sig_noisy, label="Noisy signal = Raw signal * random_noise")
axs.plot(t, sig1, label="Signal 1")
axs.plot(t, sig2, label="Signal 2")
axs.plot(t, signal, label="Raw signal = Signal 1 + Signal 2")
axs.set_xlabel("Time, s")
axs.set_ylabel("Amplitude")
axs.legend()
axs.set_ylim([-4,8])
axs.set_xlim([0,0.5])
# axs.annotate("Raw signal = signal1 + signal2", xy=(0.9,0.1), xycoords='axes fraction')

f.tight_layout()
f.savefig("figures/signal_decomposed.png")

# noisy_signal
NROWS = 1; NCOLS = 1
f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))

axs.plot(t, sig_noisy)
axs.set_xlabel("Time, s")
axs.set_ylabel("Amplitude")
axs.set_xlim([0,4])

f.tight_layout()
f.savefig("figures/noisy_signal.png")

# normalized power spectral density of noisy_signal
NROWS = 1; NCOLS = 1
f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))

axs.plot(frq,psd_noisy_norm)
axs.set_xlabel("Frequency, Hz")
axs.set_ylabel("Normalized Power, AU")

f.tight_layout()
f.savefig("figures/normalized_power_spectrum.png")

# fft_summary
NROWS = 9; NCOLS = 2
f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))
axs[0,0].plot(t, signal)
axs[0,0].set_title("Original signal. signal = sig1 + sig2")
axs[0,0].set_xlabel("Time, s")
axs[0,0].set_ylabel("Amplitude")

axs[0,1].axis('off')

axs[1,0].plot(frq, np.abs(fft))
axs[1,0].set_title("Magnitude, mag = abs(fft)")
axs[1,0].set_xlabel("Frequency, Hz")
axs[1,0].set_ylabel("Magnitude")

axs[1,1].plot(frqshftd, np.abs(fftshftd))
axs[1,1].set_title("Center shifted-nomalized magnitude, mag = abs(fftshift(fft))")
axs[1,1].set_xlabel("Frequency, Hz")
axs[1,1].set_ylabel("Magnitude")

axs[2,0].plot(frq, psd)
axs[2,0].set_title("Power spectral density from original signal")
axs[2,0].set_xlabel("Frequency, Hz")
axs[2,0].set_ylabel("Power, AU")

axs[2,1].plot(frq, psd_norm)
axs[2,1].set_title("Normalized power spectral density from original signal")
axs[2,1].set_xlabel("Frequency, Hz")
axs[2,1].set_ylabel("Normalized Power, AU")

axs[3,0].plot(t, sig_noisy)
axs[3,0].set_title("Noisy signal, signal*rand()")
axs[3,0].set_xlabel("Time, s")
axs[3,0].set_ylabel("Amplitude")

axs[3,1].axis('off')

axs[4,0].plot(frqshftd, np.abs(fftnoisy))
axs[4,0].set_title("Magnitude, mag = abs(fft)")
axs[4,0].set_xlabel("Frequency, Hz")
axs[4,0].set_ylabel("Magnitude")

axs[4,1].plot(frqshftd, np.abs(fftnoisyshftd))
axs[4,1].set_title("Center shifted-normalized magnitude, mag = abs(fftshift(fft))")
axs[4,1].set_xlabel("Frequency, Hz")
axs[4,1].set_ylabel("Magnitude")

axs[5,0].plot(frq, psd_noisy)
axs[5,0].set_title("Power spectral density from noisy signal")
axs[5,0].set_xlabel("Frequency, Hz")
axs[5,0].set_ylabel("Power")

axs[5,1].plot(frq, psd_noisy_norm)
axs[5,1].set_title("Normalized power spectral density from noisy signal")
axs[5,1].set_xlabel("Frequency, Hz")
axs[5,1].set_ylabel("Normalized Power, AU")

axs[6,0].imshow(synthetic_img, cmap='gray')

axs[6,1].axis('off')
# axs[6,1].remove()
# f.add_subplot(NROWS, NCOLS, 14, projection='3d')
# for zi in z:
#     axs[6,1].plot(xs=z,ys=psd2D[zi],zs=zi,zdir='y')

axs[7,0].plot(np.abs(fft2))
axs[7,0].set_title("Magnitude fft2")

axs[7,1].imshow(np.abs(fft2shftd))
axs[7,1].set_title("Magnitude fft2_shifted")

axs[8,0].imshow(np.log(psd2D))
axs[8,0].set_title("Power spectral density of ff2_shifted")

axs[8,1].imshow(np.log(psd2D/np.max(psd2D)))
axs[8,1].set_title("Normalized power spectral density of ff2_shifted")

f.tight_layout()
f.savefig("figures/fft_summary.png")

#############################################################################################
# Manuscrript figures

# Figure 1
NROWS = 1; NCOLS = 1
f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))

axs.plot(n, O_2N2, label="$2N^2$", color=CLRS[0])
axs.plot(n, O_2NlgN, label="$2Nlog_2(N)$", color=CLRS[1])
# axs.set_title("")
axs.set_ylabel("Number of Operations", fontsize=TITLEFONT)
axs.set_xlabel("Number of terms, N", fontsize=TITLEFONT)
ytick_values = plt.gca().get_yticks()
axs.set_yticklabels(['{:,.0f}'.format(x) for x in ytick_values])
axs.legend()
axs.tick_params(axis='both', which='major', labelsize=TICKSFONT)

# inset setup
x1, x2, y1, y2 = 0, 10, 0, 300
inset_axs = axs.inset_axes([0.125, 0.4, 0.4, 0.4],
                        xlim=[x1, x2], ylim=[y1, y2],
                        xticklabels=[str(x1),"","","","",str(x2)], yticklabels=[str(y1),"","",str(y2)])
# plot the inset
for ax in axs, inset_axs:
    ax.plot(n, O_2N2, color=CLRS[0])
    ax.plot(n, O_2NlgN, color=CLRS[1])

axs.indicate_inset_zoom(inset_axs, edgecolor="black")

f.tight_layout()
f.savefig("figures/1.computational_complexity.png")

# Figure 2
# noisy_signal and decomposed in the inset
NROWS = 1; NCOLS = 1
f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))

axs.plot(t, sig_noisy, label="Noisy raw signal")
axs.plot([],[], label="Signal 1")
axs.plot([],[], label="Signal 2")
axs.set_xlabel("Time, s", fontsize=TITLEFONT)
axs.set_ylabel("Amplitude", fontsize=TITLEFONT)
axs.set_ylim([-3,8])
axs.legend()
axs.tick_params(axis='both', which='major', labelsize=TICKSFONT)

# inset plot
inset_axs = f.add_axes([0.2, 0.56, 0.3, 0.3])
inset_axs.plot(t, sig_noisy, color=CLRS[0])
inset_axs.plot(t, sig1, color=CLRS[1])
inset_axs.plot(t, sig2, color=CLRS[2])
inset_axs.set_xlim([0,0.5])
    
# f.tight_layout()
f.savefig("figures/2.noisy_signal_decomposed_inset.png")

# Figure 3
# (a) normalized magnitude center-shifted FFT
NROWS = 1; NCOLS = 1
f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))

axs.plot(t, sig_noisy)
axs.set_title("(a)", loc='left', fontsize=TITLEFONT)
axs.set_xlabel("Time, s", fontsize=TITLEFONT)
axs.set_ylabel("Amplitude", fontsize=TITLEFONT)
axs.tick_params(axis='both', which='major', labelsize=TICKSFONT)

f.tight_layout()
f.savefig("figures/3a.original_noisy_raw_signal.png")

# (b) magnitude center-shifted FFT
NROWS = 1; NCOLS = 1
f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))

axs.plot(frqshftd, np.abs(fftnoisyshftd))
axs.set_title("(b)", loc='left', fontsize=TITLEFONT)
axs.set_xlabel("Frequency, Hz", fontsize=TITLEFONT)
axs.set_ylabel("Magnitude", fontsize=TITLEFONT)
axs.tick_params(axis='both', which='major', labelsize=TICKSFONT)

f.tight_layout()
f.savefig("figures/3b.normalized_magnitude_center-shifted_fft.png")

# (c) normalized magnitude center-shifted FFT
NROWS = 1; NCOLS = 1
f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))

axs.plot(frqshftd, psd_noisy_norm)
# justify title to the left
axs.set_title("(c)", loc='left', fontsize=TITLEFONT)
axs.set_xlabel("Frequency, Hz", fontsize=TITLEFONT)
axs.set_ylabel("Normalized Power, AU", fontsize=TITLEFONT)
axs.tick_params(axis='both', which='major', labelsize=TICKSFONT)

f.tight_layout()
f.savefig("figures/3c.normalized_power_spctral_density.png")

# Surface roughness map
NROWS = 1; NCOLS = 1

f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))
axs.scatter(Rsk, Rku, marker='^')
for i in np.arange(0, len(Rsk_vivo), 1):
    axs.scatter(Rsk_vivo[i], Rku_vivo[i], marker='o', color=CLRS[i+1])
    
axs.plot(x,y,'k--', lw=1)
axs.set_xlabel(r'Skewness, $R_{sk}$', fontsize=TITLEFONT)
axs.set_ylabel(r'Kurtosis, $R_{ku}$', fontsize=TITLEFONT)
axs.set_xlim([-3,3])
axs.set_ylim([0,10])

f.tight_layout()
f.savefig("figures/11.surface_roughness_map.png")

from mpl_toolkits.mplot3d import Axes3D

# Parameters for the 3D surface
x = np.linspace(1, 5, 50)
y = np.linspace(1, 5, 50)
X, Y = np.meshgrid(x, y)
Z = np.sin(np.sqrt(X**2 + Y**2))  # Curved surface
z = X*0 - 2

# Create a figure for the 3D plot
f = plt.figure(figsize=(FIGWIDTH, FIGHEIGHT))
ax = f.add_subplot(111, projection='3d')

# 3D surface plot
ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none', alpha=0.8)
ax.plot_surface(X, Y, z, cmap='gist_gray' , alpha=0.5)
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_zlabel("$z$")

# Adjust the z-axis to make space for the projection
ax.set_zlim(np.min(Z) - 1, np.max(Z))

# f.tight_layout()
f.savefig("figures/3D_surface_with_projection.png")