# TaoSbeechDSP

本代码模块是我在学习传统 dsp 回声消除 和 降噪 积累的代码，包括个人的改造和穿插拼凑，包括：
## 基于 lpc白化+动态规划 的人声基频跟踪模块
## 基于 porprotional NLMS 的回声消除模块
## 基于 Kalman-simpified 的回声消除模块
## 基于 Kalman-original 的回声消除模块
## 基于 speexdsp 中的降噪模块，这里使用了子带上的滤波

本模块代码纯属学习自用，代码无任何优化（甚至为了阅读，大大降低了效率，例如其中自写的 泛型vector库）。
后期可以考虑

## 相关代码学习
1. speexdsp
2. webrtc-aec3
3. xmos-dsp

## 相关论文学习
1. Speech enhancement for non-stationary noise environment
2. ON SPEECH ENHANCEMENT UNDER SIGNAL PRESENCE UNCERTAINTY
3. Multidelay Block Frequency Domain Adaptive Filters
4. MICROPHONE ARRAY POST-FILTER FOR SEPARATION OF SIMULTANEOUS NON-STATIONARY SOURCES
5. Speech Enhancement Using a- Minimum Mean-Square Error Short-Time Spectral Amplitude Estimator
6. On Adjusting the learning rate in Frequency Domain Echo Cancellation With Double-Talk
7. A Novel Psychoacoustically Motivated Audio Enhancement Algorithm Preserving Background Noise Characteristics.
8. Frequency-domain adaptive Kalman ﬁlter for acoustic echo control in hands-free telephones
9. Robust acoustic echo cancellation using Kalman ﬁlter in double talk scenario
10. State-Space Frequency-Domain Adaptive Filtering for Nonlinear Acoustic Echo Cancellation
11. Block Implementation of Adaptive Digital Filters
12. Elimination of the Musical Noise Phenomenon with the Ephraim and Malah Noise Suppressor
13. Fast Implementation of LMS adaptive filters
14. A Psychoacoustic Approach to Combined Acoustic Echo Cancellation and Noise Reduction
15. Speech Enhancement Using a Soft-Decision Noise Suppression Filter

## TODO
1. NEURALKALMAN: A LEARNABLE KALMAN FILTER FOR ACOUSTIC ECHO CANCELLATION


