- [1. 15-01 radar odometry using phase Direct Radar Odometry using Phase Correlation](#1-15-01-radar-odometry-using-phase-direct-radar-odometry-using-phase-correlation)
- [2. 15-02  2010 springer book-Radar Scan Matching SLAM Using the Fourier-Mellin Transform](#2-15-02--2010-springer-book-radar-scan-matching-slam-using-the-fourier-mellin-transform)
- [3. 15-03   A blending method for forward-looking sonar mosaicing handling intra- and inter-frame artifacts](#3-15-03---a-blending-method-for-forward-looking-sonar-mosaicing-handling-intra--and-inter-frame-artifacts)
- [4. 05-03 DiSCO: Differentiable Scan ContextWith Orientation](#4-05-03-disco-differentiable-scan-contextwith-orientation)
- [5. FMT in 3d](#5-fmt-in-3d)

## 1. 15-01 radar odometry using phase Direct Radar Odometry using Phase Correlation
a direct radar odometry method is proposed
to estimate relative pose using Fourier Mellin Transform
and local graph optimization

Park et al. [65]
applied a Fourier–Mellin transform to radar scans, generating
decoupled rotation and translation estimates from log-polar and
Cartesian radar images, respectively.Acourse-to-fine phase correlation
between Cartesian images further refines the translation
estimate.
## 2. 15-02  2010 springer book-Radar Scan Matching SLAM Using the Fourier-Mellin Transform
A more recent scan-matchingbased
SLAM approach uses the Fourier-Mellin transformation
to match consecutive Radar scans, where the power spectra
are interpreted as 360° images [3].

## 3. 15-03   A blending method for forward-looking sonar mosaicing handling intra- and inter-frame artifacts
We use the Fourier–Mellin Transform to estimate the 𝐓 (Reddy and
Chatterji, 1996), as done in Hurtós et al. (2012), Hurtos et al. (2015),
Franchi et al. (2018) and Mueller et al. (2017). To the best of our
knowledge, Fourier-based registration appears to be the most stable
and precise method for FLS image registration. This method relies on
the Fourier shift theorem for translation estimation, and it can also be
employed to estimate rotation through the polar representation of the
frequency magnitude of two images.
## 4. 05-03 DiSCO: Differentiable Scan ContextWith Orientation
Following the polar transformed bird’s-eye view representation
in the Scan Context [8], our end-to-end learnable
DiSCO architecture further achieves interpretability by adopting
the Fourier-Mellin Transform after a feature extractor. It explicitly
eliminates the rotation variance for the feature extractor,
so that the extractor only needs to learn features rather than
learning features and pose estimation. This design improves
the efficiency of the feature extractor and also constrains each
module to have clear functionality. As rotation-invariance is
converted into translation-invariance in the polar domain, we
have a spectrum-based similarity evaluation metric that can be
efficiently evaluated via Euclidean distance
## 5. FMT in 3d
[25] H. Bülowand A. Birk, “Scale-free registrations in 3D: 7 degrees of freedom
with fourier mellin soft transforms,” Int. J. Comput. Vis., vol. 126, no. 7,
pp. 731–750, 2018.
[26] H. Bülow, C. A. Mueller, A. G. Chavez, F. Buda, and A. Birk, “A divide
and conquer method for 3 d registration of inhomogeneous, partially
overlapping scans with fourier mellin soft (FMS),” in Proc. IEEE Int.
Conf. Robot. Automat., 2020, pp. 8594–8601

Bülow, H.,&Birk,A. Scale-free registrations in 3d: 7 degrees of
freedom with FourierMellin soft transforms