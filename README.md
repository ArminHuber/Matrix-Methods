# Matrix-Methods
Calculate guided wave dispersion diagrams for free, perfectly elastic, multilayered, anisotropic plates by using the transfer matrix method (TMM) or the stiffness matrix method (SMM). Using TMM, the numerical instability at high frequency-thickness products becomes obvious. By contrast, SMM remains stable. When computing SH waves (in decoupled cases), TMM is better than SMM because SMM does not find all modes and the numerical instability of TMM does not affect the guided SH waves.

Written in MATLAB R2022a.

## Files
* `MatrixMethods.m`: The general case of wave propagation along arbitrary directions in layups with arbitrary layer orientations. Motion in the sagittal plane and in the shear horizontal direction is coupled.
* `MatrixMethods_Decoupled_Lamb.m`: The decoupled case of wave propagation along axes of symmetry in single layers or in cross ply laminates. This code computes motion in the sagittal plane (pure Lamb waves).
* `MatrixMethods_Decoupled_SH.m`: Same as above, but this code computes motion in the shear horizontal direction (pure SH waves).

## Literature

### Books
* A. H. Nayfeh, *Wave Propagation in Layered Anisotropic Media with Applications to Composites* (North-Holland, Amsterdam, 1995).
* S. I. Rokhlin, D. E. Chimenti and P. B. Nagy, *Physical Ultrasonics of Composites* (Oxford University Press, Oxford, 2011).

### Journal articles
* A. H. Nayfeh, "The propagation of horizontally polarized shear waves in multilayered anisotropic media," [J. Acoust. Soc. Am.](https://doi.org/10.1121/1.398580) **86**(5), 2007–2012 (1989).
* A. H. Nayfeh, "The general problem of elastic wave propagation in multilayered anisotropic media," [J. Acoust. Soc. Am.](https://doi.org/10.1121/1.400988) **89**(4), 1521–1531 (1991).
* S. I. Rokhlin and L. Wang, "Stable recursive algorithm for elastic wave propagation in layered anisotropic media: Stiffness matrix method," [J. Acoust. Soc. Am.](https://doi.org/10.1121/1.1497365) **112**(3), 822-834 (2002).
* L. Wang and S. I. Rokhlin, "Stable reformulation of transfer matrix method for wave propagation in layered anisotropic media," [Ultrasonics](https://doi.org/10.1016/S0041-624X(01)00082-8) **39**, 413-424 (2001).
* V. G. A. Kamal and V. Giurgiutiu, "Stiffness transfer matrix method (STMM) for stable dispersion curves solution in anisotropic composites," Proc. SPIE **9064** (2014).
* A. M. A. Huber and M. G. R. Sause, "Classification of solutions for guided waves in anisotropic composites with large numbers of layers," [J. Acoust. Soc. Am.](https://doi.org/10.1121/1.5082299) **144**(6), 3236-3251 (2018).
* A. M. A. Huber, "Classification of solutions for guided waves in fluid-loaded viscoelastic composites with large numbers of layers," [J. Acoust. Soc. Am.](https://doi.org/10.1121/10.0020584) **154**(2), 1073–1094 (2023).

### PhD thesis
* A. M. A. Huber, *Numerical Modeling of Guided Waves in Anisotropic Composites with Application to Air-coupled Ultrasonic Inspection* ([University of Augsburg](https://opus.bibliothek.uni-augsburg.de/opus4/frontdoor/index/index/year/2021/docId/82760), Augsburg, 2020).
