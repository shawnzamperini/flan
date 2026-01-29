============================================================================================
Publications
============================================================================================

The following is a list of publications that involve Flan.

- `Zamperini, S. A., Bernard, T. A., Rudakov, D. L. & Boedo, J. A. Turbulent drifts of impurity ions as an explanation for anomalous radial transport in the far-SOL of DIII-D. Nucl. Fusion 64, 074002 (2024). <https://dx.doi.org/10.1088/1741-4326/ad4c78>`_
    - Not technically Flan, but this was the proof-of-principle demonstration that a code like Flan is needed. This study used a cpython implementation that transports impurities via their guiding center motion in a Gykell background. We eventually found that the guiding center approximation could be invalidated for tungsten, which ultimately triggered the full development of Flan as a standalone C++ code optimized for modern machines. 
