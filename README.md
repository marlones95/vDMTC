# Effector-Specific Neural Representations of Perceptual Decisions Independent of Motor Actions and Sensory Modalities

## Authors: 
- Marlon F. Esmeyer 
- Timo T. Schmidt
- Felix Blankenburg

## Abstract 
Neuroscientific research has shown that perceptual decision-making occurs in effector-specific brain regions that are associated with the required motor response. Recent functional magnetic resonance imaging (fMRI) studies that dissociated decisions from coinciding processes, such as motor actions partly challenge this, indicating abstract representations that might vary across stimulus modalities. However, cross-modal comparisons have been difficult since most task designs differ not only in modality but also in effectors, motor response, and level of abstraction. Here, we describe an fMRI experiment where participants compared frequencies of two sequentially presented visual flicker stimuli in a delayed match-to-comparison task, which controlled for motor actions and stimulus sequence. Using Bayesian modelling, we estimated subjective frequency differences based on the time order effect. These values were applied in support vector regression analysis of a multi-voxel pattern whole-brain searchlight approach to identify brain regions containing information on subjective decision values. Furthermore, a conjunction analysis with data from a re-analyzed analogue vibrotactile study was conducted for a cross-modal comparison. Both analyses revealed significant activation patterns in the left dorsal (PMd) and ventral (PMv) premotor cortex as well as in the bilateral intraparietal sulcus (IPS). While previous primate and human imaging research have implicated these regions in transforming sensory information into action, our findings indicate that the IPS processes abstract decision signals while PMd and PMv represent an effector-specific, but motor response independent encoding of perceptual decisions that persists across sensory domains.

## Scripts
Matlab codes for the analyses performed for the manuscript titled "Effector-Specific Neural Representations of Perceptual Decisions Independent of Motor Actions and Sensory Modalities".

### VBA
- VBA_MAIN: Batch script for the variational Bayes routine
- f_vDMTC_simple: Evolution function
- g_vDMTC_simple_bias: Observation function
- Plot_TOE: Illustration of the time-order effect
- Figure_2b_c: Script for Figures 2b and 2c
- Figure 2D: Script for Figure 2D
- Supplementary_Figure_1: Script for Supplementary Figure 1
- compute_TOE_performance: Function computing d' and empirical time-order effect
- Model_comparisons: Script performing single subject model comparison and group level Bayesian model selection

### fMRI Analyses
#### Main Analyses
- Decoding_Batch_SVR: Batch script running all pre-processing and analysis functions (first level)
- D1_glm_1stLevel_SVR: Function running the GLM
- D2_Decoding_SVR: Function running the support vector regression decoding
- Second_Level_Decoding: Script running the 2nd level one-sample t-test

### Required software packages and toolboxes: 
- SPM12: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
- The Decoding Toolbox (TDT): https://sites.google.com/site/tdtdecodingtoolbox/
- VBA toolbox: https://mbb-team.github.io/VBA-toolbox/
