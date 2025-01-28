# Effector-Dependent Neural Representations of Perceptual Decisions Independent of Motor Actions and Sensory Modalities

## Authors: 
- Marlon F. Esmeyer 
- Timo T. Schmidt
- Felix Blankenburg

## Abstract 
Neuroscientific research has shown that perceptual decision-making occurs in brain regions that are associated with the required motor response. Recent functional magnetic resonance imaging (fMRI) studies that dissociated decisions from coinciding processes, such as the motor response partly challenge this, indicating that perceptual decisions are represented in an abstract or sensory-specific manner that might vary across sensory modalities. However, comparisons across sensory modalities have been difficult since most task designs differ not only in modality but also in effectors, motor response, and level of abstraction. Here, we describe an fMRI experiment where participants compared frequencies of two sequentially presented visual flicker stimuli in a delayed match-to-comparison task, which controlled for motor responses and stimulus sequence. A whole-brain searchlight support vector machine analysis of multi voxel patterns was used to identify brain regions containing information on perceptual decisions. Furthermore, a conjunction analysis with data from a re-analyzed analogue vibrotactile study was conducted for a comparison between visual and tactile decision-making processes. Both analyses revealed significant activation patterns in the left dorsal premotor cortex (PMd) as well as in bilateral parts of the posterior parietal cortex (PPC). While previous primate and human imaging research have implicated these regions in transforming sensory information into action, our findings indicate that the PPC processes abstract decision signals while PMd represents an effector-dependent, but motor response independent encoding of perceptual decisions that persists across sensory domains.

## Scripts
Matlab codes for the analyses performed for the manuscript titled "Effector-Dependent Neural Representations of Perceptual Decisions Independent of Motor Actions and Sensory Modalities".

### fMRI Analyses
#### Main Analyses
- Decoding_Batch: Batch script running all pre-processing and analysis functions (first level)
- D1_glm_1stLevel_SVM: Function running the GLM
- D2_Decoding_SVM: Function running the support vector machine decoding
- Second_Level_Decoding_SVM: Script running the 2nd level one-sample t-test

#### Conjunction Analysis
- t_test_SVM: Batch script for the conjunction
- t_test_job_SVM: setting up the conjunction
- Estimate_job_SVM: Parameter estimation for the conjunction

#### Control Analyses
- Decoding_Batch_control: Batch for the Support Vector Machine classification for motor response and rule decoding
- Decoding_Batch_subsampling: Batch for the Support Vector Machine classification for the subsampling of motor response and task rule
- D1_glm_1stLevel_left_vs_right: Function running the GLM for motor response
- D1_glm_1stLevel_rule: Function running the GLM for rule
- D1_glm_1stLevel_sub_motor_rule: Function running the GLM for the subsampling of motor response and task rule
- D2class_Decoding: Function running the support vector machine classification
- Second_Level_Decoding_left_vs_right: Script running the 2nd level one-sample t-test for motor response
- Second_Level_Decoding_rule: Script running the 2nd level one-sample t-test for rule
- Second_Level_Decoding_sub_motor_rule: Script running the 2nd level one-sample t-test for the subsampling

#### Behavioural Analyses
- ANOVA_three_way: Script running a three-way repeated measures ANOVA (factors: rule, stimulus order, f1 frequency)
- Chi_Squared_tests: Script running the chi-squared testd (comparing motor responses)

### Required software packages and toolboxes: 
- SPM12: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
- The Decoding Toolbox (TDT): https://sites.google.com/site/tdtdecodingtoolbox/