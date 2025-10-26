🧠 fmri-fdr-analysis

Exploring brain activation patterns through statistical inference using the Haxby fMRI dataset.

Overview

This project explores how our brain responds to visual stimuli — specifically distinguishing between everyday objects like chairs and scissors — using fMRI data from the Haxby 2001 dataset.

The goal was to uncover which brain voxels (3D pixels) show significant activation patterns and to understand how statistical correction methods like Benjamini–Hochberg (BH) and q-value help control false discoveries in such high-dimensional data.

Motivation

I’ve always been fascinated by how complex brain activity can be translated into numbers, patterns, and visuals. When I discovered the Haxby dataset, I wanted to go beyond simple visualization — I wanted to statistically test which parts of the brain truly light up during different visual tasks.

This project became my deep dive into neurostatistics — connecting data science, statistical inference, and neuroscience into one cohesive exploration.

Methodology

Dataset Used: Haxby 2001 fMRI dataset (chair vs scissors conditions).

Data Preprocessing: Reshaped 4D fMRI images into 2D matrices (timepoints × voxels).

Voxel-wise Testing: Performed two-sample t-tests for each voxel to detect significant activation differences.

Multiple Comparison Correction:

Applied Benjamini–Hochberg (BH) correction for controlling the false discovery rate (FDR).

Key Results

Both BH and q-value methods identified strong activation clusters in regions associated with visual processing and object recognition.

The q-value approach detected a few extra voxels, reflecting its slightly more adaptive thresholding compared to BH.

Overlap analysis showed high consistency between both methods — confirming robustness of detected brain regions.

Insights

This analysis shows how statistical rigor can reveal meaningful neurobiological signals. Even small differences in correction methods — like BH vs. q-value — can shape our understanding of brain activity patterns.

Beyond the results, this project taught me how to blend data visualization, high-dimensional statistics, and curiosity into a single story told through voxels.
Applied q-value estimation for a more adaptive control of FDR.

Visualization: Created voxel activation maps comparing significant voxels detected by each method and examined their overlap.

Future Work

In upcoming iterations, I plan to:

Explore knockoff-based inference for variable selection in high-dimensional data.

Extend the analysis to more conditions within the Haxby dataset.

Apply dimensionality reduction and classification methods for decoding visual stimuli.

✨ Author

Sarah Shahama V P
