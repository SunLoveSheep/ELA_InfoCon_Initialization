# ELA_InfoCon_Initialization

Started on 2016.2.27 by Mike Yi SUN.

Based on the idea of Exploratory Landscape Analysis (ELA) with Information Content (InfoCon), we would like to utilize InfoCon as lower level features for functions. These features will be used to build a prediction model to classify functions to different types according to their performance in our study of initialization resource allocation.

This is a C++ program reading from benchmark functions and making samplings on each functions to calculate their information content data. The data will be outputed to .csv files first. May also add functions to predict function classifications based on Information Content data.

Sampling methods should involve Monte Carlo and Latin Hypercube

5 tracks:
1. Main framework & data structures
2. Functions for sampling
3. Calculating information contents based on samples
4. Output data to files
5. Prediction functions
