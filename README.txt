# GOBI (General ODE-based causal inference)
This is a matlab code for accurate and broadly applicable causal inference method for time-series data.

## System requirements
1. The package is based on MATLAB R2021a
2. The package is tested on macOS M1.
3. The package does not need any non-standard hardware.

## Installation guide
1. The codes can be run directly without installation.
2. No install time is needed.


# Code Description
# 1. Tutorial

## 1. Generate the ‘data.mat’ file, which contains two variables ‘t’ and ‘y’. ‘t’ is the time points at which the measurements were taken. Each column of ‘y’ should be the data for each variable at the respective time points (see Input in Supplementary Fig. 14).

## 2. Run the ‘Step0_interpolation_and_cut.m’ file. (see Interpolation & Cut in Supplementary Fig. 14). This function interpolates the data based on the interpolation method specified by the user and cuts the data into windows using moving window technique, where the window size and overlapping ratio are defined by the users as follows.

(a) Users need to specify the interpolation method using the parameter ‘method’. Specifically, ‘method = 1’ indicates ‘linear’ interpolation, ‘method = 2’ indicates ‘spline’ interpolation, and ‘method = 3’ indicates ‘fourier’ interpolation. In the case of ‘method = 3’, users also have to specify the order of ‘fourier’ interpolation (i.e., ‘num_fourier = 1 to 8’). For less noisy data, ‘spline’ method is recommended, and for highly noisy data, ‘fourier’ method with order 2 is recommended. 

(b) Users need to choose the sampling rate for the interpolation. The parameter ‘time_interval’ indicates how finely the users wants to interpolates the original time series. For example, ‘time_interval = 0.5’ indicates interpolation using a time interval twice as fine as the original time series. Selecting ‘time_interval’ to make about 100 time points per period is recommended, and please note that the low value of ‘time_interval’ (high sampling rate) makes the inference accurate, but slow as well.

(c) For the data segmentation, users need to specify the parameters for the moving window technique, i.e., window size and overlapping ratio. The parameter ‘window_size_ori’ defines the number of time points in each window. Then, along the time series, we move the window until the next window overlaps with the current window by the ratio defined in the parameter ‘overlapping_ratio’ ('overlapping_ratio = 0.1’ as a default). For oscillatory time-series data, it is recommended to choose the window size as one period. The time series in every window is saved at the variable ‘y_total’.

After interpolation and data segmentation, the data is saved in ‘data_cut.mat’ file.

## 3. Updates the ‘Step0_options.m’ file (see Options in Supplementary Fig. 14). This code integrates options that can be adjusted by the users or set via our guidelines.

(a) Users need to specify the thresholds for regulation-detection region (‘thres_R’), regulation-detection score (‘thres_S’), and total regulation score (‘thres_TRS’) as well as the critical values for ∆ test (‘p_delta = 0.01’ as defaults), and surrogate test (‘p_surrogate = 0.001’ as defaults). To assist users in selecting those values, we have provided guidelines based on the noise level of data (Supplementary Fig. 4). Thus, users can use our guidelines as defaults or make adjustments depending on whether the goal is to decrease false positive or negative predictions. 

(b) The parameter ‘type_self’ defines options for the types of self-regulation: no assumption (‘type_self = NaN’); negative self-regulation (‘type_self = -1’); no self-regulation (‘type_self = 0’); and positive self-regulation (‘type_self = 1’). Also, users can optionally incorporate another available prior knowledge into the inference (Supplementary Fig. 3). Of course, it is possible to run the inference without any prior knowledge, but incorporating such knowledge is beneficial when the amount of data is limited. 


All the options are saved atin the ‘data_with_options.mat’ file. Also, during the inference, our framework automatically gives a warning signal when the data is insufficient to run the framework. Then, users should adjust these options.

## 4. Run the codes for 1D framework (see Step 1 in Supplementary Fig. 14).

(a) Run the ‘Step1_compute_RDS_dim1.m’ function. First, this function finds all the possible 1D regulations and saves them at the variable ‘component_list_dim1’. Each row of component_list_dim1’ indicates the set of causal variable (C) and target variable (T). For each pair (C and T), regulation-detection region and score are computed for all the regulation types (+ and -) using time-series data, and they are saved at the variables ‘R_total_list’ and ‘S_total_list’. Those values are saved in the ‘RDS_dim1.mat’ file.

(b) Run the ‘Step1_compute_TRS_dim1.m’ function. Using the ‘thres_R’ and ‘thres_S’ that users specified, Total Regulation Score (TRS) is computed for each possible 1D regulation. As a result, the heatmap of TRS is displayed, and the exact values of TRS are saved at the variable ‘TRS_total’ in the ‘TRS_dim1.mat’ file. In this heatmap, each row indicates the possible 1D regulation (C and T) and each column indicates the regulation type (+ and -). Using the ‘thres_TRS’ that users specified, 1D regulations are inferred.

(c) Run the ‘Check1D.m’ function. This function checks whether the data is sufficient to confidently infer 1D regulations. If the warning signal comes out, then users are recommended to either stop the inference or adjust the options.

## 5. Run the codes for 2D framework (see Step 2 in Supplementary Fig. 14).

(a) Run the ‘Step2_compute_RDS_dim2.m’ function. First, this function finds all the possible 2D regulations and saves them at the variable ‘component_list_dim2’. Each row of component_list_dim2’ indicates the set of two causal variables (C_1 and C_2), and target variable (T). For each set (C_1, C_2, T), regulation-detection region and score are computed using time-series data for all the regulation types ((+,+), (+,-), (-,+), and (-,-)) and they are saved at the variables ‘R_total_list’ and ‘S_total_list’. Theose values are saved in the ‘RDS_dim2.mat’ file.

(b) Run the ‘Step2_compute_TRS_dim2.m’ function. Using the ‘thres_R’ and ‘thres_S’ that users specified, Total Regulation Score (TRS) is computed for each possible 2D regulation. As a result, the heatmap of TRS is displayed, and the exact values of TRS are saved at the variable ‘TRS_total’ in the ‘TRS_dim2.mat’ file. In this heatmap, each row indicates the possible 2D regulation (C_1,C_2,T) and each column indicates the regulation type ((+,+), (+,-), (-,+), and (-,-)). Using the ‘thres_TRS’ that users specified, candidates for 2D regulations are inferred.

(c) Run the ‘Step2_Delta_test_dim2.m’ function. For every candidate for 2D regulations (Inferred from 5-(b)), this function performs the ∆ test for each causal variable (C_1 and C_2). If the number of data is smaller than 25, then this function tests whether the signs of regulation-delta functions are non-negative or not. If the number of data is larger than 25, this function performs the Wilcoxon signed ranked test. The result of ∆ test is saved at the variable ‘delta_list’ in the ‘Delta_dim2.mat’ file. Each row of the ‘delta_list’ represents the candidate for 2D regulation, and two columns of ‘delta_list’ represents the results of ∆ test of for C_1 and C_2. Using the ‘p_delta’ that users specified, candidates for 2D regulations are inferred.

(d) Run the ‘Step2_Surrogate_test_dim2.m’ function. For every candidate for 2D regulations (inferred from 5-(c)), this function performs the surrogate test for each causal variable. Users need to specify the number of bootstrapping (‘num_boot = 100’ as defaults) for the surrogate test. During the simulation, for every time-series data, one of causal variables (C_1 or C_2) is shuffled ‘num_boot’ times, and regulation-detection scores are computed. Then the p-value is computed for each data and causal variable. Those p-values are combined using Fisher’s method. The results of the surrogate test are saved in the variable ‘surrogate_list’ in the ‘Surrogate_dim2.mat’ file. Each row of ‘surrogate_list’ represents the candidate for 2D regulations. The first two columns of ‘surrogate_list’ represents the results of surrogate test for C_1 and C_2. The third and fourth columns of ‘surrogate_list’ represents the thresholds for combined p-values (combine ‘p_surrogate’ for all the data). Finally, by using these thresholds, 2D regulations are inferred.

(e) Run the ‘Check2D.m’ function. This function checks whether the data is sufficient to confidently infer 2D regulations. If the warning signal comes out, then users are recommended to either stop the inference or adjust the options.

Theose steps are continued until the ‘max_D’-dimensional framework. After that, run the function ‘Merging_regulations.m’ to infer a network structure by merging all the inferred regulations. Since GOBI involves multi-dimensional inferences, it is possible to detect various dimensional regulations for a single target. In this case, GOBI infers the regulation with the highest value of TRS.  Here, we illustrate up to the 2D framework, but users can easily expand this approach to include higher dimensions as needed (see Github codes for Fig. 4). 


# 4. Figures
The matlab codes used for drawing all the figures in the paper.
