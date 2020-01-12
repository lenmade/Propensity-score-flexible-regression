Files in this folder include:

1. main_psReg.R: main code to perform the flexible regression approach with different propensity score weighting methods;
2. func_psReg.R: supporting function for main_psReg.R;
3. func_glimmix.txt: SAS code to estimate penalty parameter for spline fitting using PROC GLIMMIX;
4. sampe_data.csv: a data example;
5. data_Z0.csv, data_Z1.csv: data generate from sample_data.csv to be used in func_glimmix.txt to calculate penalty parameters for treated and control groups, respectivley.