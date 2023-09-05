# code_SWD_HTE
This repository includes multiple R files that contain the functions and implementation code for all the simulation results presented in the paper [Li F, Chen X, Tian Z, Wang R, Heagerty PJ, Planning Stepped Wedge Cluster Randomized Trials to Detect Treatment Effect Heterogeneit. Under Review.]. Their descriptions are as follows. For questions or comments about the code, please contact Zizhong Tian at <zqt5121@psu.edu>.

I. Supporting Files: These supporting files are sourced in the corresponding main files that reproduce the simulation tables in the main manuscript as well as the supplementary material.

1) functions_calc_ss.R = 4 functions about sample size calculations for both cross-sectional and closed-cohort stepped wedge designs and both HTE and ATE tests;
2) functions_gendata.R = 2 functions about data generation, respectively under the cross-sectional design and the closed-cohort design;
3) functions_empirical.R = 6 functions about calculating the empirical type I error or empirical power under either cross-sectional or closed-cohort design and either HTE or ATE test (also includes 2 functions for unadjusted fitting of ATE).

II. Main Files: These main files are used to reproduce the simulation results in the main manuscript as well as the web appendix.

4) sim_CS_HTE.R = reproduce the simulation results for testing HTE under cross-sectional design (Table 1, Web Table 1);
5) sim_CC_HTE.R = reproduce the simulation results for testing HTE under closed-cohort design (Web Table 3);
6) sim_CS_ATE.R = reproduce the simulation results for testing ATE under cross-sectional design (Table 2, Web Table 2);
7) sim_CC_ATE.R = reproduce the simulation results for testing ATE under closed-cohort design (Web Table 4);
8) sim_CS_ATE_unadj.R = reproduce the comparison simulation results for testing adjusted and unadjusted ATE under cross-sectional design (Web Table 5-8);
9) sim_CC_ATE_unadj.R = reproduce the comparison simulation results for testing adjusted and unadjusted ATE under closed-cohort design (Web Table 9, 10).

III. Software 

Analyses were conducted with R, version 4.1.1 (https://www.r-project.org/)
The calculations used R packages nlme (version 3.1-152) and lme4 (version 1.1-27.1).

IV. R commands for the installation of R packages 

install.packages(c("nlme", "lme4")) 
