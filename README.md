DK_fitting Protocol

Computes the diffusion tensor and kurtosis tensor values for a data set. Results are saved into the file name specified. Final results include lambda1, lambda2, lambda3, vec1, vec2, vec3, K1, K2, K3, MK, exitLog, resLog, R2Log.

Files Needed:
•	Diffusion weighted images - .rec (float)
•	Gradient direction file - .txt


Output Files
•	.mat – all DTI and DKI metrics
•	.axdiff – axial diffusion
•	.raddiff – radial diffusion
•	.md – mean diffusion
•	.fa – fractional anisotropy
•	.axkurt – axial kurtosis
•	.radkurt – radial kurtosis
•	.mk – mean kurtosis

Locate and Place Necessary Files into Folders
1.	Open Matlab
2.	Go to directory with current DKI version
3.	Place registered DWI into folder recFiles
54.	Place gradient file into folder gradFiles


Create DKI Options File
1.	Open genDKIopt.m
2.	Enter in option parameters
a.	numMatPool: number of Matlab Pools to open for parallelization
b.	dim: dimensions of file (w,h,slice,grad)
c.	f1name: first file name (.rec float)
i.	add f2name, f3name, etc. for additional files
d.	gradientFile: gradient file name (.txt)
e.	b_value1: first file b-value
i.	add b_value2, b_value3, etc. for additioni b values
g.	S0thresh: intensity threshold for zero-gradient image
h.	gradZeroIdx: specifies which gradient index is the zero gradient
i.	outputFile: name of the output file (.mat)

options.numMatPool = 7;
options.dim = [256 256 65 32];
options.numB = 2;
options.f1name = 'recFiles/3T5324_1000.rec';
options.f2name = 'recFiles/3T5324_2500.rec';
options.gradientFile = 'gradFiles/Jones30_gradient.txt';
options.b_value1 = 1000;
options.b_value2 = 2500;
options.S0thresh = 10;
options.gradZeroIdx = 31;
options.outputFile = 'outputFiles/DKinfo_adult1.mat';
save('optFiles/DKIopt_adult1.mat','options')

3.	Run genDKIopt – press F5 

Run DKI Calculation
1.	Call DKI_compute

DKI_compute('optFiles/DKIopt_adult1.mat');
