options ls=78 NODATE;
 
data AR1;
  infile "C:\Users\arielle.thibeault\Desktop\STAT 580\Project 1\humanAR_1mbData.csv" firstobs=2 delimiter=',';
  input chr $ start end insstd delstd substd microsatstd state GC CpG nCGm LINE SINE NLp;
  run;

PROC PRINT data = AR1;
RUN;


PROC SGSCATTER DATA=AR1;                                                                                                                 
	MATRIX start end insstd delstd substd microsatstd state GC CpG nCGm LINE SINE NLp /
	DIAGONAL = (histogram normal);                                                                                                                                                                                  
RUN;
/* Subsetting our dataset to just the four mutation rates
 * of interest in our analysis.
 */

data AR1_sub;
SET AR1;
DROP chr start end state GC CpG nCGm LINE SINE NLp;
RUN;

/* Some EDA on covariance of the bivariate pairs.
 */
PROC CORR DATA = AR1_sub OUTP = pearson_corr COV;
  VAR insstd delstd substd microsatstd;
RUN;

PROC CORR DATA = AR1_sub OUTP = pearson_corr;
  VAR insstd delstd substd microsatstd;
RUN;

PROC SGSCATTER DATA=AR1_sub;                                                                                                                 
	MATRIX insstd delstd substd microsatstd /
	DIAGONAL = (histogram normal);                                                                                                                                                                                  
RUN;

/* Anderson Darling test for Normality: two variables pass,
 * the rest are close enough that empirically we will consider
 * them to be normal.
 */
PROC UNIVARIATE data=AR1_sub NORMAL;
VAR insstd delstd substd microsatstd; 
HISTOGRAM insstd delstd substd microsatstd;
RUN;

 /* The princomp procedure performs pca on the AR1_sub data.
  * The cov option specifies results are calculated from the covariance
  * matrix, instead of the default correlation matrix.
  */
 
proc princomp data=AR1_sub cov out=a;
  var insstd delstd substd microsatstd;
  run;
 
 /* The cov procedure is used to calculate pairwise correlations
  * between the first 2 principal components and the original variables.
  */
 
proc corr data=a;
  var prin1 prin2 insstd delstd substd microsatstd;
  run;
 
 /* The gplot procedure is used to plot the first 2 principal components.
  * axis1 and axis2 options set the plotting window size,
  * and these are then set to vertical and horizontal axes, respectively.
  */
 
proc gplot data=a;
  axis1 length=5 in;
  axis2 length=5 in;
  plot prin2*prin1 / vaxis=axis1 haxis=axis2;
  run;

/*Biplot*/
proc prinqual data=AR1_sub plots=(MDPref) 
              n=2       /* project onto Prin1 and Prin2 */
              mdpref=1; /* use COV scaling */
   transform identity(insstd delstd substd microsatstd);  /* identity transform */
   ods select MDPrefPlot;
run;

proc factor data = AR1_sub cov scree ev method = principal;
var insstd delstd substd microsatstd;
run;
