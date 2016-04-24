libname negbin 'E:\Project Dataset\Negbin\Lkeq';
libname pearson2 'E:\Project Dataset\Negbin\Lkeq\Pearson_Chi';
data nega.lkeq;
seed=42850190;
 rate0=6;
 rate1=6;
 alpha0=log(rate0);
 alpha1=log(rate1);
 phi=0.5;
 alpha=1/phi;
 beta=phi;
 /* set number of blocks  */
 N_Block=8;
 /* set variance for random effects */
 blk_var=1;

 /* generate data for 1000 experiments  */
do expt=1 to 1000;
 do Block=1 to N_Block;
  B_j=sqrt(blk_var)*rannor(seed);  
  do treatment=1 to 4;
   eta=((treatment=1)+(treatment=2))*alpha0+((treatment=3)+(treatmetn=4))*alpha1+B_j; 
   lamda_ij=exp(eta);                               
   unit_ij=beta*rangam(seed,Alpha);
    count_ij=ranpoi(seed,lamda_ij*unit_ij);                
   output;
 end;
end;  
end;
run;

libname poiunit 'G:\Simulated SAS dataset\Poisson-Norm\Lkeq';
libname negbin'G:\Simulated SAS dataset\Nega-Norm\Lkeq';
libname pearson 'G:\Simulated SAS dataset\Poisson-Norm\Pearson_Chi';
libname estimt 'G:\Simulated SAS dataset\Poisson-Norm\Estimators';
ods results off;
ods html exclude all;
%macro toexcel(dataset);
%do i=1 %to &dataset;
ods csv body="G:\Simulated dataset\Nega-Norm\Lkeq\&i..csv";
proc print data=negbin.lkeq_&i;
var expt block treatment count_ij;
run;
ods csv close;
%end;
%mend;
%toexcel(1000);


