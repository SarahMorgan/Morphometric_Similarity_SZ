# Morphometric_Similarity_SZ

This repo contains code to analyse morphometric similarity matrices from patients with schizophrenia and healthy control subjects, as reported in... Please cite that paper if you find any of this code helpful in your own analyses. The code is written in [MATLAB](https://uk.mathworks.com/products/matlab.html).

## Importing the data:

Information about where the data is stored and how to import it.

```
nregs=308; % number of regions
nsubs=115; % number of subjects
```

## Calculate the morphometric similarity matrices:

```
```

## Global differences in morphometric similarity:

## Regional differences in morphometric similarity:

The following code calculates a t-statistic for regional differences in morphometric similarity:

```
clear meanMSN_regional
for subj=1:nsubs
    meanMSN_regional(subj,:)=sum(subj_MSN_7{1,subj})./(nregs-1);
end

dummy=meanMSN_regional;

clear mytstat mypval
for region=1:nregs
  tbl = table(age(mask_conpat),gender(mask_conpat),group,dummy(mask_conpat,region));
  tbl.Var2 = categorical(tbl.Var2);
  tbl.group = categorical(tbl.group);
  lm = fitlm(tbl,'Var4~Var1*Var2+group');
  mytstat(region)=lm.Coefficients{4,3};
  mypval(region)=lm.Coefficients{4,4};
end

dlmwrite('mytstat_Maast.dat',mytstat)
dlmwrite('mypval_Maast.dat',mypval)
```

This code was run to calclulate t-statistics and p-values for all three datasets (mytstat_Maast, mytstat_Dublin, mytstat_Cobre, mypval_Maast, mypval_Dublin and mypval_Cobre). The p-values were then combined using Fisher's method:

```
clear pcomb
for region=1:nregs
  pcomb(region)=pfast([mypval_Maast,mypval_Dublin,mypval_Cobre]);
end

pvalue_fdr = mafdr(mypval,'BHFDR',1); % FDR corrected p-values
sigregs=find(pvalue_fdr<0.05) % list of the statistically significant regions

meantstat = (mytstat_Maast + mytstat_Dublin + mytstat_Cobre)./3;

```

## Gene expression PLS:

## Gene enrichments:

