This file provides code to perform the main morphometric similarity analyses (essentially Figures 1 and 2 of the paper).

## Importing the data:

Information about where the data is stored and how to import it.

Regional values for cortical thickness (CT), grey matter volume (GM), surface area (SA), mean curvature (MC), Gaussian curvature (GC), fractional anisotropy (FA) and mean diffusivity (MD) are provided for each dataset and can be loaded in MATLAB using commands such as ```load('PARC500_CT.dat')```. For each dataset the files 'group.dat', 'age.dat' and 'sex.dat' provide information about the group which the subject belongs to (control subject=1, patient=2), their age and their sex (male=1, female=2).

```
load('PARC500_GC.dat')
load('PARC500_MC.dat')
load('PARC500_MD.dat')
load('PARC500_SA.dat')
load('PARC500_GM.dat')
load('PARC500_CT.dat')
load('PARC500_FA.dat')

load('group.dat')
load('age.dat')
load('sex.dat')

nregs=308; % number of regions
nsubs=length(group); % number of subjects- 151 for Maastricht, 115 for Dublin and 146 for Cobre
```
If you are using these data for your own analyses, please note that there is one outlier subject in the Cobre MD and FA values (as discussed in the SI of the paper).

## Calculate the morphometric similarity matrices:

```
% z-score the inputs:
PARC500_CT_zscore=zscore(transpose(PARC500_CT));
PARC500_SA_zscore=zscore(transpose(PARC500_SA));
PARC500_GM_zscore=zscore(transpose(PARC500_GM));
PARC500_MC_zscore=zscore(transpose(PARC500_MC));
PARC500_GC_zscore=zscore(transpose(PARC500_GC));
PARC500_FA_zscore=zscore(transpose(PARC500_FA));
PARC500_MD_zscore=zscore(transpose(PARC500_MD));

% Create a cell for each subject with all of the required inputs:
clear subj_features7
for subj=1:nsubs
    subj_features7{1,subj}(:,1)=PARC500_CT_zscore(:,subj);
    subj_features7{1,subj}(:,2)=PARC500_SA_zscore(:,subj);
    subj_features7{1,subj}(:,3)=PARC500_GM_zscore(:,subj);
    subj_features7{1,subj}(:,4)=PARC500_MC_zscore(:,subj);
    subj_features7{1,subj}(:,5)=PARC500_GC_zscore(:,subj);
    subj_features7{1,subj}(:,6)=PARC500_FA_zscore(:,subj);
    subj_features7{1,subj}(:,7)=PARC500_MD_zscore(:,subj);
end

% Calculate the MS matrices by correlating all inputs and set the diagonal to zero:
for subj=1:nsubs
    subj_MSN_7{1,subj}=corr(transpose(subj_features7{1,subj}));
    subj_MSN_7{1,subj}(logical(eye(size(subj_MSN_7{1,subj})))) = 0;
end
```

## Global differences in morphometric similarity:

```
clear meanMS_regional
for subj=1:nsubs
    meanMS_regional(subj,:)=sum(subj_MSN_7{1,subj})./(nregs-1);
end

x1=age;
x2=sex;
X = [ones(size(x1)) x1 x2 x1.*x2];

% Calculate regional residuals:
clear myresid_region
for region=1:nregs
    y = meanMS_regional(:,region);
    [b,bint,resid]=regress(y,X);
    YFIT = b(1) + b(2)*x1 + b(3)*x2 + b(4)*x1.*x2;
    myresid_region(region,:)=y-YFIT;
end

cons=find(group==1);
pats=find(group==2);

x = reshape(myresid_region(:,cons),[nregs*length(cons),1]);
y = reshape(myresid_region(:,pats),[nregs*length(pats),1]);

% Plot histograms of regional residuals for control subjects and patients:
figure
h1 = histogram(x);
hold on
h2 = histogram(y);
h1.Normalization = 'probability';
h1.BinWidth = 0.01;
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
legend('Controls','Patients')

% Calculate mean MS:
for subj=1:nsubs
    meanMS(subj)=mean(meanMS_regional(subj,:));
end

tbl = table(age,sex,group,transpose(meanMS));
tbl.sex = categorical(tbl.sex);
tbl.group = categorical(tbl.group);
lm = fitlm(tbl,'Var4~age*sex+group');
p_mean=lm.Coefficients{4,4} % p-value for the effect of group on mean MS

% Plot box plot:

for subj=1:nsubs
    myresid_region_mean(subj)=mean(myresid_region(:,subj));
end

figure
boxplot(myresid_region_mean,group,'notch','on','Labels',{'Controls','Patients'})

```

## Regional differences in morphometric similarity:

The following code calculates a t-statistic for regional differences in morphometric similarity:

```
dummy=meanMS_regional;

clear mytstat mypval
for region=1:nregs
  tbl = table(age,sex,group,dummy(:,region));
  tbl.sex = categorical(tbl.sex);
  tbl.group = categorical(tbl.group);
  lm = fitlm(tbl,'Var4~age*sex+group');
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

The code below calculates the mean regional MS for all controls, from all three datasets.

First you need to extract meanMS_regional (as calculated above) for the control subjects only from all three datasets, e.g. for Maastricht:
```
meanMS_regional_Maast_con=meanMS_regional(find(group==1),:);
```

Then simply average them:

```
meanMS_con=mean(vertcat(meanMS_regional_Maast_con,meanMS_regional_Dublin_con,meanMS_regional_Cobre_con),1);
```

Plot the correlation between the regional mean control MS and the mean t-statistic:
```
figure
scatter(meanMS_con,meantstat,'x')
refline(0,0)
hold on
plot([0 0], ylim)
plot([0 0], ylim,'-b')
xlabel('Mean control MS')
ylabel('Mean t-statistic')

% Calculate the percentage of scatter points in each quadrant:

a=0;
b=0;
c=0;
d=0;

xvalues=meanMS_con;
yvalues=meantstat;

clear a b c d
for ind=1:nregs
    xval=xvalues(ind);
    yval=yvalues(ind);
    if ((xval<0)&&yval>0) % a is top left quadrant
        a=a+1;
        elseif ((xval>0)&&yval>0) % b is top right quadrant
        b=b+1;
        elseif ((xval<0)&&yval<0) % c is bottom left quadrant
        c=c+1;
        elseif ((xval>0)&&yval<0) % d is bottom right quadrant
        d=d+1;
    end
end

a=a/nregs
b=b/nregs
c=c/nregs
d=d/nregs
```

## Regional differences in MS- lh only

For the gene expression analyses, calculate the t-statistic for the lh only:

```
nregs_lh=152;
for subj=1:nsubs
    meanMS_regional_lh(subj,:)=sum(subj_MSN_7{1,subj}(1:nregs_lh,1:nregs_lh))./(nregs_lh-1);
end

dummy=meanMS_regional_lh;
clear mytstat mypval
for region=1:nregs_lh
  tbl = table(age,sex,group,dummy(:,region));
  tbl.sex = categorical(tbl.sex);
  tbl.group = categorical(tbl.group);
  lm = fitlm(tbl,'Var4~age*sex+group');
  mytstat(region)=lm.Coefficients{4,3};
  mypval(region)=lm.Coefficients{4,4};
end

dlmwrite('mytstat_Maast_lh.dat',mytstat)
dlmwrite('mypval_Maast_lh.dat',mypval)
```

## von Economo/Yeo networks:

The code below assesses whether there are differences in MS within particular von Economo classes or Yeo networks. You will need to import the lists of which Yeo network/von Economo class each of the 308 cortical regions belongs to, which can be found in the files 'Yeo_500_overlap.txt' and 'vonEcon_500_overlap.txt'. The Yeo network mapping was performed by Dr Jakob Seidlitz as part of the paper [Váša et al, Cereb Cortex. 2018](https://doi.org/10.1093/cercor/bhx249). The von Economo class mapping was performed by [Konrad Wagstyl](https://github.com/kwagstyl) and [Dr Kirstie Whitaker](https://github.com/kirstiejane) as part of the paper [Whitaker and Vértes, PNAS 2016](https://doi.org/10.1073/pnas.1601745113).

```
networks=Yeo500overlap; % or vonEcon500overlap

clear myt myp
for class=1:7
    myregions=find(networks==class);
    for subj=1:nsubs
        classreg(subj)=sum(meanMS_regional(subj,myregions));
    end

    tbl = table(age,sex,group,transpose(classreg));
    tbl.sex = categorical(tbl.sex);
    tbl.group = categorical(tbl.group);
    lm = fitlm(tbl,'Var4~age*sex+group');
    myt(class)=lm.Coefficients{4,3};
    myp(class)=lm.Coefficients{4,4};
end

class=4; % insert which class you're interested in here
myregions=find(networks==class);
for subj=1:nsubs
    classreg(subj)=sum(meanMS_regional(subj,myregions));
end

x1=age;
x2=sex;
X = [ones(size(x1)) x1 x2 x1.*x2];
y=transpose(classreg);
[b,bint,resid]=regress(y,X);
YFIT = b(1) + b(2)*x1 + b(3)*x2 + b(4)*x1.*x2;
myresid=y-YFIT;
[r p]=ttest2(myresid(cons),myresid(pats));

figure
boxplot(myresid,group,'notch','on','labels',{'Controls','Patients'})
```
