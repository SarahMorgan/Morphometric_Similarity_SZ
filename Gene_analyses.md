This file provides code to perform the gene analyses and the gene enrichment calculations.

## Gene expression PLS:

This code was written by Dr Petra Vértes and is taken from [Whitaker and Vértes, PNAS 2016](http://www.pnas.org/content/113/32/9105), please cite that paper if you this code in your own work.

```
load('gene_regional_expression_zscored.mat')

% Note that here we use the left hemisphere only
nregs=308;
nregs_lh=152;

X=gene_regional_expression(1:nregs_lh,:); % Predictors
Y=horzcat(mytstat_Maast_lh,mytstat_Dublin_lh,mytstat_Cobre_lh); % Response variable

% z-score:
X=zscore(X);
Y=zscore(Y);

%perform full PLS and plot variance in Y explained by top 15 components
%typically top 2 or 3 components will explain a large part of the variance
%(hopefully!)
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y);
dim=15;
plot(1:dim,cumsum(100*PCTVAR(2,1:dim)),'-o','LineWidth',1.5,'Color',[140/255,0,0]);
set(gca,'Fontsize',14)
xlabel('Number of PLS components','FontSize',14);
ylabel('Percent Variance Explained in Y','FontSize',14);
grid on

% plot samples in PLS space (using top 2 components)
%might want to colour dots by their score for learning style
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim); % no need to do this but it keeps outputs tidy

%%%% 2D plot of subjects in PLS component space %%%% 
figure
plot(XS(:,1), XS(:,2),'ok','MarkerSize',8,'MarkerFaceColor','r');
xlabel('PLS component 1','FontSize',14);
ylabel('PLS component 2','FontSize',14);
grid on

%%% plot correlation of PLS component 1 with t-statistic (from Cobre as an example):
figure
plot(XS(:,1),mytstat_Cobre_lh,'r.')
[R,p]=corrcoef(XS(:,1),mytstat_Cobre_lh) 
xlabel('XS scores for PLS component 1','FontSize',14);
ylabel('Cobre t-statistic- lh','FontSize',14);
grid on

% permutation testing to assess significance of PLS result as a function of
% the number of components (dim) included:

rep=1000;
for dim=1:8
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
temp=cumsum(100*PCTVAR(2,1:dim));
Rsquared = temp(dim);
    for j=1:rep
        %j
        order=randperm(size(Y,1));
        Yp=Y(order,:);

        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Yp,dim);

        temp=cumsum(100*PCTVAR(2,1:dim));
        Rsq(j) = temp(dim);
    end
dim
R(dim)=Rsquared
p(dim)=length(find(Rsq>=Rsquared))/rep
end
figure
plot(1:dim, p,'ok','MarkerSize',8,'MarkerFaceColor','r');
xlabel('Number of PLS components','FontSize',14);
ylabel('p-value','FontSize',14);
grid on
```

Bootstrap to get the gene list:

```
genes=genes20647; % this needs to be imported first
geneindex=1:20647;

%number of bootstrap iterations:
bootnum=1000;

% Do PLS in 2 dimensions (with 2 components):
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

%store regions' IDs and weights in descending order of weight for both components:
[R1,p1]=corr([XS(:,1),XS(:,2)],mytstat_Dublin_lh);

%align PLS components with desired direction for interpretability 
if R1(1,1)<0  %this is specific to the data shape we were using - will need ammending
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0 %this is specific to the data shape we were using - will need ammending
    stats.W(:,2)=-1*stats.W(:,2);
    XS(:,2)=-1*XS(:,2);
end

[PLS1w,x1] = sort(stats.W(:,1),'descend');
PLS1ids=genes(x1);
geneindex1=geneindex(x1);
[PLS2w,x2] = sort(stats.W(:,2),'descend');
PLS2ids=genes(x2);
geneindex2=geneindex(x2);

%print out results
csvwrite('PLS1_ROIscores.csv',XS(:,1));
csvwrite('PLS2_ROIscores.csv',XS(:,2));

%define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
PLS2weights=[];

%start bootstrap
for i=1:bootnum
    i
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled subjects
    Yr=Y(myresample,:); % define X for resampled subjects
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
      
    temp=stats.W(:,1);%extract PLS1 weights
    newW=temp(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights=[PLS1weights,newW];%store (ordered) weights from this bootstrap run
    
    temp=stats.W(:,2);%extract PLS2 weights
    newW=temp(x2); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS2weights=[PLS2weights,newW]; %store (ordered) weights from this bootstrap run    
end

%get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');
PLS2sw=std(PLS2weights');

%get bootstrap weights
temp1=PLS1w./PLS1sw';
temp2=PLS2w./PLS2sw';

%order bootstrap weights (Z) and names of regions
[Z1 ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);
[Z2 ind2]=sort(temp2,'descend');
PLS2=PLS2ids(ind2);
geneindex2=geneindex2(ind2);

%print out results
% later use first column of these csv files for pasting into GOrilla (for
% bootstrapped ordered list of genes) 
fid1 = fopen('PLS1_geneWeights.csv','w')
for i=1:length(genes)
  fprintf(fid1,'%s, %d, %f\n', PLS1{i}, geneindex1(i), Z1(i));
end
fclose(fid1)

fid2 = fopen('PLS2_geneWeights.csv','w')
for i=1:length(genes)
  fprintf(fid2,'%s, %d, %f\n', PLS2{i},geneindex2(i), Z2(i));
end
fclose(fid2)
```

## Gene enrichments:

This code calculates enrichments of the known gene lists (used for Gandal, DISEASES and GAD). The input gene list must be a list of entrez IDs.

```
knowngenes=gad_entrez; % List of known genes (entrez IDs)
names=pls1_entrez; %flip(PLS1);

gene_pos=zeros(1,length(knowngenes));

for index=1:length(knowngenes)
    result=find(names==knowngenes(index));
    if (isempty(result)==0) %
        gene_pos(index)=result;
    else
        gene_pos(index)=nan;
    end
end

gene_pos_median=nanmedian(gene_pos)

% Check whether enriched:
lengthnonan=size(find(isnan(PGC_pos)==0),2);
bground=probeInformation_maxi.EntrezID;
bground_red=intersect(bground,pls1_entrez);

randno=10000; % no. of randomisations to perform

% Randomise:
clear myRpos myRposmedian
for rand=1:randno
    randgenes=bground_red(randperm(length(bground_red), lengthnonan)); % pick lengthnonan genes from background list at random
    clear positionR
    for index=1:lengthnonan
        positionR(index)=find(names==randgenes(index));
    end
myRposmedian(rand)=median(positionR);
end

myRposmedian_sort=sort(transpose(myRposmedian),'descend');
p_pos_median=length(find(myRposmedian_sort<=gene_pos_median))/randno
```
