This file provides code to perform the gene analyses, essentially Figure 3 of the paper and the gene enrichment calculations.

## Gene expression PLS:

This code was written by Dr Petra Vértes and introduced in [Whitaker and Vértes, PNAS 2016](http://www.pnas.org/content/113/32/9105), please cite that paper if you this approach in your own work.

```
X=gene_regional_expression; % Predictors
Y=meanMaastconpatcontrol; % Response variable

% z-score:

X=zscore(X); % Not strictly necessary (doesn't change results) because it's already been done
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

```

## Gene enrichments:

