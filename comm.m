function solution = comm(model,exp,txt1,condition,threashold)
% Use this algorithm to integrate a co-expression network and genome scale 
% metabolic model. This algorithm calculate the reaction flux distribution
% in each condition by applying quadratic programming.
%
% USAGE:
%
%    solution = comm(model,exp,txt1,threashold,condition)
%
% INPUTS:
%
%    model:               input model (COBRA model structure)
%    exp:                 expression profile corresponding gene names that
%                         extract from gene expression profile file
%    txt1:                text that extract from gene expression profile
%                         file
%
% OPTIONAL INPUT:
%    threshold:           value of correlation coefficient for construct 
%                         co-expression network (default value 0.9)
%    condition:           row vector of number of condition corresponding 
%                         condition in exp
%                         (default value - 1:size(exp,2))
%
% OUTPUTS:
%    solution:             flux distribution table corresponding reaction
%                          flux names
%
% 
if nargin< 4 || isempty(condition)
      condition =1:size(exp,2);
end
if nargin< 5 ||isempty(threashold)
      threashold=0.9;
end

% construct the template model
model_n=model;
model_n.lb(model_n.lb>=0)=0;
model_n.lb(model_n.lb<0)=-1000;
model_n.ub(model_n.ub<=0)=0;
model_n.ub(model_n.ub>0)=1000;

%convert to irreversible format
[modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(model_n);

%check gene name in gene expression profile agree gene name in metabolic
%model
geneinMet=zeros(size(txt1,1)-1,1);
for i=2:size(txt1,1)
    for j=1:size(modelIrrev.genes,1)
        if string(txt1(i))==string(modelIrrev.genes(j))
            geneinMet(i-1)=1;
        end
    end
end
if sum(sum(geneinMet))==0
    disp("gene name in gene expression data is not agree gene name in metabolic model")
end

%Choose only gene expression profile agree gene in metabolic model
exp1=exp(geneinMet==1,:);
txt2=txt1([0;geneinMet]==1,1);

%remove missing data
[exp1,TF]=rmmissing(exp1,1);
txt2=txt2(~TF,1);

%Construct co-expression network 
cor=corr((exp1)');
cor1=cor>=threashold;

%find gene pair that have hihg correlation in co-expression network
ind_cor_gene=zeros(1,2);
l=1;
for i =1:size(cor,1)
    for j =i:size(cor,1)
        if (cor1(i,j)~=0)&&(i~=j)
            ind_cor_gene(l,1)=i;
            ind_cor_gene(l,2)=j;
            l=l+1;
         end
    end
end
co_gene=cell(size(ind_cor_gene,1),2);
co_gene(:,1)=txt2(ind_cor_gene(:,1),1);
co_gene(:,2)=txt2(ind_cor_gene(:,2),1);

%find reaction that satisfies gene 
Name_cor_rxn={};
for i=1:size(co_gene,1)
    [z1,Name_cor_rxn{i,1}]=findRxnsFromGenes(modelIrrev, co_gene(i,1),1,1);
    [z2,Name_cor_rxn{i,2}]=findRxnsFromGenes(modelIrrev, co_gene(i,2),1,1);
end
f=0;
h=0;
ind_cor_rxn=zeros(10,2);
for i=1:size(Name_cor_rxn,1)
    if ~isempty(Name_cor_rxn{i,1}) && ~isempty(Name_cor_rxn{i,2})
    for j=1:size(Name_cor_rxn{i,1},1)
        f=f+size(Name_cor_rxn{i,2},1);
        ind_cor_rxn(h+1:h+size(Name_cor_rxn{i,2},1),1)=findRxnIDs(modelIrrev,Name_cor_rxn{i,1}{j,1}); 
            for w=1:size(Name_cor_rxn{i,2},1)
                ind_cor_rxn(h+w,2)=findRxnIDs(modelIrrev,Name_cor_rxn{i,2}{w,1}); 
            end
            h=f;
    end
    end
end

%construct reaction pair matrix
R=zeros(length(modelIrrev.rxns),length(modelIrrev.rxns));
for i=1:size(ind_cor_rxn,1)
    if ind_cor_rxn(i,1)~=ind_cor_rxn(i,2)
        R(ind_cor_rxn(i,1),ind_cor_rxn(i,2))= 1;       
    end
end

%construct irreversible reaction matrix that are decomposed from  the same reversible reaction 
Re=zeros(size(R));
for i=1:length(model.rxns)
    if length(rev2irrev{i,1})==2
        Re(i,rev2irrev{i,1}(2))=1/2;
        Re(rev2irrev{i,1}(2),i)=1/2;
    end
end

%process gene expression data based on GPR association
gene_exdat=zeros(length(modelIrrev.genes),1);
nupb=zeros(length(modelIrrev.rules),size(exp,2));
for ch=1:size(exp,2)
for i=1:length(modelIrrev.genes)
    cc=0;
    for j=1:length(txt1(:,1))-1
        if string(modelIrrev.genes(i))==string(txt1(j+1,1))
            gene_exdat(i)=exp(j,ch);
            cc=cc+1;
        end
    end
    if gene_exdat(i)==0 && cc==0
        gene_exdat(i)=1000;
    end
end

for i = 1:length(modelIrrev.rules)
    rulegene=modelIrrev.rules{i};
    if isempty(rulegene)
        rsum=1000;
    else
        newrule = split(rulegene,"|");
        nrule=size(newrule,1);
        rsum=0;
        for j=1:nrule
            newrule1=newrule{j};
            nnrule=length(newrule1);
            rmin=inf;
            for k=1:nnrule
                if newrule1(k)=='x'
                    r1=k+2;
                end
                if newrule1(k)==')' && newrule1(k-1)~=' ' && newrule1(k-1)~=')'
                    rmin=min(rmin,gene_exdat(str2num(convertCharsToStrings((newrule1(r1:k-1))))));
                end    
            end
            rsum=rsum+rmin;
        end 
    end
    nupb(i,ch)=rsum;
end

end
pos_nupb=nupb>=1000;
nupb(pos_nupb)=max(nupb(~pos_nupb));
nupb=filloutliers(nupb,'clip','mean',1);

%construct table for reporting the result
solution=table(model.rxns);
 n=0;
 
%construct model and calculate flux distribution
for ch=condition
    n=n+1;
disp("Condition :"); disp(ch);
model3=changeRxnBounds(modelIrrev,modelIrrev.rxns, nupb(:,ch),'u');
solution1 = optimizeCbModel(model3);

Trans0=zeros(size(modelIrrev.mets,1),size(modelIrrev.rxns,1));
Trans2=-1*eye(size(modelIrrev.rxns,1));
S2=zeros(size(modelIrrev.rxns,1));
for i=1:length(modelIrrev.rxns)
    S2(i,i)=1/max(nupb(:,ch));
end

Obj4=[modelIrrev.c' zeros(1,size(modelIrrev.rxns,1))] ;

lob=[model3.lb; (-1)*inf*ones(size(modelIrrev.rxns,1),1)];
upb=[model3.ub; inf*ones(size(modelIrrev.rxns,1),1)];

O=[zeros(size(R)) zeros(size(R)) ; zeros(size(R)) R]; 
Aeq=[modelIrrev.S Trans0 ; S2 Trans2 ; Obj4]; 
beq=[zeros(size(modelIrrev.mets,1),1); (-1)*ones(size(modelIrrev.rxns,1),1) ; 0.99*solution1.f];

model2=struct;
model2.lb=lob;
model2.ub=upb; 
model2.A=sparse(Aeq);
model2.sense = [char('='* ones(size(model2.A,1)-1, 1)) ; char('>')]; 
model2.rhs=beq;
model2.modelsense = 'max'; 
numrxn=[1:length(modelIrrev.rxns)]; 
j=1;
for i=1:length(model.rxns)
    if length(Re(Re(i,:)>0))==1
model2.quadcon(j).Qrow = i;
model2.quadcon(j).Qcol = numrxn(Re(i,:)>0);
model2.quadcon(j).Qval = 1.0;
model2.quadcon(j).rhs = 0.0;
model2.quadcon(j).q = sparse(zeros(2*size(Re,1),1));
model2.quadcon(j).sense = '=';
j=j+1;
    end
end
model2.Q = sparse(O);
params.NonConvex = 2;
result = gurobi(model2, params);
x=result.x;
sol_flux=[];
for i=1:length(model.rxns)
    if length(rev2irrev{i,1})==1
        sol_flux(i,1)=x(i,1);
    else
        sol_flux(i,1)=(x(rev2irrev{i,1}(1,1))-x(rev2irrev{i,1}(1,2)));
    end
end
solution(:,n+1)=table(sol_flux);
end
filename = 'result.csv';
writetable(solution,filename)
end
