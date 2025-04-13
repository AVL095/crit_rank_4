from scipy import stats
from more_itertools import distinct_permutations as dp
import pandas as pd
import numpy as np
from pandas.core.common import flatten
from time import time

############################################################################
def Mood(a,k,n,m):
    msmal=2
    if (m[0] < m[1]): msmal = 1
    r= 0
    for i in range(n):
        if (a[i]==msmal): r=r+(i+1-(n+1)/2)**2
    return(r)
###################################################
def Wilcoxon(a,k,n,m):
    msmal=2
    if (m[0] < m[1]): msmal = 1
    r= 0
    for i in range(n):
        if (a[i]==msmal): r=r+i+1
    return(r)
###################################################
def Kruskal(a,k,n,m):
    s=0
    for j in range(k):
        r=0
        for i in range(n):
            if(a[i]==j+1): r=r+i+1
        s=s+r*r/m[j]
    return(12.*s/(1.*n*(n+1.))-3.*(1.*n+1.))
####################################################
def Leman(a,k,n,m):
    r1 = 0; r2 = 0; k1 = 0; k2 = 0
    for i in range(n):
        if(a[i]==1):
            r1=r1+(i+1-k1)**2
            k1=k1+1
        if(a[i]==2):
            r2=r2+(i+1-k2)**2
            k2=k2+1
    zleman=(r1*m[0]+r2*m[1]+m[0]*m[1]/6)/(m[0]*m[1])**2-2/3
    return(zleman)
#######################################################################
def David(a,k,n,m):
    msmal=2
    if(m[0]<m[1]): msmal = 1
    r=0
    for i in range(n):
        if a[i]==msmal:
            if i<=n/2: r+=n/2-i
            if i>n/2: r+=i-n/2-1
    return(r)
########################################################################
def Ansari(a,k,n,m):
    r=0
    for i in range(n):
        if (a[i]==1): r=r+((n+1)/2-abs(i+1-(n+1)/2))
    return(r)
 ################################Чтение данных###################################################       
def rank_exact(crit,m):
    k=len(m)
    n=sum(m)
    dc={'Kruskal':Kruskal,'Wilcoxon':Wilcoxon,'Mood':Mood,'Ansari':Ansari,'Leman':Leman,'David':David}
    a=list(flatten([[i+1]*m[i] for i in range(k)]))
    now=time()
    h=pd.Series(map(lambda x:dc[crit](x,k,n,m),dp(a)))
    kx=len(h)
    x=h.value_counts().sort_index()
    ksize=len(x)
    w=x.index
    w=np.array(w)
    pw=x.values/sum(x.values)
    pw=np.array(pw)
    s1=0
    for i in range(ksize):
        s1+=pw[i]; pw[i]=s1
    stime=time()-now
    return(kx,ksize,stime,w,pw)


################Точное распределение критерия Уилкоксона###########################################  

def Wilcoxon_Analitic(m):
    work=[]
    w=[]
    pw=[]
    wrange=[]
    minmn=min(m)
    maxmn=max(m)
    mn=m[0]*m[1]+1
    n1=maxmn+1
    work=[0]*(mn+2)
    w = [0]*(mn + 2)
    for i in range(1,n1+1): w[i] = 1
    for i in range(n1+1,mn+1): w[i] = 0
    pw=[0]*(mn)
    wrange=[0]*(mn)
    work[1]=0
    inx=maxmn
    for i in range(2,minmn+1):
        work[i] = 0
        inx+=maxmn
        n1=inx+2
        kk=1+inx//2
        k=i
        for j in range(1, kk+1):
            k+=1
            n1-=1
            sum=w[j]+work[j]
            w[j]=sum
            work[k]=sum-w[n1]
            w[n1]=sum

    for i in range(mn):
        pw[i]=float(w[i+1])
    sum=0
    for i in range(mn):
        sum+=pw[i]
        wrange[i]=float(i+minmn*(minmn+1)/2)
        pw[i]=sum
    for i in range(mn):
        pw[i]=pw[i]/sum
    stime=0 
    return(mn,mn,stime,wrange,pw)

#########################################################################

def Ansari_Analitic(m):
    x1=stats.norm.rvs(loc=0,scale=1,size=m[0])
    x2=stats.norm.rvs(loc=0,scale=1,size=m[1])
    AB,pval,astart,a1=stats.ansari(x1,x2)
    hnum=len(a1)
    ksize=hnum
    total=a1.sum()
    pw=[astart+i for i in range(hnum)]
    w=[sum(a1[0:i+1]/total) for i in range(hnum)]
    stime=0
    return(hnum,ksize,stime,pw,w)

###########################################################

crit="Leman"

txt="Inp1/"+crit+".inp"
finp=open(txt)
st=finp.readline()
m=tuple(map(int,finp.readline().split()))
finp.close()

if crit=="Ansari_Analitic":
    [hnum,ksize,stime,pw,w]=Ansari_Analitic(m)
elif crit=="Wilcoxon_Analitic" :
    [hnum,ksize,stime,pw,w]=Wilcoxon_Analitic(m)
else:
    [hnum,ksize,stime,pw,w]=rank_exact(crit,m)

txt="Out/"+crit+".out"
fout=open(txt,'w')
print("Time: ",stime,file=fout) 
print("Criterion: ",crit,file=fout) 
print("Samples: ",m,file=fout)
print("Sum: ",sum(m),file=fout)
print("Variants: ",hnum,file=fout)
print("Size: ",ksize,file=fout)
print("Statistics",pw,file=fout)
print("P-value",w,file=fout)
fout.close()
