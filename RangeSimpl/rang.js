let num,znaki=7;
//******************************************************************
function main(crit) {
  let k,n,w2,z,knum,i,nn,pw=[],wrange=[],m=[];
 
  if(crit=="VanDerVarden" || crit=="Klotz" || crit=="Fisher" || crit=="Capon" || crit=="Ansari" || crit=="Mood" || crit=="Leman" || crit=="Wilcoxon"){
      k=2;m[0]=5;m[1]=7;
      nn=permrepeat(crit,k,m,wrange,pw);
   }
  if(crit=="Kruskal") {
      k=4;m[0]=2;m[1]=2;m[2]=2;m[3]=2;
      nn=permrepeat(crit,k,m,wrange,pw);
  }  
  if(crit=="Series") {
     k=2;m[0]=10;m[1]=12;
     nn=seriestest(m,wrange,pw);
   }
  if(crit=="Ansari-analitic") {
     k=2;m[0]=5;m[1]=7;
     nn=ansarianalitic(m,wrange,pw);
   }
   if(crit=="Wilcoxon-analitic" || crit=="Sigel") {
     k=2;m[0]=5;m[1]=7;
     nn=wilcoxonanalitic(m,wrange,pw);
   }
    critprint(crit,k,m,nn,wrange,pw);
}
//********************************************************
let SelectCrit={
    Mood: function(a,kx,n,m) {
    let i, msmal,r;
    msmal=2;
    if (m[0] < m[1]) msmal = 1;
    r= 0;
    for (i = 0; i <n; i++)  if (a[i]==msmal) r=r+Math.pow((i+1.- (n + 1.) / 2.),2);
    return r;
    },
//******************************************************************
    Leman: function(a,kx,n,m) {
    let i, k1,k2,zleman,r1,r2;
    r1 = 0.; r2 = 0.; k1 = 0; k2 = 0;
    for (i = 0; i < n; i++) {
        if (a[i] == 1) {
            r1=r1+(i-k1)**2;k1++;
        }
        if (a[i]==2) {
            r2+=(i-k2)**2;k2++;
        }
    }
    zleman=(r1*m[0]+r2*m[1]+m[0]*m[1]/6.)/(m[0]*m[1])**2-2./3.;
    return zleman;
    },
//*******************************************************
    Kruskal: function(a,kx,n,m) {
    let i,j;
    let s,r;
    s=0.;
     for (j=0;j<kx;j++) {
        r=0.;
        for (i=0;i<n;i++) if (a[i]==j+1) r=r+i+1;
        s=s+r*r/m[j];
     }
     return 12.*s/(1.0*n*(n+1.))-3.*(1.0*n+1.0);
    },

//***********************************************************
    Wilcoxon: function(a,kx,n,m) {
        r=0;
         for(i=0;i<n;i++) if(a[i]==1) r+=i+1;
        return r;
    },
//***********************************************************
    Ansari: function(a,kx,n,m) {
      r=0;
      for (i=0;i<n;i++) if (a[i]==1) r+=((n+1)/2-Math.abs(i+1-(n+1)/2));
      return r;
    },
//*************************************************************
    Fisher: function(a,kx,n,m) {
      r=0;
       for (i=0;i<n;i++) { 
         vorder=ordern(n,(i+1)/(n+1),(i+1)/(n+1));
         er=vorder[0];
         if (a[i]==1) r=r+er;
        }
         return r; 
     },
//**************************************************************
    VanDerVarden: function(a,kx,n,m) {
      r=0;
       for (i=0;i<n;i++) { 
         er=invnormaldistribution((i+1)/(m[0]+m[1]+1));
         if (a[i]==1) r=r+er;
        }
         return r; 
     },
//**************************************************************
    Capon: function(a,kx,n,m) {
      r=0;
       for (i=0;i<n;i++) { 
        vorder=ordern(n,(i+1)/(n+1),(i+1)/(n+1));
        er=vorder[0];
        cr=vorder[1];
        er=cr+er*er;
        if (a[i]==1) r=r+er;
       }
         return r; 
     },
//**************************************************************
 Klotz: function(a,kx,n,m) {
      r=0;
       for (i=0;i<n;i++) { 
          er=invnormaldistribution((i+1)/(m[0]+m[1]+1));
          if (a[i]==1) r=r+er*er;
        }
        
   return r;
  }
}
//*******************************************************************
      function permrepeat(crit,k,m,wrange,pw) {

     let n,i,j,l,r,a=[],h=[],nn,nnew,s,kk;

     n=0;
     for(i=0;i<k;i++) {
     for(j=0;j<m[i];j++) a[j+n]=i+1;
     n=n+m[i];
     }

     num=1;
     h[num]=SelectCrit[crit](a,k,n,m);
     
     while(true) {
        j=n-2;
        while (j != -1 && a[j] >= a[j + 1]) j--;
        if (j==-1) break;
        kk=n-1;
        while (a[j]>= a[kk]) kk--;
        swap(a,j,kk);
        l=j+1;r=n-1;
        while(l<r) swap(a,l++,r--);
        num++;
        h[num]=SelectCrit[crit](a,k,n,m);
      }

      h.sort(function(a,b){return a-b});
 
      nn=0;nnew=0;
      for (i = 0;i<num; i++) {
        if (parseFloat(h[i]).toFixed(5)==parseFloat(h[i+1]).toFixed(5)) {
          nnew++;
        }
        else {
            pw[nn]=nnew+1;
            wrange[nn]=h[i];
            nn++;
            nnew=0;
        }
    }
      s=0;
    for (i=0;i<nn;i++) {
        s += pw[i]/parseFloat(num);
        pw[i] = s;
    }
    return nn;
  } 
//*******************************************************************
   function swap(a,i,j)  {
        let s;
        s=a[i];a[i]=a[j];a[j]=s;
    }

//***********Точное распределение критерия Ансари-Брэдли*********
function ansarianalitic(m,astat,pw) {
  let min_val,max_val,i,ir,nrows,sum,a0,s;
  let w=[];

  min_val=m[0];max_val=m[1];
  nrows=gscale(min_val,max_val,w);
  //nrows=P.length;
  a0=Math.floor((min_val+1)/2)*(1+Math.floor(min_val/2));
   sum=0;
   for(i=0;i<nrows;i++) {
     astat[i]=a0+i;
     sum=sum+w[i+1];
    }
    
    s=0;
    for(i=0;i<nrows;i++) {
       s=s+w[i+1];
       pw[i]=s/sum;
    }
    return nrows;
}
//**********************************************************************
  function gscale(test,other,pw) {

/*
  algorithm as 93 appl. statist. (1976) vol.25, no.1
  from the sizes of two samples the distribution of the
  ansari-bradley test for scale is generated in array pw.
*/
	let ai,one,i,symm,m,lres,mm1,nm1,nm2,ier,mnow,ks,j,z,ndo,l1out;
        let n,ln1,ln2,nc,l2out,n2b1,n2b2,ln3,kk;
        let fadd=[],fimply=[],a2=[],a3=[];

//pw.length=1+Math.floor((m*n)/2);  

          one = 1;
          m = Math.min(test, other);
          if(m < 0) return(0);
          n=Math.max(test,other);
          lres = 1 + Math.floor((m * n)/2);
          symm=false;
          z=(m+n)%2;
          if(z==0) symm=true;
          mm1 = m - 1;
//*****************************************************         
       if(m<=2) {
          if(mm1<0) {
            pw[1]=1;return(lres);
          }
          
          if(mm1==0) ln1=start1(n, pw);
          if(mm1>0)  ln1=start2(n, pw);
            if(symm || (other > test)) return(lres);
            j = lres;
            ndo =Math.floor(lres/2);
          for(i=1;i<=ndo;i++) {
              ai = pw[i];
              pw[i]=pw[j];
              pw[j]=ai;
              j = j - 1;
          }
            return(lres);
       }
//***********************************************************
          nm1=n-1;nm2=n-2;mnow=3;nc=3;ier=0;
//**************************************************************
   while(true) {
    if(ier==0) {
        if((n % 2)!=1) {
          n2b1 = 3;
          n2b2 = 2;
          ln1=start2(n,pw);
          ln3=start2(nm2,a3);
          ln2=start1(nm1,a2);
//***********************************************************************          
          fadd=frqadd(a2,a3,ln2,l2out,ln3,n2b2);
          l2out=fadd[1];n2b2=fadd[2];
          ln2 = ln2 + nm1;
          fimply=imply(a2,a3,l2out,ln2,j,nc);
          ln2=fimply[1];j=fimply[2];nc=fimply[3];
//***************************************************************************
          nc = nc + 1;
          if(mnow==m) break;
          mnow = mnow + 1;
         } else {
          n2b1 = 2;
          n2b2 = 3;
          ln1=start1(n,pw);
          ln2=start2(nm1,a2);
        }
     }
//*******************************************************************
          fadd=frqadd(pw,a2,ln1, l1out,ln2, n2b1);
          l1out=fadd[1];n2b1=fadd[2];
          ln1 = ln1 + n;
          fimply=imply(pw,a3,l1out, ln1,ln3, nc);
          ln1=fimply[1];ln3=fimply[2];nc=fimply[3];
          nc = nc + 1;
          if(mnow==m) break;
          mnow = mnow + 1;
          fadd=frqadd(a2,a3, ln2, l2out,ln3, n2b2);
          l2out=fadd[1];n2b2=fadd[2];
          ln2 = ln2 + nm1;
          fimply=imply(a2,a3,l2out, ln2,j, nc);
          ln2=fimply[1];j=fimply[2];nc=fimply[3];
//***************************************************************************
          nc = nc + 1;
          if(mnow==m) break;
          mnow = mnow + 1;
          ier = 1;
  }
//********************************************************************
          if(symm) return(lres);
          ks = Math.floor((m + 3)/ 2);
          j = 1;
          for(i = ks;i<=lres;i++) {
           if(i>ln1) {
            pw[i]=a2[j];
           } else {
            pw[i]=pw[i]+a2[j];
          }
             j = j + 1;
          }
          if(other < test) return(lres);
          j = lres;
          ndo = Math.floor(lres/2);
          for(i=1;i<=ndo;i++) {
              ai = pw[i];
              pw[i]=pw[j];
              pw[j]=ai;
              j = j - 1;
          }
 return(lres);
	}
//********************************************************************
function start1(n,f) {
/*
	  algorithm as 93.1 appl. statist. (1976) vol.25, no.1
	  generates a 1,n ansari-bradley distribution in f.
*/
	let i,lout;
	lout=Math.floor(1+n/2);
	for(i=1;i<=lout;i++) f[i]=2;
	if ((n % 2)==0) f[lout]=1;
        return(lout);
  }
//************************************************************
  function start2(n,f) {

/*
	  algorithm as 93.2 appl. statist. (1976) vol.25, no.1
	  generates a 2,n ansari-bradley distribution in f.
*/
	let one,two,three,four,i,j,a,b,lt1,ndo,nu,lout;

	one=1;two=2;three=3;four=4;
	nu=n-n % 2;
	j=nu+1;lout=j;lt1=lout+1;
        ndo=Math.floor(lt1/2);
        a=one;b=three;
      for (i=1;i<=ndo;i++) {
	f[i]=a;f[j]=a;j=j-1;a=a+b;b=four-b;
      }
	if(nu==n) return(lout);
	nu=ndo+1;
	for(i=nu;i<=lout;i++) f[i]=f[i]+two;
	f[lt1]=two;lout=lt1;
	return(lout);
  }
//************************************************************
  function frqadd(f1,f2,l1in,l1out,l2,nstart) {
         

/*
	  algorithm as 93.3 appl. statist. (1976) vol.25, no.1
	  array f1 has twice the contents of array f2 added into it
	  starting with elements nstart and 1 in f1 and f2 respectively.
*/

	let i1,i2,nxt,fadd=[];

        i2=1;
          for(i1=nstart;i1<=l1in;i1++) {
              f1[i1]=f1[i1]+2*f2[i2];
              i2=i2+1;
          }
          nxt=l1in+1;
          l1out=l2+nstart-1;
          for(i1=nxt;i1<=l1out;i1++) {
              f1[i1]=2*f2[i2];
              i2=i2+1;
          }
          nstart=nstart+1;
          fadd[1]=l1out;fadd[2]=nstart;
          return(fadd);
 }
//*******************************************************************
  function imply(f1,f2,l1in,l1out,l2,noff) {
/*
	  algorithm as 93.4 appl. statist. (1976) vol.25, no.1
	  given l1in elements of an array f1, a symmetrical
	  array f2 is derived and added onto f1, leaving the
	  first noff elements of f1 unchanged and giving a
	  symmetrical result of l1out elements in f1.
*/
	let sum,diff,i2,i1,j2,j1,j2min,ndo,fimply=[];

    i2=1-noff;j1=l1out;j2=l1out-noff;l2=j2;
    j2min=Math.floor((j2 + 1)/2);
    ndo=Math.floor((l1out+1)/2);

          for(i1=1;i1<=ndo;i1++) {
              if(i2>0) {
                  sum=f1[i1]+f2[i2];
                  f1[i1]=sum;
              }
              else {
                  sum=f1[i1];
              }
              i2=i2+1;
              if(j2>=j2min) {
                  if(j1<=l1in) {
                      diff=sum-f1[j1];
                   }
                  else {
                     diff=sum;
                   }

              f2[i1]=diff;f2[j2]=diff;j2=j2-1;
             }
             f1[j1]=sum;j1=j1-1;
          }
            fimply[1]=l1out;fimply[2]=l2;fimply[3]=noff;
  return(fimply);
}

//*************************************
function wilcoxonanalitic(m,range,pw) {

/*
  AS 62 generates the frequencies for the Mann-Whitney U-statistic.
  Users are much more likely to need the distribution function.
  Code to return the distribution function has been added at the end
  of AS 62 by Alan Miller
*/


 let work=[],minmn,maxmn,inx,i,w=[];
 let n1,l,k,j,sum,mn;
  
  mn=m[0]*m[1]+1;maxmn=m[1];minmn=m[0];n1=maxmn+1;
 
     for(i=1;i<=n1;i++) w[i]=1;
      n1=n1+1;
      for(i=n1;i<=mn;i++) w[i]=0;
  
 work[1]=0;inx=maxmn;

    for(i=2;i<=minmn;i++) {
        work[i]=0;inx=inx+maxmn;n1=inx+2;l=1+inx/2;k=i;

       for(j=1;j<=l;j++) {
          k=k+1;n1=n1-1;
          sum=w[j]+work[j];
          w[j]=sum;work[k]=sum-w[n1];w[n1]=sum;
       }
    }
  
  //w[i] - Frequencies
     for(i=0;i<mn;i++) pw[i]=w[i+1];
    
    sum=0;
     for(i=0;i<mn;i++) {
        sum=sum+pw[i];
        range[i]=i+minmn*(minmn+1)/2; //Wilcoxon
        //range[i]=i; //Mann-Whitney
        pw[i]=sum;
     }
    for(i=0;i<mn;i++) pw[i]=pw[i]/sum;

    return(mn);
 }
//************ Критерий серий точное распределение********
function seriestest(m,wrange,p) {
   let i,kcur,s,m1,m2,n,k,pc,w2;

   m1=m[0];m2=m[1];n=m1+m2;
   w2=cnm(n,m1);
   s=0;
  for (i=0;i<=5000;i++) {
    kcur=2+i;pc=0;
   if (kcur%2==0) {
    pc=2*cnm(m1-1,kcur/2-1)*cnm(m2-1,kcur/2-1);
   } else {
    pc=cnm(m1-1,(kcur-1)/2)*cnm(m2-1,(kcur-1)/2-1)+cnm(m1-1,(kcur-1)/2-1)*cnm(m2-1,(kcur-1)/2);
   }
   s+=pc;
   //w[i]=parseFloat(pc);
   p[i]=parseFloat(s/w2);
   wrange[i]=parseFloat(kcur);
   if(p[i]>=1) break;
   }
   return i;
}
//*********************************************************
function critprint(crit,k,m,nn,wrange,pw) {
  let i,c;

  c=window.open("",crit,"toolbar=yes,menubar=yes,scroolbar=yes,width=650,height=400,left=400, top=150");
  c.document.write(crit);
  c.document.write("<br>");
  c.document.write("Size="+nn);
  c.document.write("<br>");
  for (i=0;i<k;i++) c.document.write(m[i]+"  ");
    c.document.write("<br>");
    for (i=0;i<nn;i++) {
     c.document.write((i+1)+": "+wrange[i].toFixed(znaki)+"  "+pw[i].toFixed(znaki));
     c.document.write("<br>");  
    }
 }
//***********************************************************

