let s2pi = 2.50662827463100050242;

function invnormaldistribution(y0)  {

let expm2 = 0.13533528323661269189;
let maxrealnumber=Math.pow(10,300);
let minrealnumber=Math.pow(10,-300);

            let result = 0;
            let x = 0;
            let y = 0;
            let z = 0;
            let y2 = 0;
            let x0 = 0;
            let x1 = 0;
            let code = 0;
            let p0 = 0;
            let q0 = 0;
            let p1 = 0;
            let q1 = 0;
            let p2 = 0;
            let q2 = 0;
            let zz;
            
            if(y0<=0) {
                result = -maxrealnumber;
                return result;
            }
            if(y0>=1) {
                result = maxrealnumber;
                return result;
            }
            code = 1;
            y = y0;
            if(y>1.0-expm2) {
                y = 1.0-y;
                code = 0;
                }
            if(y>expm2) {
                y = y-0.5;
                y2 = y*y;
                p0 = -59.9633501014107895267;
                p0 = 98.0010754185999661536+y2*p0;
                p0 = -56.6762857469070293439+y2*p0;
                p0 = 13.9312609387279679503+y2*p0;
                p0 = -1.23916583867381258016+y2*p0;
                q0 = 1;
                q0 = 1.95448858338141759834+y2*q0;
                q0 = 4.67627912898881538453+y2*q0;
                q0 = 86.3602421390890590575+y2*q0;
                q0 = -225.462687854119370527+y2*q0;
                q0 = 200.260212380060660359+y2*q0;
                q0 = -82.0372256168333339912+y2*q0;
                q0 = 15.9056225126211695515+y2*q0;
                q0 = -1.18331621121330003142+y2*q0;
                x = y+y*y2*p0/q0;
                x = x*s2pi;
                result = x;
                return result;
            }
            zz=Math.log(y);
            x = Math.sqrt(-(2.0*zz));
            x0 = x-Math.log(x)/x;
            z = 1.0/x;
            if(x<8.0) {
                p1 = 4.05544892305962419923;
                p1 = 31.5251094599893866154+z*p1;
                p1 = 57.1628192246421288162+z*p1;
                p1 = 44.0805073893200834700+z*p1;
                p1 = 14.6849561928858024014+z*p1;
                p1 = 2.18663306850790267539+z*p1;
                p1 = -(1.40256079171354495875*0.1)+z*p1;
                p1 = -(3.50424626827848203418*0.01)+z*p1;
                p1 = -(8.57456785154685413611*0.0001)+z*p1;
                q1 = 1;
                q1 = 15.7799883256466749731+z*q1;
                q1 = 45.3907635128879210584+z*q1;
                q1 = 41.3172038254672030440+z*q1;
                q1 = 15.0425385692907503408+z*q1;
                q1 = 2.50464946208309415979+z*q1;
                q1 = -(1.42182922854787788574*0.1)+z*q1;
                q1 = -(3.80806407691578277194*0.01)+z*q1;
                q1 = -(9.33259480895457427372*0.0001)+z*q1;
                x1 = z*p1/q1;
            }
            else  {
                p2 = 3.23774891776946035970;
                p2 = 6.91522889068984211695+z*p2;
                p2 = 3.93881025292474443415+z*p2;
                p2 = 1.33303460815807542389+z*p2;
                p2 = 2.01485389549179081538*0.1+z*p2;
                p2 = 1.23716634817820021358*0.01+z*p2;
                p2 = 3.01581553508235416007*0.0001+z*p2;
                p2 = 2.65806974686737550832*0.000001+z*p2;
                p2 = 6.23974539184983293730*0.000000001+z*p2;
                q2 = 1;
                q2 = 6.02427039364742014255+z*q2;
                q2 = 3.67983563856160859403+z*q2;
                q2 = 1.37702099489081330271+z*q2;
                q2 = 2.16236993594496635890*0.1+z*q2;
                q2 = 1.34204006088543189037*0.01+z*q2;
                q2 = 3.28014464682127739104*0.0001+z*q2;
                q2 = 2.89247864745380683936*0.000001+z*q2;
                q2 = 6.79019408009981274425*0.000000001+z*q2;
                x1 = z*p2/q2;
            }
            x = x0-x1;
            if( code!=0 ) x = -x;
            result = x;
            return result;
          }
//********Нормальные порядковые статистики**************************************

function ordern(n,pr,ps) {

/*
  Computes David-Johnson approximation for mean and covariance between rth
  and sth order statistics from the normal dist. for a sample size n.
  pr=r/(n+1);ps=s/(n+1);
  vorder[0] - mean;
  vorder[1] - covariance between rth and sth order statistics 
  if r=s then vorder[1] - variance
*/
  let qr,qs,xr,xs,dr,ds;
  let xr1,xr2,xr3,xr4,xr5,xr6;
  let xs1,xs2,xs3,xs4,xs5,xs6;
  let z1,z2,z3,z4,z5,z6,z7;
  let er,crs;
  let vorder=[];

 qr=1-pr;qs=1-ps; 
xr=invnormaldistribution(pr);
xs=invnormaldistribution(ps);
dr=s2pi*Math.exp(xr*xr/2.);
ds=s2pi*Math.exp(xs*xs/2.);
xr1=dr;
xr2=xr*Math.pow(dr,2);
xr3=(2*xr*xr+1)*Math.pow(dr,3);
xr4=(6*xr*xr*xr+7*xr)*Math.pow(dr,4);
xr5=(24*Math.pow(xr,4)+46*xr*xr+7)*Math.pow(dr,5);
xr6=(120*Math.pow(xr,5)+326*Math.pow(xr,3)+127*xr)*Math.pow(dr,6);
xs1=ds;
xs2=xs*Math.pow(ds,2);
xs3=(2*xs*xs+1)*Math.pow(ds,3);
xs4=(6*Math.pow(xs,3)+7*xs)*Math.pow(ds,4);
xs5=(24*Math.pow(xs,4)+46*Math.pow(xs,2)+7)*Math.pow(ds,5);
xs6=(120*Math.pow(xs,5)+326*Math.pow(xs,3)+127*xs)*Math.pow(ds,6);

er=xr+pr*qr*xr2/(2*(n+2))+pr*qr*((qr-pr)*xr3/3+pr*qr*xr4/8)/Math.pow(n+2,2)+pr*qr*(-(qr-pr)*xr3/3+(Math.pow(qr-pr,2)-pr*qr)*xr4/4+qr*pr*(qr-pr)*xr5/6+Math.pow(qr*pr,2)*xr6/48)/Math.pow(n+2,3);

z1=(qr-pr)*xr2*xs1+(qs-ps)*xr1*xs2+pr*qr*xr3*xs1/2+ps*qs*xr1*xs3/2+pr*qs*xr2*xs2/2;
z1=z1*pr*qs/Math.pow(n+2,2);
z2=-(qr-pr)*xr2*xs1-(qs-ps)*xr1*xs2+(Math.pow(qr-pr,2)-pr*qr)*xr3*xs1;
z3=(Math.pow(qs-ps,2)-ps*qs)*xr1*xs3+(1.5*(qr-pr)*(qs-ps)+0.5*ps*qr-2*pr*qs)*xr2*xs2;
z4=(5/6)*pr*qr*(qr-pr)*xr4*xs1+(5/6)*ps*qs*(qs-ps)*xr1*xs4+(pr*qs*(qr-pr)+0.5*pr*qr*(qs-ps))*xr3*xs2;
z5=(pr*qs*(qs-ps)+0.5*ps*qs*(qr-pr))*xr2*xs3 + (1 /8)*Math.pow(pr*qr,2)*xr5*xs1+(1/8)*Math.pow(ps*qs,2)*xr1*xs5;
z6=0.25*Math.pow(pr,2)*qr*qs*xr4*xs2+0.25*pr*ps*Math.pow(qs,2)*xr2*xs4+(2*Math.pow(pr*qs,2)+3*pr*qr*ps*qs)*xr3*xs3/12;
z7=z2+z3+z4+z5+z6;
crs=z1+pr*qs*z7/Math.pow(n+2,3)+pr*qs*xr1*xs1/(n+2);

vorder[0]=er;
vorder[1]=crs;
return(vorder);
}
//***********cnm=n!/m!*(n-m)!********************************
function cnm(n,m) {
  var s1,s2,i;
  s1=0; s2=0;
  for (i=m+1;i<=n;i++) s1=s1+Math.log(i);
    for (i=1;i<=n-m;i++)   s2=s2+Math.log(i);
return Math.exp(s1-s2);
}

