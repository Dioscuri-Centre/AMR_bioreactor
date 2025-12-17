// this program simulates a population of cells growing in a turbidostat, using the simplistic tau leaping algorithm
// cells can either replicate or mutate with probability mu
// antibiotic is injected at some time t_ant which causes cessation of growth of the wild type but a mutant can still reproduce
// simulation is stopped when either OD 0 or od_crit is reached
// the simulation is then repeated some number of times 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SQR(x) (x)*(x)
//#define SWAPD(x, y) tempd = (x); (x) = (y); (y) = tempd
//#define SWAP(x, y) temp = (x); (x) = (y); (y) = temp
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

typedef unsigned long long int QWORD ;

const float vol=25 ; // volume [ml]
float gama_res=1e-9 ; // mutation probability per replication   3e-11 gives 3/20 successes
//const QWORD capacity=3e9 ; // cfu/ml
const double n_to_OD=2e8 ; // cfu/ml which gives OD_linear=1
const QWORD n_ini=QWORD(n_to_OD*0.001*vol) ; // initial number of cells
const float t_ant=5 ; 
float tdil=0.5;
const float dil=0.75 ; // dilution factor
const float od_crit=0.1 ;
float tmax=30 ;
char NAME[]="example_bs_2.0_br_1.5" ;

float gr[]={2.,1.5} ; // estimated from turbidostat OD curves
//#define RUN_ONCE // if defined, run only once
#define SAVE_MANY_TRAJECTORIES  // if defined, runs samples but it saves all data to separate files
const int samples=10 ; //96 ;


void err(char *reason)
{
  cout <<reason<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif  
  exit(0) ;
}

void err(char *reason, int a)
{
  cout <<reason<<": "<<a<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif    
  exit(0) ;
}

void err(char *reason, char *a)
{
  cout <<reason<<": "<<a<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif    
  exit(0) ;
}

void err(char *reason, double a)
{
  cout <<reason<<": "<<a<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif    
  exit(0) ;
}

static long long unsigned int _x=0x000100010001LL, _mul=0x0005deece66dLL, _add=0xbLL ;
double _drand48(void)  // works only on compilers with long long int!
{
  _x=_mul*_x+_add ; _x&=0xffffffffffffLL ;
  return (_x/281474976710656.0) ;
}

void _srand48(int a) { _x=a ; }

double gauss()
{
 static int iset=0;
 static float gset;
 float gasdev,v1,v2,r,fac;
 if (iset==0)
 {
  do
  {
	v1=2*_drand48()-1.0 ; v2=2*_drand48()-1.0 ;
	r=v1*v1+v2*v2;
  } while (r>1) ;
  fac=sqrt(-2.0*log(r)/r) ;
  gset=v1*fac ; gasdev=v2*fac ; iset=1 ;
 }
 else
 {
  gasdev=gset ; iset=0 ;
 }
 return (gasdev) ;
}

const unsigned int poisson(const double lambda)
{
  int k=0;
  const double target=exp(-lambda);
  double p=_drand48();
  while (p>target)
  {
    p*=_drand48();
    k++;
  }
  return k;
}

QWORD binomial_2(double p, QWORD n)
{
//  int i,j=0;
//  if (n<0LL) err("n<0LL") ;
  if (p<=1e-14 || n<=0) return 0 ;
  if (p>=1) return n ; if (p*n>=100.) {  
    QWORD m=(QWORD)(p*n+sqrt(p*(1-p)*n)*gauss()) ;
    if (m<0) m=0 ; if (m>n) m=n ;
    return m ;
  } else {
    QWORD x=0LL ;
    double q=-log(1-p), sum=0, e ;
    if (q<=0) return n ; //err("q", q) ;
    do {
      e=-log(1-_drand48()) ; sum+=e/(n-x) ; x++ ;
      if (sum<0) err("sum<0, p=",p) ;
    } while (sum<=q) ;
    if (x-1>n) err("x-1-n= ",int(x-1-n)) ;
    return x-1 ;
  }
}



double tt, tt0 ;
int nsam=0;
QWORD run_once()
{
  const int L=2;
  int i,j,k,l;  
  tt=tt0=0 ;
  QWORD n[L] ;
  for (i=0;i<L;i++) n[i]=0LL ;
  QWORD ntot=0 ;

  const double dt=1./int(1<<7) ;
  n[0]=n_ini ; // initial number of cells
  float c=0 ;  // initial antibiotic concentration

  float od ;
    
#ifdef RUN_ONCE
  ofstream data("turbidostat_cipro.dat") ;
#endif
#ifdef SAVE_MANY_TRAJECTORIES
	char txt[256] ;
	sprintf(txt,"turbidostat_%s_n=%d.dat",NAME,nsam) ;
	ofstream data(txt) ;
#endif
    
  do {      // main loop 
    tt+=dt ; //cout <<genotypes.size()<<endl ;

    if (tt>tt0+tdil || od>od_crit) { // dilution
      tt0=tt ;
      for (j=0;j<L;j++) n[j]=binomial_2(dil,n[j]) ;
    }

    for (j=0;j<L;j++) {
			double d=0.01 ;
			double b=gr[j] ;
			if (tt>t_ant && j==0) b=0 ;

      QWORD dtot=binomial_2((b+d)*dt,n[j]) ;
      QWORD db=binomial_2(b/(b+d),dtot), dd=dtot-db ;
//      s-=1.*db/(vol*capacity) ; if (s<0) s=0 ;
      if (db>0) {         // replication
        int dbmut ; 
        if (j<L-1) dbmut=binomial_2(gama_res,db) ; else dbmut=0 ;
        if (dbmut==0) {
          n[j]+=db ; //genotypes[c->gen]->number+=db ; 
        } else {
//            cout <<"dm="<<dbmut<<"\t" ;
          if (dbmut>db) err("dbmut") ;
          n[j]+=db-dbmut ; //c->number+=db-dbmut ; 
//              cout<<j<<"->"<<c->gen<<","<<c->n<<" " ;
          n[j+1]+=dbmut ;
              //cout <<nc.gen<<":" ;
//              cout <<"x"<<c->gen<<endl ;
        }
      }
      if (dd>0) {         // death
        n[j]-=dd ; if (n[j]<0) err("n[j]<0");
      } 
    }
    ntot=0LL ;
    for (j=0;j<L;j++) ntot+=n[j] ;
    od=ntot/(n_to_OD*vol) ;
//    for (j=0;j<L;j++) cout <<n[j]<<"\t" ; cout <<endl ;
#if (defined(RUN_ONCE) || defined(SAVE_MANY_TRAJECTORIES))
    static int ttt=0 ;
    if (ttt%1==0) {
      data<<tt<<"\t" ; for (j=0;j<L;j++) data <<1.*n[j]<<"\t" ; data <<ntot<<"\t"<<od<<endl ;
    }
    ttt++ ;
#endif
  } while (ntot>0LL && tt<tmax) ;
#if (defined(RUN_ONCE) || defined(SAVE_MANY_TRAJECTORIES))
  data.close() ;
#endif
  return ntot ;
}


int main()
{
  _srand48(24) ;
#ifdef RUN_ONCE  
  run_once() ;
#elif defined(SAVE_MANY_TRAJECTORIES)
	for (nsam=0;nsam<samples;nsam++) run_once() ;
#else
  ofstream data("out.dat") ;
  for (gama_res=1e-10;gama_res<1e-8;gama_res*=pow(10,0.2)) {
    int nsurv=0 ;
    double tevol=0 ;
    for (int i=0;i<samples;i++) {
      QWORD n=run_once() ; 
  //    data <<n<<"\t"<<tt<<endl ;
      cout <<n<<" " ; 
      if (n>1e6) { nsurv++ ; tevol+=tt ; }
  //    data<<size<<"\t"<<nsurv<<"\t"<<samples<<endl ;
    }
    cout <<"gama="<<gama_res<<" nsurv="<<nsurv<<" "<<1.*tevol/nsurv<<endl ;
    data <<gama_res<<"\t"<<nsurv<<"\t"<<1.*tevol/nsurv<<endl ;
  }
  data.close() ;
  system("pause") ;  
#endif

}
