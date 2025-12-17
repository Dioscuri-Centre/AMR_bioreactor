// this program simulates a population of cells growing in a turbidostat, using the simplistic tau leaping algorithm
// it simulates a standing variation experiment
// cells can only replicate
// antibiotic is injected at some time t_ant which causes cessation of growth of the wild type but a mutant can still reproduce
// simulation is stopped when either OD 0 or od_crit is reached
// the simulation is then repeated some number of times 
// OD returned takes into account a correction for RIF effect on OD to make it consistent with how the experimental data are processed
//		(mind that this is different to the long-term RIF simulation code)

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

int samples=1000;
float vol=20 ; // volume [ml]
//const QWORD capacity=3e9 ; // cfu/ml
const double n_to_OD=2e8 ; // cfu/ml which would give OD_linear=1 in the turbidostat had we grown bacteria beyond OD=0.1
QWORD n_ini ; // initial number of cells
float t_ant ; 
const float tdil=0.5;
float dil=0.75 ; // dilution factor - will be overriden
const float od_crit=0.1 ;
float tmax ;
char *NAME ;
float birthS,birthRbefore,birthRafter ; // estimated from turbidostat OD curves

float gr[2][2] ; // growth rate, [before=0/after=1][type]  
float t_regrowth, od_regrowth;


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

vector <float> ods,tts;
vector <QWORD> ns,nr ;
float RIF_lambda = 0.65, RIF_rif0 = 0.018 ; // for modelling RIF, will be overridden


double tt, tt0 ;
int nsam=0;
QWORD run_once(QWORD nini, float res_frac)
{
	ods.clear() ;
	ns.clear() ; 
	nr.clear() ;
	tts.clear() ;
  const int L=2;
  int i,j,k,l;  
  tt=tt0=0 ;
  QWORD n[L] ;
  for (i=0;i<L;i++) n[i]=0LL ;
  QWORD ntot=0 ;
  float od=0 ;

  n[0]=nini ; // initial number of cells
  n[1]=res_frac*n_ini ; // initial number of mutants
	// growth rates without the antibiotic
	gr[0][0]=birthS ; //*(0.98+_drand48()*0.01) ;
	gr[0][1]=birthRbefore ; //*(0.975+_drand48()*0.05) ; 	
	// growth rates in the antibiotic
	gr[1][0]=0 ; 
	gr[1][1]=birthRafter ; //*(0.975+_drand48()*0.05) ;   // error of the R growth rate
  const double dt=1./int(1<<8) ;

	int c=1 ; // for modelling RIF
	float rif=0,dec=0 ; // for modelling RIF

        
  do {      // main loop 
    tt+=dt ; //cout <<genotypes.size()<<endl ;

		//if (tt>t_ant) od_real+=0.038 ; // correction for RIF
    if (tt>tt0+tdil || od>od_crit) { // dilution
      tt0=tt ;
      for (j=0;j<L;j++) n[j]=binomial_2(dil,n[j]) ;
      if (tt>t_ant && c>1) { rif = rif*dil + RIF_rif0*(1 - dil); dec = dec*dil; c++ ; }
    }

    for (j=0;j<L;j++) {
			//double d=0 ; //0.01 ;
			double b ;
			if (tt<t_ant) b=gr[0][j] ; else b=gr[1][j] ;
			//if (tt>t_ant && j==0) b=gr[0][0]*exp(-10*(tt-t_ant)) ; // small lag - does not help much in shifting T_regrowth to longer times...

      QWORD db=binomial_2(b*dt,n[j]) ;
      n[j]+=db ; //genotypes[c->gen]->number+=db ; 
    }
    ntot=0LL ;
    for (j=0;j<L;j++) ntot+=n[j] ;
    od=ntot/(n_to_OD*vol) ;
    //if (tt>t_ant) od+=rif_od0+(tt-tt0)*rif_slope ; // correction for RIF
    // new, mechanistic correction for RIF:
    if (tt>t_ant) {
    	if (c==1) { rif=RIF_rif0 ; dec=0 ; c=2 ; }
 			if (c>1) {
 				od += rif + 3.3*dec; // in this simulation we add RIF contribution to the OD to make it consistent with how the experimental data are processed
 				rif -= RIF_lambda*rif*dt;
 				dec += RIF_lambda*rif*dt;
 			}
    }

    tts.push_back(tt) ; ods.push_back(od) ;
    ns.push_back(n[0]) ; nr.push_back(n[1]) ;
//    for (j=0;j<L;j++) cout <<n[j]<<"\t" ; cout <<endl ;
  } while (ntot>0LL && tt<tmax) ;
  return ntot ;
}


int main(int argc, char *argv[])
{
	double res_frac,res_frac_cv,n_ini_cv;
  int i,j,k;
	if (argc!=18) { err(" Error:: arguments needed: NAME vol od_ini gr0 gr1before gr1after t_ant t_max t_regrowth od_regrowth nsamples res_frac res_frac_cv n_ini_cv RIF_rif0 RIF_lambda dil\nProgram terminated. \n"); } 
  else { 
    NAME=argv[1] ;
    vol=atof(argv[2]) ;
    float od_ini=atof(argv[3]) ; n_ini=QWORD(n_to_OD*od_ini*vol) ;
    cout <<"od_ini="<<od_ini<<endl;
    cout <<"n_ini="<<n_ini<<endl ;
    birthS=atof(argv[4]) ; cout <<"birthS="<<birthS<<endl;
    birthRbefore=atof(argv[5]) ; cout <<"birthRbefore="<<birthRbefore<<endl;
    birthRafter=atof(argv[6]) ; cout <<"birthRafter="<<birthRafter<<endl;
    t_ant=atof(argv[7]) ; cout <<"t_ant="<<t_ant<<endl ;
    tmax=atof(argv[8]) ; cout <<"tmax="<<tmax<<endl ;
    t_regrowth=atof(argv[9]) ; cout<<"t_regrowth="<<t_regrowth<<endl;
    od_regrowth=atof(argv[10]) ; cout<<"od_regrowth="<<od_regrowth<<endl ;
    samples=atoi(argv[11]) ; 
    res_frac=atof(argv[12]) ; cout<<"res_frac="<<res_frac<<endl ;
    res_frac_cv=atof(argv[13]) ; cout<<"res_frac_cv="<<res_frac_cv<<endl ;
    n_ini_cv=atof(argv[14]) ; cout<<"n_ini_cv="<<n_ini_cv<<endl ;
    RIF_rif0=atof(argv[15]) ; cout<<"RIF_rif0="<<RIF_rif0<<endl ;
    RIF_lambda=atof(argv[16]) ; cout<<"RIF_lambda="<<RIF_lambda<<endl ;
    dil=atof(argv[17]) ; cout<<"dil="<<dil<<endl ;
  }
  _srand48(24) ;

	char txt[256];
	sprintf(txt,"std_var_distribution_%s.dat",NAME) ;
  ofstream data(txt) ;
	sprintf(txt,"std_var_n_%s.dat",NAME) ;
  ofstream ndata(txt) ;
  
	for (int sam=0;sam<samples;sam++) {
		if (sam%100==0) cout<<sam<<" ";

		//cout <<gama_res<<" " ;
		float rfrac=res_frac*exp(res_frac_cv*gauss()) ;
		QWORD nini=n_ini*exp(n_ini_cv*gauss()) ;
	  QWORD n=run_once(nini,rfrac) ; 
    //cout <<n<<" " ; 
    if (n>n_to_OD*vol*od_crit*0.5) { 
  		//cout <<"t=" ;  
  		for (i=tts.size()-1;i>0;i--) if (ods[i]<od_regrowth) break ;
  		for (j=0;j<tts.size();j++) if (ods[j]>od_crit) break ;
  		for (k=j;k<tts.size();k++) if (ods[k]<0.01) break ;
  		//cout <<tts[i] ;
  		//	time_to_regrowth	time_to_first_od_crit  time_to_kill	res_frac	nini
  		data<<tts[i]<<" "<<tts[j]<<" "<<tts[k]<<" "<<rfrac<<" "<<nini<<endl ;
  		
  		if (fabs(tts[i]-t_regrowth)<0.1 || tmax<t_ant || t_regrowth==100) { // save the data in case of a good fit
				cout<<"*" ;
				for (i=0;i<tts.size();i+=4) ndata<<tts[i]<<"\t"<<ns[i]<<"\t"<<nr[i]<<"\t"<<ods[i]<<endl ;
			}
    } else {
    	for (j=0;j<tts.size();j++) if (ods[j]>od_crit) break ;
    	for (k=j;k<tts.size();k++) if (ods[k]<0.01) break ;
      data<<"infty "<<tts[j]<<" "<<tts[k]<<" "<<rfrac<<" "<<nini<<endl;
    }
    //cout<<endl;
//    cout <<"gama="<<gama_res<<" nsurv="<<nsurv<<" "<<1.*tevol/nsurv<<endl ;
//    data <<gama_res<<"\t"<<nsurv<<"\t"<<1.*tevol/nsurv<<endl ;
  }
  data.close() ;
  ndata.close() ;
  //system("pause") ;  

}
