// this program simulates a population of cells growing in a turbidostat, using the simplistic tau leaping algorithm
// cells can either replicate or mutate with probability mu
// antibiotic is injected at some time t_ant which causes cessation of growth of the wild type but a mutant can still reproduce
// simulation is stopped when either OD 0 or od_crit is reached
// the simulation is then repeated some number of times 
// this version allows one to use more than one initial strain, and each is simulated individually, 
// in addition, each strain can mutate into a variable number of mutants, each with a different growth rate
// each simulation is provided with an individual file of growth rates
// 		in the case of RIF, its contribution to the OD is calculated but not included, since the simulation is compared with RIF-corrected data

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

enum AB_TYPE {RIF=0,CIP=1,STR=2,GENERIC=3} ;
AB_TYPE AB_type ; // antibiotic type, will be selected from the params file

int samples=1000;
float vol=20 ; // volume [ml]
float mu_total=1e-7 ; // total mutation probability per replication
//const QWORD capacity=3e9 ; // cfu/ml
const double n_to_OD=2e8 ; // cfu/ml which would give OD_linear=1 in the turbidostat had we grown bacteria beyond OD=0.1
QWORD n_ini ; // initial number of cells
float t_ant ; 
const float tdil=0.5;
float dil ; // dilution factor - will be obtained from command line
const float od_crit=0.1 ;
float tmax ;
char *NAME,*name_input ;

int strains ; // no. of strains in the initial mix
int alleles ;	// how many resistant alleles can be generated
float *gr[2] ; // growth rate, [before=0/after=1][type]  
float t_regrowth, od_regrowth;

// S: 0..strains-1
// R: strains..strains*(alleles+1)-1

void err(const char *reason)
{
  cout <<reason<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif  
  exit(0) ;
}

void err(const char *reason, int a)
{
  cout <<reason<<": "<<a<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif    
  exit(0) ;
}

void err(const char *reason, char *a)
{
  cout <<reason<<": "<<a<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif    
  exit(0) ;
}

void err(const char *reason, double a)
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

vector <float> ods,tts,rifs,decs;
vector <QWORD> *ns,*nr ;
vector <QWORD> finaln ; // final occupations for all strains and alleles

float RIF_lambda = 0.65, RIF_rif0 = 0.018 ; // for modelling RIF, will be overridden

// for modelling CIP
double n_to_OD_filaments ; // used to simulate filaments, will be initialized at the start of each simulation
float CIP_time_to_stop_filamenting=0.0 ; // for modelling CIP, will be overridden
float CIP_lambda ;


void init()
{
	char txt[256];
	FILE *f=fopen(name_input,"r") ;
	fscanf(f,"%s",txt) ;
	switch (txt[0]) {
		case 'R': cout<<"RIF simulation"<<endl ; AB_type=RIF ; fscanf(f,"%f",&RIF_lambda) ; fscanf(f,"%f",&RIF_rif0) ; break ;
		case 'C': cout<<"CIP simulation"<<endl ; AB_type=CIP ; fscanf(f,"%f", &CIP_time_to_stop_filamenting) ; break ;
		case 'S': cout<<"STR simulation"<<endl ; AB_type=STR ; break ;
		case 'G': cout<<"generic antibiotic"<<endl ; AB_type=GENERIC ; break ;
		default: err("unknown antibiotic") ;
	}
	fscanf(f,"%d",&strains) ; fscanf(f,"%d",&alleles) ; 
	cout<<strains<<" "<<alleles<<endl ;
	if (strains<1 || strains>2 || alleles<1 || alleles>100) err("params wrong") ;
	gr[0]=new float[strains*(1+alleles)] ;
	gr[1]=new float[strains*(1+alleles)] ;
	ns=new vector<QWORD>[strains] ;
	nr=new vector<QWORD>[strains] ;
	for (int i=0;i<strains*(1+alleles);i++) fscanf(f,"%f",&gr[0][i]) ; 
	if (ferror(f) || feof(f)) err("error 1 when reading parameters") ;
	for (int i=0;i<strains*(1+alleles);i++) fscanf(f,"%f",&gr[1][i]) ; 
	if (ferror(f) || feof(f)) err("error 2 when reading parameters") ;
	fclose(f) ;
}

double tt, tt0 ;
int nsam=0;
QWORD run_once(int sam)
{
  int i,j,k,l;  
	ods.clear() ;
	finaln.clear() ;
	for (i=0;i<strains;i++) { ns[i].clear() ; nr[i].clear() ; }
	tts.clear() ;
  const int L=strains*(alleles+1);
  tt=tt0=0 ;
  QWORD n[L] ;
  for (i=0;i<L;i++) n[i]=0LL ;
  QWORD ntot=0 ;
  float od=0,od_add=0 ;

	// for strains=2:
	// n[0] = strain 1_S, n[1] = strain 2_S,
	// n[2] = strain 1_R1, n[3] = strain 2_R1,
	// n[4] = strain 1_R2, n[5] = strain 2_R3,
	// etc.

	QWORD n0=n_ini*(0.5+1.5*_drand48()) ; // initial number of cells
  for (i=0;i<strains;i++) {
		n[i]=n0/strains ;
		if (sam==0) {
			cout <<gr[0][i]<<" " ; //=birthS*(0.975+_drand48()*0.05) ;
			cout <<gr[1][i]<<" " ; //=0 ; 
		}
	}

	//cout <<"grr: " ;
	for (i=strains;i<strains*(alleles+1);i++) {
		if (sam==0) {
			cout<<gr[0][i]<<" " ; 
			cout<<gr[1][i]<<" " ; 
		}
		//cout <<gr[0][i]<<" " ;
	}
	if (sam==0) cout <<endl ;


  const double dt=1./int(1<<6) ;

	int c=1 ; // for modelling RIF
	float rif=0,dec=0 ; // for modelling RIF
	rifs.clear() ; decs.clear() ;
  
  // for modelling CIP
	n_to_OD_filaments=n_to_OD ; 
	CIP_lambda=1 ;

	      
  do {      // main loop 
    tt+=dt ; //cout <<genotypes.size()<<endl ;

    if (tt>tt0+tdil || od+od_add>od_crit) { // dilution
      tt0=tt ;
      for (j=0;j<L;j++) n[j]=binomial_2(dil,n[j]) ;
      if (AB_type==RIF && tt>t_ant && c>1) {  
				rif = rif*dil + RIF_rif0*(1 - dil); dec = dec*dil; c++ ;
			}
    }

    for (j=0;j<L;j++) { // no death in this model
			double b ;
			if (tt<t_ant) b=gr[0][j] ; else b=gr[1][j] ;
			if (AB_type==CIP && tt>t_ant && CIP_time_to_stop_filamenting>0 && j==0 && CIP_lambda>0) {
				n_to_OD_filaments/=1+dt*CIP_lambda*gr[0][0] ; 
				CIP_lambda=1-SQR((tt-t_ant)/CIP_time_to_stop_filamenting) ;
			}


      QWORD db=binomial_2(b*dt,n[j]) ;
      if (db>0) {         // replication
        int dbmut ; 
        if (j<strains) dbmut=binomial_2(mu_total,db) ; else dbmut=0 ;
        if (dbmut==0) {
          n[j]+=db ; //genotypes[c->gen]->number+=db ; 
        } else {
//            cout <<"dm="<<dbmut<<"\t" ;
          if (dbmut>db) err("dbmut") ;
          n[j]+=db-dbmut ; //c->number+=db-dbmut ; 
//              cout<<j<<"->"<<c->gen<<","<<c->n<<" " ;
					// j=0 should mutate into 2,4,6...
					// j=1 should mutate into 3,5,7,...

          n[strains+j+strains*int(_drand48()*alleles)]+=dbmut ;
              //cout <<nc.gen<<":" ;
//              cout <<"x"<<c->gen<<endl ;
        }
      }
    }
    ntot=0LL ;
    
    od=od_add=0 ;
    for (j=0;j<L;j++) {
			ntot+=n[j] ;
    	if (j>=strains) od+=n[j]/(n_to_OD*vol) ;
    	else od+=n[j]/(n_to_OD_filaments*vol) ; // the sensitive strain can filament
		}

    
    if (AB_type==RIF && tt>t_ant) {
    	if (c==1) { rif=RIF_rif0 ; dec=0 ; c=2 ; }
 			if (c>1) {
 				// od += rif + 3.3*dec;  // do not add to the exported data since the model will be compared with corrected data, but the RIF contribution will be still outputted for plotting example trajectories
 				od_add = rif + 3.3*dec; // also, OD from RIF is added to the OD when deciding when to dilute
 				rif -= RIF_lambda*rif*dt;
 				dec += RIF_lambda*rif*dt;
 			}
    }
    
    tts.push_back(tt) ; ods.push_back(od) ;
    rifs.push_back(rif) ; decs.push_back(dec) ; 
//    QWORD ntest=0 ;
		for (i=0;i<strains;i++) { 
			ns[i].push_back(n[i]) ; 
			QWORD nnr=0 ;
			for (j=0;j<alleles;j++) nnr+=n[strains+i+j*strains] ;
			nr[i].push_back(nnr) ; 
//			ntest+=nr[i][nr[i].size()-1]+ns[i][ns[i].size()-1] ;
		}
//		if (ntest!=ntot) err("!") ;
//    for (j=0;j<L;j++) cout <<n[j]<<"\t" ; cout <<endl ;
  } while (ntot>0LL && tt<tmax) ;
  
  //for (i=0;i<L;i++) cout <<n[i]<<" " ; cout <<endl ;
  for (i=0;i<L;i++) finaln.push_back(n[i]) ; 
  return ntot ;
}


int main(int argc, char *argv[])
{
	double mu_min,mu_max,ttol=0.01;
  int i,j,k;
	if (argc<14) { err(" Error:: arguments needed: NAME vol od_ini gr0 gr1 t_ant t_max t_regrowth od_regrowth nsamples mu_min mu_max dilf\nProgram terminated. \n"); } 
  else { 
    NAME=argv[1] ;
    vol=atof(argv[2]) ;
    float od_ini=atof(argv[3]) ; n_ini=QWORD(n_to_OD*od_ini*vol) ;
//    birthS=atof(argv[4]) ;
//    birthR=atof(argv[5]) ;
    t_ant=atof(argv[4]) ;
    tmax=atof(argv[5]) ;
    t_regrowth=atof(argv[6]) ;
    od_regrowth=atof(argv[7]) ;
    samples=atoi(argv[8]) ;
    mu_min=atof(argv[9]) ;
    mu_max=atof(argv[10]) ;
    dil=atof(argv[11]) ;
    name_input=argv[12] ;
  }
  if (argc<=14) { ttol=atof(argv[13]) ; cout<<"ttol="<<ttol<<endl ; }
  _srand48(24) ;

	char txt[256];
	sprintf(txt,"distribution_%s.dat",NAME) ;
  ofstream data(txt) ;
	sprintf(txt,"n_%s.dat",NAME) ;
  ofstream ndata(txt) ;
	sprintf(txt,"strains_%s.dat",NAME) ;
  ofstream sdata(txt) ;
  
  init() ;
  
	for (int sam=0;sam<samples;sam++) {
		if (sam%100==0) cout<<sam<<" ";
		mu_total=mu_min*pow(mu_max/mu_min,_drand48()) ;

		//cout <<gama_res<<" " ;
	  QWORD n=run_once(sam) ; 
    //cout <<n<<" " ; 
    if (n>n_to_OD*vol*od_crit*0.5) { 
  		//cout <<"t=" ;  
  		for (i=tts.size()-1;i>0;i--) if (ods[i]<od_regrowth) break ;
  		for (j=0;j<tts.size();j++) if (ods[j]>od_crit) break ;
  		for (k=j;k<tts.size();k++) if (ods[k]<0.01) break ;
  		//cout <<tts[i] ;
  		// gamma	time_to_regrowth	time_to_first_od_crit
  		data<<mu_total<<" "<<tts[i]<<" "<<tts[j]<<" "<<tts[k] ;
			for (j=0;j<strains;j++) data<<" "<<ns[j][ns[j].size()-1]<<" "<<nr[j][nr[j].size()-1] ;
			data<<endl;
  		
  		if (fabs(tts[i]-t_regrowth)<ttol || tmax<t_ant || t_regrowth==100) { // save the data in case of a good fit
				cout<<"*" ;
				for (i=0;i<tts.size();i++) {
					ndata<<tts[i]<<"\t" ; 
					for (j=0;j<strains;j++) ndata<<ns[j][i]<<"\t"<<nr[j][i]<<"\t" ;
					ndata<<ods[i] ; 
					ndata<<"\t"<<rifs[i]<<"\t"<<decs[i] ; 
					
					ndata<<endl ;
				}
			}
			
			for (j=0;j<finaln.size();j++) sdata<<finaln[j]<<" " ; sdata<<endl ;
			
    } else {
    	for (j=0;j<tts.size();j++) if (ods[j]>od_crit) break ;
    	for (k=j;k<tts.size();k++) if (ods[k]<0.01) break ;
      data<<mu_total<<" infty "<<tts[j]<<" "<<tts[k] ; 
			for (j=0;j<strains;j++) data<<" "<<ns[j][ns[j].size()-1]<<" "<<nr[j][nr[j].size()-1] ;
			data<<endl;
    }
    //cout<<endl;
//    cout <<"gama="<<gama_res<<" nsurv="<<nsurv<<" "<<1.*tevol/nsurv<<endl ;
//    data <<gama_res<<"\t"<<nsurv<<"\t"<<1.*tevol/nsurv<<endl ;
  }
  data.close() ;
  ndata.close() ;
  sdata.close() ;
  //system("pause") ;  

}
