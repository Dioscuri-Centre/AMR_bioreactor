#include <windows.h>
#include <stdio.h>
#include <tchar.h>
#include <string>
using namespace std;

int _tmain( int argc, TCHAR *argv[] )
{


    if( argc < 5 )
    {
        printf("Usage: %s [1]=[executable] [2]=[no. of replicates] [3]=[no. of args per replicate] [4...]=[all args. for all replicates]  \n", argv[0]);
        return -1;
    }
    
    int nrep=atoi(argv[2]) ;
    int n=atoi(argv[3]) ;
    
    char *cmdline=GetCommandLineA();
    printf("command line: %s\n",cmdline);
    if (argc!=4+nrep*n) {
    	printf("no of arguments does not make sense\n") ;
    	return -1 ;
    }
    
    string *args=new string[nrep] ;
    for (int i=0;i<nrep;i++) {
    	args[i]=argv[1] ;
			for (int j=4+i*n;j<4+i*n+n;j++) args[i]=args[i]+" "+argv[j] ;
			printf("args to be passed to replicate %d: %s\n",i,args[i].c_str());
		}
		
//		return -1 ;
    STARTUPINFO si[64];
    PROCESS_INFORMATION pi[64];

		for (int i=0;i<nrep;i++) {
			
    	ZeroMemory( &si[i], sizeof(si[i]) );
    	si[i].cb = sizeof(si[i]);
    	ZeroMemory( &pi[i], sizeof(pi[i]) );


	  // Start the child process. 
	    if( !CreateProcess( NULL,   // No module name (use command line)
	        args[i].c_str(),        // Command line
	        NULL,           // Process handle not inheritable
	        NULL,           // Thread handle not inheritable
	        FALSE,          // Set handle inheritance to FALSE
	        0,              // No creation flags
	        NULL,           // Use parent's environment block
	        NULL,           // Use parent's starting directory 
	        &si[i],            // Pointer to STARTUPINFO structure
	        &pi[i] )           // Pointer to PROCESS_INFORMATION structure
	    ) 
	    {
	        printf( "CreateProcess failed (%d).\n", GetLastError() );
	        return -1;
	    }
  	}

    // Wait until child processes exit
    for (int i=0;i<nrep;i++) {
    	WaitForSingleObject( pi[i].hProcess, INFINITE );

    // Close process and thread handles. 
	    CloseHandle( pi[i].hProcess );
  	  CloseHandle( pi[i].hThread );
  	}
    
    return 0 ;
}
