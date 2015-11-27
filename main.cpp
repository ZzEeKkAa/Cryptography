#include "config.h"
#include "help_functions.h"
#include "factorization.h"
#include "discrete_logarithm.h"


int main(){
    srand(time(0));
    vector<ll> prime;
    BuildSieveOfEratosthenes(35765,prime);

    /** Discrete Logarithms */

    ll p,s,n,a,b;

    a=254,b=871;
    p=7,s=3; n=pow(p,s)+1; //344//
    n=prime[234];

    //a=2,b=12;
    //n=19;

    //a=41,b=76;
    //p=2,s=8; n=pow(p,s)+1; //257//



    /*cout<<"Generators"<<endl;
    int i,j,t;
    for(i=2; i<n-1; ++i){
        if(pow(i,2*3*13,n)!=1 &&
           pow(i,2*3*19,n)!=1 &&
           pow(i,2*19*13,n)!=1 &&
           pow(i,3*19*13,n)!=1)
            cout<<i<<endl;
    }
    cout<<"End"<<endl;*/

    cout<<"log_"<<a<<"("<<b<<") = "<<GetLogarithmPrimitive(a,n,b)        <<" (mod "<<n<<")"<<endl;
    cout<<"log_"<<a<<"("<<b<<") = "<<GetLogarithmByShenks(a,n,b)         <<" (mod "<<n<<")"<<endl;
    cout<<"log_"<<a<<"("<<b<<") = "<<GetLogarithmIndexing(a,n,b)         <<" (mod "<<n<<")"<<endl;
    cout<<"log_"<<a<<"("<<b<<") = "<<GetLogarithmByPollard(a,n,b)        <<" (mod "<<n<<")"<<endl;
    //cout<<"log_"<<a<<"("<<b<<") = "<<GetLogarithmByPohligHellman(a,b,p,s)<<" (mod "<<n<<")"<<endl;/**/


    /** Factorization */

    n = prime[234]*prime[123]*prime[783]*prime[543]*16;
    //n=63214;
    //n=115249ll * 112909ll;

    vector<ll> factors;

    cout<<n<<endl;

    /*cout<<GetFactorByPollard(n)<<endl;
    cout<<GetFactorByDixon(n)<<endl;
    cout<<GetFactorQuadraticSieve(n)<<endl;*/


    GetAllFactors(n,GetFactorFerma,factors);
    //GetAllFactors(n,GetFactorByLenstra,factors);
    sort(factors.begin(),factors.end());
    for(auto& e: factors)
        cout<<e<<" ";
    cout<<endl;/**/

    return 0;
}
