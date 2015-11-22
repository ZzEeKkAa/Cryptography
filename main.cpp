#include "config.h"
#include "help_functions.h"
#include "factorization.h"
#include "discrete_logarithm.h"


int main(){
    srand(time(0));
    //vector<ll> p;
    //BuildSieveOfEratosthenes(10000,p);
    //cout<<324<<"x = 1 (mod "<<p[300]<<")"<<endl;
    //cout<<inverse(16,18)<<endl;

    //cout<<GetLogarithmPrimitive(2,19,9)<<endl;
    //cout<<GetLogarithmByShenks(2,19,9)<<endl;
    //for(int i=0; i<40; ++i)
    ll p,s,n,a,b;

    a=23,b=87;
    p=7,s=3; n=pow(p,s)+1; //344*/

    /*a=23,b=46;
    p=2,s=8; n=pow(p,s)+1; //257*/

    cout<<"log_"<<a<<"("<<b<<") = "<<GetLogarithmPrimitive(a,n,b)        <<" (mod "<<n<<")"<<endl;
    cout<<"log_"<<a<<"("<<b<<") = "<<GetLogarithmByShenks(a,n,b)         <<" (mod "<<n<<")"<<endl;
    cout<<"log_"<<a<<"("<<b<<") = "<<GetLogarithmByPollard(a,n,b)        <<" (mod "<<n<<")"<<endl;
    cout<<"log_"<<a<<"("<<b<<") = "<<GetLogarithmByPohligHellman(a,b,p,s)<<" (mod "<<n<<")"<<endl;

    //srand(time(0));
    //vector<ll> p;
    //BuildSieveOfEratosthenes(35765,p);

    //cout<<p[232]*p[350]<<endl;
    //cout<<GetFactorQuadraticSieve(p[632]*p[750])<<endl;
    //cout<<GetFactorQuadraticSieve(p[232]*p[350])<<endl;
    //cout<<GetFactorByDixon(p[632]*p[750])<<endl;

    //cout<<GetFactorByPollard(p[632]*p[750])<<endl;

    //cout<<GetFactorByPollard(p[232]*p[350])<<endl;

    //cout<<GetFactorByDixon(p[232]*p[350])<<endl;
    //cout<<GetFactorQuadraticSieve(65237)<<endl;
    //cout<<GetFactorQuadraticSieve(15347)<<endl;

    return 0;
}
