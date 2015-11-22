#include "help_functions.h"

ll sqrt(ll a){
    return ll(sqrtl(ld(a)));
}

ll gcd(ll a, ll b){
    while(a>0 && b>0)
        if(a>b) a%=b; else b%=a;
    return a+b;
}

ll Random(ll a,ll b){
    ll r=rand()*(INT_MAX>>4)+rand()%(INT_MAX>>4);
    return a+r%(b-a);
}

ll pow(ll a, ll x, ll p){
    ll ans=1, t=a;
    for(;x;x=x>>1,t*=t,t%=p)
        if(x&1)
            ans*=t, ans%=p;
    return ans;
}

ll gcd(ll a, ll b, ll&x, ll&y){
    if(a<b) return gcd(b,a,y,x);
    // ax+by=g // a>b
    // ((a/b)*b + a%b) x + by = g
    // a%b x + b(y+ (a/b) *x) = g
    // x1 = x
    // y1 = y+ (a/b) x
    // y  = y1-(a/b) x

    if(b==0){
        x=1, y=0;
        return a;
    }
    ll x1,y1,g;
    g=gcd(b,a%b,y1,x1);
    x=x1;
    y=y1-(a/b)*x1;
    return g;
}

ll inverse(ll a, ll p){
    ll x, y;
    gcd(a,p,x,y);
    return x;
}

bool isQuadraticResidue(ll a, ll p){
    // a^((p-1)/2) == 1 mod p
    ll tp=(p-1)>>1;
    ll ans=1;
    for(;tp; tp=tp>>1){
        if(tp&1){
            ans*=a;
            ans%=p;
        }
        a*=a;
        a%=p;
    }
    return ans==1;
}

void BuildSieveOfEratosthenes(int to, vector<ll>& p){
    // Prime number generator

    vector<int> sieve;
    sieve.assign(to,1);

    for(int i=2; i<to; ++i){
        if(sieve[i]){
            p.push_back(i);
            for(int j=i+i; j<to; j+=i){
                sieve[j]=0;
            }
        }
    }
}

void SolveGauss(vector<vector<int> >& matrix, vector<int>& m_c, vector<int>& m_r){
    if(!matrix.size()) return;
    m_c.assign(matrix[0].size(),-1);   // return main element [row]=column
    m_r.assign(matrix.size(),-1); // return main element [column]=row
    for(int i=0; i<matrix[0].size(); ++i){ // rows
        for(int j=0; j<matrix.size(); ++j){ // columns
            if(matrix[j][i]%2!=0){
                m_r[j]=i;
                m_c[i]=j;
                for(int i0=0; i0<matrix[0].size(); ++i0) if(i!=i0 && matrix[j][i0]%2!=0)
                    for(int j0=0; j0<matrix.size(); ++j0)
                        matrix[j0][i0]+=matrix[j0][i];
                break;
            }
        }
    }
}

bool isSmooth(ll a,vector<ll>& B, vector<ll>& powers){
    powers.assign(B.size(),0);
    ll d;
    for(int i=0; i<B.size(); ++i){
        d=a/B[i];
        if(d*B[i]==a){
            ++powers[i];
            a=d;
            --i;
        }
    }
    return a==1;
}



