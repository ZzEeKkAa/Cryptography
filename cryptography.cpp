#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef long double ld;

ll GetFactorFerma(ll n);
ll GetFactorQuadraticSieve(ll n);
inline ll sqrt(ll);
bool isQuadraticResidue(ll,ll);
ll gcd(ll,ll);
void BuildSieveOfEratosthenes(int);
static vector<ll> p;
void SolveGauss(vector<vector<int> >& matrix, vector<int>& m_c, vector<int>& m_r);

int main(){
    BuildSieveOfEratosthenes(35765);

    //cout<<GetFactorQuadraticSieve(9524242539L)<<endl;
    //cout<<GetFactorQuadraticSieve(95832425L)<<endl;
    cout<<GetFactorQuadraticSieve(15347)<<endl;

    return 0;
}



ll GetFactorFerma(ll n){
    ll x,y;
    ll d;
    ll to=(n+1)/2;


    for(x=sqrt(n), d=x*x-n; x<to; d+=(x<<1)+1,++x){
        y=sqrt(d);
        if(y*y==d){
            return x-y;
        }
    }

    return 1;
}

ll GetFactorQuadraticSieve(ll n){
    // 1. Factor base
    static vector<ll> B;

    ll sq_n=sqrt(n);
    /*ld lg=log(n);
    ld L=exp(sqrt(lg*log(lg)));

    int M=sqrt(L);*/
    /* !! */const int h = 4;
    const int v= 5;

    // Creating factor base
    int r;
    for(int i=0; B.size()<h && i<p.size(); ++i){
        if(isQuadraticResidue(n,p[i])){
            B.push_back(p[i]);
            cout<<p[i]<<" ";
        }
        r=i;
    }
    cout<<"total checked: "<<r<<endl;


    vector<vector<int> > matrix, bkp;
    vector<int> vect;
    vector<ll> bs;

    ll a,b,d,t;
    b=sq_n;
    ll i=0;
    do{
        //2. Generating b
        b+=i;
        if(i>0) ++i; else --i; i*=-1;

        //3. Calculating a
        a=(b*b) - n;
        vect.assign(h,0);

        //4. Checking for smooth
        for(int i=0; i<B.size(); ++i){
            d=a/B[i];
            if(d*B[i]==a){
                a=d;
                ++vect[i];
                --i;
            }
        }
        if(d==0){
            std::cout<<"Matrix size: "<<matrix.size()<<"// i="<<i<<" // b="<<b<<std::endl;
            for(auto it=vect.begin(); it<vect.end(); ++it)
                if(*it%2) {
                    matrix.push_back(vect);
                    bs.push_back(b);
                    cout<<"good!"<<endl;
                    break;
                }
        }
    }while(matrix.size()<v); //5. repeating while size of matrix less than h+1
    std::cout<<"Gauss"<<std::endl;


    for(int i=0; i<vect.size(); ++i){ // rows
        for(int j=0; j<matrix.size(); ++j){ // columns
            std::cout<<matrix[j][i]%2<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;


    // 6. Gauss
    bkp.assign(matrix.begin(), matrix.end());
    vector<int> m_r,m_c;
    SolveGauss(matrix, m_c,m_r);


    for(int i=0; i<vect.size(); ++i){ // rows
        for(int j=0; j<matrix.size(); ++j){ // columns
            std::cout<<matrix[j][i]%2<<" ";
        }
        std::cout<<std::endl;
    }


    ll x,y;
    std::vector<int> depended;
    for(int i=0; i<matrix.size(); ++i){
        if(m_r[i]==-1){
            depended.resize(0);
            for(int j=0; j<vect.size(); ++j){
                if(matrix[i][j]%2!=0){
                    depended.push_back(m_c[j]);
                }
            }

            for(auto it=depended.begin(); it<depended.end(); ++it) cout<<*it<<" "; cout<<endl;

            if(depended.size()){
                vect.assign(bkp[i].begin(), bkp[i].end());

                x=bs[i];
                for(auto it=depended.begin(); it<depended.end(); ++it){
                    x=(x*bs[*it])%n;
                    for(int i=0; i<vect.size(); ++i)
                        vect[i]+=bkp[*it][i];
                }

                y=1;
                for(int i=0; i<vect.size(); ++i)
                    for(int j=0; j<vect[i]; j+=2)
                        y=(y*B[i])%n;
                cout<<"here"<<endl;


                if(x!=y && x!=n-y){
                    ll ans=gcd(x+y,n);
                    if(ans!=1)
                        return gcd(x+y,n);
                }
            }
        }
    }

    return 1;
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

ll gcd(ll a, ll b){
    cout<<"gcd: "<<a<<" "<<b<<endl;
    while(a>0 && b>0)
        if(a>b) a%=b; else b%=a;
    cout<<"gcd: "<<a+b<<endl;
    return a+b;
}

inline ll sqrt(ll a){
    return ll(sqrtl(ld(a)));
}

void BuildSieveOfEratosthenes(int to){
    // Prime number generator

    vector<int> sieve;
    sieve.assign(to,1);

    for(int i=2; i<to; ++i){
        if(sieve[i]){
            //cout<<i<<" ";
            p.push_back(i);
            for(int j=i+i; j<to; j+=i){
                sieve[j]=0;
            }
        }
    }

    cout<<endl;

}
