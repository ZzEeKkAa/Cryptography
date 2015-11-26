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
    x%=p; if(x<0) x+=p;
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


point GroupAdd(ll a, ll b, point P, point Q, ll n){
    point O = make_pair(LLONG_MAX,LLONG_MAX);
    if(P==O) return Q;
    if(Q==O) return P;
    if(P.first == Q.first && P.second==-Q.second) return O;

    ll s,d;
    if(P.first == Q.first){
        if(P.second == Q.second){
            s=((3*P.first)%n)*P.first + a;
            d=2*P.second;
        }else
            return make_pair(LLONG_MIN, LLONG_MIN);

    } else{
        s=Q.second-P.second;
        d=Q.first-P.first;
    }


    if(d<0) s*=-1, d*=-1;

    ll x, y;
    ll g = gcd(d,n,x,y);

    if(g!=1) return make_pair(LLONG_MIN, g);

    s*=x; s%=n; if(s<0) s+=n;

    x=s*s-P.first-Q.first; x%=n; if(x<0) x+=n;
    y=P.second+s*(x-P.first);

    return make_pair(x,y);

}

point GroupMul(ll a, ll b, ll k, point P, ll n){
    point t=P;
    point ans=make_pair(LLONG_MAX,LLONG_MAX);

    while(k) {
        if(k&1){
            ans=GroupAdd(a,b,ans,t,n);
            if(ans.first==LLONG_MIN) return ans;
        }
        t=GroupAdd(a,b,t,t,n);
        if(t.first==LLONG_MIN) return t;
        k=k>>1;
    }

    return ans;
}

bool AddToGaussSystem(vector<vector<ll> >& matrix, vector<int>& m_c, vector<int>& m_r, vector<ll>const&vect, ll p){
    const int n = matrix.size()+1;
    const int m = vect.size();

    ll d;

    int i;

    #ifdef DEBUG
    cout<<"__________________________"<<endl;

    for(auto& row: matrix){
        for(auto& el: row)
            printf("%4lld ",el);
        cout<<endl;
    }

    for(auto& el: m_r)
        printf("%4d ",el);
    cout<<endl;
    #endif // DEBUG


    matrix.push_back(vect);

    #ifdef DEBUG
    cout<<endl<<"Candidate:"<<endl;
    for(auto& el: matrix[n-1])
        printf("%4d ",el);
    cout<<endl;
    #endif // DEBUG

    // Removing existed base elements from new vector
    for(int j=0; j<m-1; ++j){
        if(m_r[j]!=-1 && matrix[n-1][j]!=0){
            i=m_r[j];

            d = p - matrix[n-1][j];

            for(int j0=0; j0<m; ++j0){
                matrix[n-1][j0]+=d*matrix[i][j0]; matrix[n-1][j0]%=p;
            }
        }
    }

    #ifdef DEBUG
    cout<<endl<<"After removed basis"<<endl;

    for(auto& e: matrix[n-1])
        printf("%4lld ",e);
    cout<<endl;
    #endif // DEBUG

    int jb=-1;

    // Checking if new vector is not repeat and if yes normalizing basis element
    for(int j=0; j<m-1; ++j){
        if(matrix[n-1][j]!=0 && m_r[j]==-1){
            jb=j;
            m_c.push_back(j);
            m_r[j]=n-1;

            ll x, y;

            ll g = gcd(matrix[n-1][j],p,x,y);
            if(g!=1) {
                m_r[j]=-1;
                jb=-1;
                break;
            }

            x%=p; if(x<0) x+=p;

            d=x;

            for(int j0=0; j0<m; ++j0){
                matrix[n-1][j0]*=d; matrix[n-1][j0]%=p;
            }

            break;
        }
    }

    #ifdef DEBUG
    cout<<"After inversing"<<endl;

    for(auto& e: matrix[n-1])
        printf("%4lld ",e);
    cout<<endl;
    #endif // DEBUG

    // Removing new base element from others vectors
    if(jb!=-1){
        for(int i=0; i<n-1; ++i){
            if(matrix[i][jb]!=0){
                d=p-matrix[i][jb];
                for(int j=0; j<m; ++j){
                    matrix[i][j]+=d*matrix[n-1][j];
                    matrix[i][j]%=p;
                }
            }
        }

    }else matrix.pop_back();

    return jb!=-1;

}

