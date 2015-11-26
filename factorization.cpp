#include "factorization.h"
#include "help_functions.h"

ll GetFactorFerma(ll n){
    if(n==2) return 1;
    if(n%2==0) return 2;
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

ll GetFactorByPollard(ll n){
    const int MAX_ITERATIONS = 10000;
    ll x = Random(1,n);

    ll y = 1,g;
    for(int i = 0, stage = 2; i<MAX_ITERATIONS; ++i){
        #ifdef DEBUG
        std::cout<<i<<std::endl;
        #endif // DEBUG
        g = gcd(n, abs(x-y));
        if(g>1 && g<n) return g;
        if (i == stage ){
            y = x;
            stage = stage<<1;
        }
        x = (x*x + 1)%n;
    }

    return gcd(n, abs(x-y));
}

ll GetFactorQuadraticSieve(ll n){
    vector<ll> p;
    BuildSieveOfEratosthenes(1000, p);
    // 1. Factor base
    vector<ll> B;

    ll sq_n=sqrt(n);
    const int h = 12;
    const int v= 25;

    // Creating factor base
    int r;
    for(int i=0; B.size()<h && i<p.size(); ++i){
        if(isQuadraticResidue(n,p[i])){
            B.push_back(p[i]);
            #ifdef DEBUG
            cout<<p[i]<<" ";
            #endif // DEBUG
        }
        r=i;
    }
    #ifdef DEBUG
    cout<<"total checked: "<<r<<endl;
    #endif // DEBUG

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

        //4. Checking for B-smooth
        for(int i=0; i<B.size(); ++i){
            d=a/B[i];
            if(d*B[i]==a){
                a=d;
                ++vect[i];
                --i;
            }
        }
        if(a==1){
            for(auto it=vect.begin(); it<vect.end(); ++it)
                if(*it%2) {
                    matrix.push_back(vect);
                    bs.push_back(b);
                    #ifdef DEBUG
                    cout<<"Adding to matrix "<<b<<endl;
                    #endif //DEBUG
                    break;
                }
        }
    }while(matrix.size()<v); //5. repeating while size of matrix less than v

    #ifdef DEBUG

    std::cout<<"Gauss"<<std::endl;

    for(int i=0; i<vect.size(); ++i){ // rows
        for(int j=0; j<matrix.size(); ++j){ // columns
            std::cout<<matrix[j][i]%2<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;

    #endif // DEBUG

    // 6. Gauss
    bkp.assign(matrix.begin(), matrix.end());
    vector<int> m_r,m_c;
    SolveGauss(matrix, m_c,m_r);


    #ifdef DEBUG

    for(int i=0; i<vect.size(); ++i){ // rows
        for(int j=0; j<matrix.size(); ++j){ // columns
            std::cout<<matrix[j][i]%2<<" ";
        }
        std::cout<<std::endl;
    }

    #endif // DEBUG


    ll x,y;
    vector<int> depended;
    for(int i=0; i<matrix.size(); ++i){
        if(m_r[i]==-1){
            depended.resize(0);
            for(int j=0; j<vect.size(); ++j){
                if(matrix[i][j]%2!=0){
                    depended.push_back(m_c[j]);
                }
            }

            //for(auto it=depended.begin(); it<depended.end(); ++it) cout<<*it<<" "; cout<<endl;

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
                //cout<<"here"<<endl;


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

ll GetFactorByDixon(ll n){
    // https://ru.wikipedia.org/wiki/%D0%90%D0%BB%D0%B3%D0%BE%D1%80%D0%B8%D1%82%D0%BC_%D0%94%D0%B8%D0%BA%D1%81%D0%BE%D0%BD%D0%B0

    // 1. Factor base
    vector<ll> p;
    BuildSieveOfEratosthenes(1000,p);

    ll sq_n=sqrt(n);
    double lg=log(ld(n));
    double L=exp(sqrt(lg*log(lg)));

    int M=sqrt(L);
    #ifdef DEBUG
    cout<<"M = "<<M<<endl;
    #endif // DEBUG

    int m,l=1,r=p.size()-1,h;
    while(l<r){
        m=(l+r+1)/2;
        if(p[m]>M) r=m-1;
        else       l=m;
    }
    h=l+1;

    #ifdef DEBUG
    cout<<"h: "<<h<<endl;
    #endif // DEBUG

    vector<vector<int> > matrix, bkp;
    vector<int> vect;
    vector<ll> bs;

    ll a,b,d,t;
    b=sq_n;
    do{
        //2. Generating b
        b=Random(sq_n,n);
        //b=b+BigInteger("1",b.GetSystem());
        //std::cout<<"checking "<<b<<std::endl;
        //3. Calculating a
        a=(b*b) % n;
        vect.assign(h,0);

        //4. Checking for smooth
        for(int i=0; i<h; ++i){
            d=a/p[i];
            if(d*p[i]==a){
                a=d;
                ++vect[i];
                --i;
            }
        }

        if(a==1){
            for(auto it=vect.begin(); it<vect.end(); ++it)
                if(*it%2) {
                    matrix.push_back(vect);
                    bs.push_back(b);
                    #ifdef DEBUG
                    cout<<"Adding to matrix "<<b<<endl;
                    #endif
                    break;
                }
        }
    }while(matrix.size()<h); //5. repeating while size of matrix less than h

    // 6. Gauss
    #ifdef DEBUG
    std::cout<<"Gauss"<<std::endl;

    for(int i=0; i<vect.size(); ++i){ // rows
        for(int j=0; j<matrix.size(); ++j){ // columns
            std::cout<<matrix[j][i]%2<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    #endif // DEBUG


    bkp.assign(matrix.begin(), matrix.end());
    std::vector<int> m_c, m_r;

    SolveGauss(matrix,m_c,m_r);

    #ifdef DEBUG
    for(int i=0; i<vect.size(); ++i){ // rows
        for(int j=0; j<matrix.size(); ++j){ // columns
            std::cout<<matrix[j][i]%2<<" ";
        }
        std::cout<<std::endl;
    }
    #endif // DEBUG

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

            if(depended.size()){
                vect.assign(bkp[i].begin(), bkp[i].end());

                x=bs[i];
                for(auto it=depended.begin(); it<depended.end(); ++it){
                    x=x*bs[*it];
                    x%=n;
                    for(int i=0; i<vect.size(); ++i)
                        vect[i]+=bkp[*it][i];
                }
                x%=n;

                y=1;
                for(int i=0; i<vect.size(); ++i)
                    for(int j=0; j<vect[i]; ++j)
                        y*=p[i], y%=n;
                y%=n;


                ll s=(x+y)%n; if(s<0) s+=n;
                ll g=gcd(s,n);

                if(x!=y && x!=n-y && g>1)
                     return g;
            }
        }
    }

    return GetFactorByDixon(n);
}

ll GetFactorByLenstraBase(ll a, ll b, ll x, ll y, ll w, ll v, ll n, vector<ll> const& p){
    #ifdef DEBUG
    #define DEBUG

    cout<<"a = "<<a<<"; ";
    cout<<"b = "<<b<<"; ";
    cout<<"x = "<<x<<"; ";
    cout<<"y = "<<y<<"; ";
    cout<<"w = "<<w<<"; ";
    cout<<"v = "<<v<<"; ";
    cout<<endl;

    #endif // DEBUG

    ll m;
    ll rp;
    vector<ll> k;

    for(auto r =p.begin(); *r<w; ++r){
        m=log(v)/log(*r);
        rp=pow(*r,m);
        while(rp<=v) rp*=*r;
        while(rp>v) rp/=*r;
        k.push_back(rp);
    }

    sort(k.begin(),k.end());

    point P = make_pair(x,y);
    for(auto& ki: k){
        P=GroupMul(a,b,ki,P,n);
        if(P.first==LLONG_MIN) return P.second;
    }

    return 1;
}

ll GetFactorByLenstra(ll n){
    vector<ll> p;
    BuildSieveOfEratosthenes(10000,p);

    ll a,b,x,y,v,w,lg,fac;
    int i=0;
    do{
        a=rand()%n;
        x=rand()%n;
        y=rand()%n;
        b=(((y*y)%n) - ((((x*x)%n)*x)%n) - ((a*x)%n))%n; if(b<0) b+=n;

        lg=log(n);
        w = sqrt(exp(sqrt(lg*log(lg))))+1;
        v = 1000; // !!!!!!

        fac = GetFactorByLenstraBase(a,b,x,y,w,v,n,p);
        ++i;
    }while( !(fac<n && fac>1) && i<100  );
    //cout<<fac<<endl<<LLONG_MAX<<endl<<endl
    if(fac<n && fac>1) return fac;
    return 1;
}

void GetAllFactors(ll n, ll(*GetFactor)(ll), vector<ll>&factors){
    ll a = (*GetFactor)(n);
    if(a==1){
        factors.push_back(n);
    }else{
        GetAllFactors(a,GetFactor,factors);
        GetAllFactors(n/a,GetFactor,factors);
    }
}


