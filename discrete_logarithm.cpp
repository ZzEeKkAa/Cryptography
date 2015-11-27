#include "discrete_logarithm.h"
#include "help_functions.h"

ll GetLogarithmPrimitive(ll g, ll n, ll b){
    ll x=0, t=1;
    for(; x<n; ++x,t*=g,t%=n)
        if(b==t) return x;

    return -1;
}

ll GetLogarithmByShenks(ll g, ll n, ll b){
    ll a=sqrt(n)+1;
    ll ga=pow(g,a,n);
    ll t;

    #ifdef DEBUG
    cout<<"a="<<a<<endl;
    #endif // DEBUG

    vector<ll> l1(1,1),l2,ans(a+1);

    t=ga;
    for(int i=0; i<a; ++i){
        l1.push_back(t); t*=ga; t%=n;
    }

    t=b;
    for(int i=0; i<a; ++i){
        l2.push_back(t); t*=g; t%=n;
    }

    sort(l1.begin(), l1.end());
    sort(l2.begin(), l2.end());

    #ifdef DEBUG
    for(auto& e:l1) cout<<e<<" "; cout<<endl;
    for(auto& e:l2) cout<<e<<" "; cout<<endl;
    #endif // DEBUG

    auto it=set_intersection(l1.begin(), l1.end(), l2.begin(), l2.end(), ans.begin());
    ans.resize(it-ans.begin());

    if(ans.size()) {
        ll x=0;
        t=ga;
        for(int i=0; i<a; ++i){
            if(t==ans[0]) {x+=(i+1)*a; break;}
            t*=ga; t%=n;
        }

        t=b;
        for(int i=0; i<a; ++i){
            if(t==ans[0]) {x-=i; break;}
            t*=g; t%=n;
        }
        x%=n-1; if(x<0) x+=n-1;
        return x;
    }
    return -1;
}

ll GetLogarithmByPollard(ll a, ll n, ll b){
    vector<ll> x(1,1),pa(1,0),pb(1,0);
    int i;
    for(i=0; !(i && i%2==0 && x[i]==x[i>>1]); ++i){
        if(pa[i]>(LLONG_MAX>>4) || pb[i]>(LLONG_MAX>>4) ) return GetLogarithmByPollard(a, n, b);

        switch(rand()%3){
        case(0):
            x.push_back((a*x[i])%n);
            pa.push_back(pa[i]+1);
            pb.push_back(pb[i]);
            break;
        case(1):
            x.push_back((x[i]*x[i])%n);
            pa.push_back(pa[i]<<1);
            pb.push_back(pb[i]<<1);
            break;
        case(2):
            x.push_back((b*x[i])%n);
            pa.push_back(pa[i]);
            pb.push_back(pb[i]+1);
            break;
        }
        //cout<<i+1<<": "<<x[i+1]<<endl;
        //cout<<(i+1)%2<<" "<<x[ (i+1)>>1 ]<<" "<<x[i+1]<<endl;
    }

    #ifdef DEBUG
    cout<<x[i]<<" "<<x[i>>1]<<endl;
    for(auto& e:x) printf("%5lld",e); printf("\n");
    for(auto& e:pa) printf("%5lld",e); printf("\n");
    for(auto& e:pb) printf("%5lld",e); printf("\n");
    #endif // DEBUG
    // a^pa[i>>1] * b^pb[i>>1] = a^pa[i] * b^pb[i]
    // a^(pa[i>>1] - pa[i]) = b^(pb[i]-pb[i>>1])

    // a^t1=b^t2
    // t1=t2 *l (mod)
    // g = gcd(t2,n-1)
    // g|t1

    ll t1=-(pa[i>>1]-pa[i]);
    ll t2=-(pb[i]-pb[i>>1]);

    #ifdef DEBUG
    cout<<t2<<"x = "<<t1<<" (mod "<<n-1<<")"<<endl;
    #endif // DEBUG

    t1%=(n-1); if(t1<0) t1+=n-1;
    t2%=(n-1); if(t2<0) t2+=n-1;

    ll g = gcd(t2,n-1),m;
    ll ans;

    if( t1!=0 && t2!=0 && (t1/g)*g == t1 ){
        t2/=g;
        t1/=g;
        m=(n-1)/g;

        #ifdef DEBUG
        cout<<t2<<"x = "<<t1<<" (mod "<<m<<")"<<endl;
        cout<<inverse(t2,m)<<endl;
        cout<<"x = "<<(t1*inverse(t2,m))%m<<" (mod "<<m<<")"<<endl;
        #endif // DEBUG
        ans = (t1*inverse(t2,m))%m; if(ans<0) ans+=m;
        for(; ans<n-1; ans+=m){
            #ifdef DEBUG
            cout<<a<<"^"<<ans<<" (mod "<<n<<" )= "<<pow(a,ans,n)<<"// b = "<<b<<endl;
            #endif // DEBUG
            if(b==pow(a,ans,n)) break;
        }

        ans%=n-1; if(ans<0) ans+=n-1;
        return ans;
    } else return GetLogarithmByPollard(a,n,b);
}

ll GetLogarithmIndexing(ll a, ll n, ll b){
    vector<ll> B,p;
    BuildSieveOfEratosthenes(10,B);
    int t = B.size();

    vector<vector<ll> > matrix;
    vector<int> m_r(t+1,-1), m_c;
    vector<ll> vect;

    ll k;
    while(matrix.size()<t){
        k=rand()%n;
        if(isSmooth(pow(a,k,n),B,vect)){
            vect.push_back(k);
            AddToGaussSystem(matrix,m_c,m_r,vect,n-1);
        }
    }

    for(bool smooth = false; !smooth; ){
        k=rand()%n;
        smooth=isSmooth((b*pow(a,k,n))%n,B,vect);
    }

    ll x = -k;
    for(int i=0; i<t; ++i){
        /*cout<<vect[i]<<endl;
        cout<<m_r[i]<<endl;
        cout<<matrix[m_r[i]][t]<<endl;*/
        x+=vect[i]*matrix[m_r[i]][t];
    }

    x%=n-1; if(x<0) x+=n-1;
    return x;
}

ll GetLogarithmByPohligHellman(ll g, ll h, ll p, ll s){
    ll n=1; for(ll i=0; i<s; ++i) n*=p; ++n; // n = p^s +1
    ll x=0;
    ll tp=(n-1)/p;     //    p^(s-1)
    ll gp=pow(g,tp,n); // g^(p^(s-1))
    ll tg=1;
    ll th=h;
    ll xi;

    vector<ll> revgs(p,1);
    revgs[p-1]=pow(g,n-p,n);
    for(int i=p-2; i>=0; --i) revgs[i]=(revgs[i+1]*g)%n;

    #ifdef DEBUG
    cout<<tp<<endl;
    cout<<gp<<endl;
    #endif // DEBUG

    unordered_map<ll,ll> xv;
    for(int i=0; i<p; ++i){
        #ifdef DEBUG
        cout<<"tg: "<<tg<<endl;
        #endif // DEBUG
        xv[tg]=i;
        tg*=gp, tg%=n;
    }

    #ifdef DEBUG
    for(auto& e: revgs) cout<<e<<" "; cout<<endl;
    for(auto& e: xv) cout<<"("<<e.first<<", "<<e.second<<") "; cout<<endl;
    #endif // DEBUG

    ll pw=1;
    for(int i=0; i<s; ++i){
        xi=xv[pow(th,tp,n)];
        x+=xi*pw;
        th*=pow(g,n-1-pw*xi,n), th%=n;
        #ifdef DEBUG
        cout<<"h: "<<pow(th,tp,n)<<endl;
        cout<<xi<<" "<<endl;
        cout<<"P: "<<pw*xi<<endl;
        cout<<"th: "<<th<<endl;
        #endif // DEBUG
        pw*=p;
        tp/=p;
    }

    x%=n-1; if(x<0) x+=n-1;
    return x;
    /*vector<ll> gp(s,g), hp(s,h);
    for(ll i=1; i<s; ++i){
        gp[]
    } p[i]=(p[i-1]*p)%n;

    tg=pow()*/

}
