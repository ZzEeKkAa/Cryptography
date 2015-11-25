#ifndef _HELP_FUNCTIONS_
#define _HELP_FUNCTIONS_

#include "config.h"

ll sqrt(ll);
ll gcd(ll,ll);
ll gcd(ll,ll,ll&,ll&);
ll inverse(ll,ll);
ll pow(ll,ll,ll);
ll Random(ll,ll);

point GroupAdd(ll a, ll b, point P, point Q, ll n);
point GroupMul(ll a, ll b, ll k, point P, ll n);

bool isQuadraticResidue(ll,ll);
void BuildSieveOfEratosthenes(int, vector<ll>&);
void SolveGauss(vector<vector<int> >& matrix, vector<int>& m_c, vector<int>& m_r);
bool isSmooth(ll a,vector<ll>& B, vector<ll>& powers);
bool AddToGaussSystem(vector<vector<ll> >& matrix, vector<int>& m_c, vector<int>& m_r, vector<ll>const & vect, ll p);

#endif // _HELP_FUNCTIONS_
