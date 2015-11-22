#ifndef _HELP_FUNCTIONS_
#define _HELP_FUNCTIONS_

#include "config.h"

ll sqrt(ll);
ll gcd(ll,ll);
ll gcd(ll,ll,ll&,ll&);
ll inverse(ll,ll);
ll pow(ll,ll,ll);
ll Random(ll,ll);

bool isQuadraticResidue(ll,ll);
void BuildSieveOfEratosthenes(int, vector<ll>&);
void SolveGauss(vector<vector<int> >& matrix, vector<int>& m_c, vector<int>& m_r);
bool isSmooth(ll a,vector<ll>& B, vector<ll>& powers);


#endif // _HELP_FUNCTIONS_
