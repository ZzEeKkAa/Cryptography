#ifndef _FACTORIZATION_
#define _FACTORIZATION_

#include "config.h"

ll GetFactorFerma(ll);
ll GetFactorQuadraticSieve(ll);
ll GetFactorByDixon(ll);
ll GetFactorByPollard(ll);
ll GetFactorByLenstra(ll);
void GetAllFactors(ll n, ll(*GetFactor)(ll), vector<ll>&factors);

#endif // _FACTORIZATION_
