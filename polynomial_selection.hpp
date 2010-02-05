
#ifndef __GGNFS_POLYNOMIAL_SELECTION_HPP__
#define __GGNFS_POLYNOMIAL_SELECTION_HPP__


#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/ZZX.h>

#include <gnfs.hpp>

class Polynomial
{
	public:
   NTL::ZZX f;
	std::vector<NTL::RR> roots;
	NTL::ZZ m;
	int d;
};


void find_roots(Polynomial &polynomial);

bool is_reducible(Polynomial &polynomial, Target &target);

NTL::RR polynomial_goodness(Polynomial &polynomial, const char *primes);

NTL::ZZ F(NTL::ZZX &poly, NTL::ZZ &x);

NTL::RR FR(const NTL::ZZX &poly, NTL::RR &x);

NTL::ZZ dF(NTL::ZZX &poly, NTL::ZZ &x);

NTL::RR Newton(Polynomial &polynomial, NTL::RR start);

bool has_roots(Polynomial &polynomial);

NTL::ZZX get_base_m_expansion(int degree, NTL::ZZ &m, NTL::ZZ &n);

void polynomial_selection(Polynomial &polynomial, Target &target, 
	const char *primes);


#endif

