
#include <iostream>
#include <fstream>

#include <polynomial_selection.hpp>



// {{{ is_reducible()
bool is_reducible(Polynomial &polynomial, Target &target)
{
   NTL::ZZ div, root, val, product;

   root = SqrRoot(abs(polynomial.f.rep[0])) + 1;
   for(div = 1; div < root; div++)
   {
      if(target.n % div == 0)
      {
         val = div;
         product = 1;
         product *= F(polynomial.f, val);
         val = 0 - div;
         product *= F(polynomial.f, val);
         val = target.n / div;
         product *= F(polynomial.f, val);
         val = 0 - (target.n / div);
         product *= F(polynomial.f, val);
         if(product == 0)
         {
            return 1;
         }
      }
   }
   return 0;
}
// }}}

// {{{ polynomial_goodness()
NTL::RR polynomial_goodness(Polynomial &polynomial, const char *primes)
{
   NTL::RR result;
   NTL::ZZ temp, pZ, jZ;
   int p, j, count;
   std::ifstream fin;

	temp = 0;
	for(int i=0; i<polynomial.d; i++)
	{
		temp += NTL::coeff(polynomial.f, i) * NTL::coeff(polynomial.f, i);
	}	

   temp = NTL::SqrRoot(temp);
   result = NTL::to_RR(temp) / NTL::to_RR(polynomial.m);


   fin.open(primes);
   fin >> p;

   count = 0;
   while(p < 100)
   {
      pZ = p;
      for(j=0; j<p; j++)
      {
         jZ = j;
         if(F(polynomial.f, jZ) % pZ == 0)
            count++;
      }
      fin >> p;
   }
   fin.close();


   pZ = count;
   jZ = 72;   // Max possible roots modulo primes up to 100
   result -= NTL::to_RR(pZ) / NTL::to_RR(jZ);


   return result;
}
// }}}

// {{{ F()
// Polynomial function
NTL::ZZ F(NTL::ZZX &poly, NTL::ZZ &x)
{
   NTL::ZZ temp;

   temp = power(x, NTL::deg(poly)+1);
   for(int i = 0; i < NTL::deg(poly)+1; i++)
      temp += NTL::coeff(poly, i) * NTL::power(x, i);

   return temp;
}
// }}}

// {{{ FR()
NTL::RR FR(const NTL::ZZX &poly, NTL::RR &x)
{
   NTL::RR temp;

   temp = power(x, NTL::deg(poly)+1);
   for(int i = 0; i < NTL::deg(poly)+1; i++)
      temp += NTL::to_RR(NTL::coeff(poly, i)) * NTL::to_RR(NTL::power(x, i));

   return temp;
}
// }}}

// {{{ Newton()
NTL::RR Newton(Polynomial &polynomial, NTL::RR start)
{
	//TODO: degree 3 specific
   NTL::ZZ val;
   NTL::RR error = NTL::to_RR(0.001);
	NTL::RR prev;
	NTL::RR next;
	
	prev = 0;
	next = start;
	
	while(NTL::abs(next-prev) > error)
	{
		prev = next;

      NTL::RR val = FR(polynomial.f, prev) / 
			(NTL::to_RR(3) * NTL::to_RR(NTL::power(prev,2)) + 
          NTL::to_RR(2) * NTL::to_RR(polynomial.f.rep[2]) * prev + 
			 NTL::to_RR(polynomial.f.rep[1]));

      /*
      NTL::RR val2 = FR(polynomial.f, prev) / 
         (NTL::to_RR(3)*NTL::to_RR(NTL::power(prev,2) + 
          FR(NTL::diff(polynomial.f), prev)));

      std::cout << polynomial.f << std::endl;
      std::cout << "3x^2 + 2·" << polynomial.f.rep[2] << "x + " 
                << polynomial.f.rep[1] << " val=" << val << std::endl;
      std::cout << NTL::diff(polynomial.f) << "val=" << val2 << std::endl << std::endl;
      */

		next = prev - val;
   }

	return next;
}
// }}}

// {{{ has_roots()
bool has_roots(Polynomial &polynomial)
{
	// TODO: dregree 3 especific

   NTL::RR root = NTL::to_RR(-12.345);

   root = Newton(polynomial, root);

   if(NTL::power(NTL::MakeRR(polynomial.f.rep[2], 0) + root, 2) > 
		(-4 * NTL::MakeRR(polynomial.f.rep[0], 0) / root))
      return true;
   else
      return false;

}
// }}}

// {{{ get_base_m_expansion()
NTL::ZZX get_base_m_expansion(int degree, NTL::ZZ &m, NTL::ZZ &n)
{
   NTL::ZZX f;
   NTL::ZZ N;
   
   f.rep.SetLength(degree);

   N = n-NTL::power(m, degree);
   for(int i=degree-1; i>0; i--) 
   {
      f.rep[i] = N / NTL::power(m, i);
      N -= NTL::power(m, i) * f.rep[i];;
   }


   f.rep[0] = N;
   for(int i=0; i<degree-1; i++) 
   {
      if(f.rep[i] > m/2) 
      {
         f.rep[i] -= m;
         f.rep[i+1] ++;
      }
   }

   return f;
}
// }}}

// {{{ polynomial_selection()
void polynomial_selection(Polynomial &polynomial, Target &target, 
	const char *primes)
{
   int i;
   NTL::ZZ numZ;
   NTL::ZZ goodM;
   NTL::RR numberS;

   // Set m, root of polynomial modulo n
   bool stop = false;
   numZ = polynomial.d;
   polynomial.m = NTL::to_ZZ(NTL::pow(NTL::to_RR(target.n),1/NTL::to_RR(numZ)));

   numberS = 10;
   goodM = polynomial.m;

   
   // test m at random points
   NTL::ZZ width, start;
   start = polynomial.m;
   width = (int)std::pow((double)10, (double)target.digits / 6);
   start -= width;

   for(i=0; i<target.digits * target.digits ; i++)
   {
      polynomial.m = NTL::RandomBnd(width) + start;
      polynomial.f = get_base_m_expansion(polynomial.d, polynomial.m, target.n);

      if(has_roots(polynomial))
      {
			NTL::RR goodness = polynomial_goodness(polynomial, primes);
   		if(goodness < numberS)
         {
            numberS = goodness;
            goodM = polynomial.m;
         }
		}

   }
   polynomial.m = goodM;
   polynomial.f = get_base_m_expansion(polynomial.d, polynomial.m, target.n);

   int loopTimes;
   if(target.digits > 28)
      loopTimes = (int)(pow((double)10, (double)3));
   
   else
      loopTimes = (int)(std::pow((double)10, (double)target.digits / 10)) + 10;
   
   polynomial.m = FloorToZZ(pow(MakeRR(target.n, 0), 1 / MakeRR(numZ, 0)));
   while(!stop && loopTimes-- > 0)
   {
      polynomial.f = get_base_m_expansion(polynomial.d, polynomial.m, target.n);

      if(has_roots(polynomial))
      {
			NTL::RR goodness = polynomial_goodness(polynomial, primes);
   		if(goodness < numberS)
         {
            numberS = goodness;
            goodM = polynomial.m;
         }
		}

      if(2 * polynomial.m * polynomial.m * polynomial.m < target.n)
         stop = true;
      
      polynomial.m--;
   }
   polynomial.m = goodM;
   polynomial.f = get_base_m_expansion(polynomial.d, polynomial.m, target.n);



   // f(x) reducible?
   if(is_reducible(polynomial, target))
   {
      std::cout << "\tf(x) is reducible: " << polynomial.f << std::endl;
      exit(1);
   }
   

}
// }}}


