
#include <sieve.hpp>
#include <linear_algebra.hpp>

#include <NTL/GF2.h>
#include <NTL/vec_GF2.h>
#include <NTL/mat_GF2.h>

// {{{ block_lanczos()
bool
block_lanczos(const NTL::mat_GF2 & A, NTL::vec_GF2 & x, const NTL::vec_GF2 & y)
{
   using namespace NTL;

   long nrows = A.NumRows(); 
   long ncols = A.NumCols();

   std::cout << "solving " << nrows << "x" << ncols << " matrix" << std::endl;

   if (nrows<ncols) 
   {
      std::cout << "error: nrows("<<nrows<<") < ncols("<<ncols<<")"<<std::endl;
      return false;
   }

   /* temporaries */
   GF2 t1;
   vec_GF2 ts1;
   vec_GF2 tv1;

   /* D is a (nrows x nrows) diagonal matrix with random entries (squared) */
   vec_GF2 D;

   D.SetLength(nrows);
   for(long i = 0; i < nrows; ++i)
   {
      do
      {
         D[i] = random_GF2();
      }
      while(IsZero(D[i]));
      //sqr(D[i], D[i]);
   }

   /* w0 = w = transpose(A)*D*y */
   vec_GF2 w;

   w.SetLength(ncols);

   for(long j = 0; j < nrows; ++j)
   {
      mul(t1, D[j], y[j]);
      mul(ts1, A[j], t1);
      add(w, w, ts1);
   }

   vec_GF2 w0(w);

   /* v1 = transpose(A)*D*A*w0 */
   vec_GF2 v1;

   v1.SetLength(ncols);
   for(long j = 0; j < nrows; ++j)
   {
      InnerProduct(t1, A[j], w0);
      mul(t1, D[j], t1);
      mul(ts1, A[j], t1);
      add(v1, v1, ts1);
   }

   /* dv = |v1|^2 */
   GF2 dv;

   InnerProduct(dv, v1, v1);

   /* dvw = (w0,v1) */
   GF2 dvw;

   InnerProduct(dvw, w0, v1);

   /* check if dvw is zero */
   if(IsZero(dvw))
   {
      std::cerr << "Lanczos() (w0,v1)=0\n";
      //return false;
      dvw=1;
   }

   /* w1 = v1 - (dv/dvw)*w0 */
   vec_GF2 w1;

   div(t1, dv, dvw);
   mul(tv1, t1, w0);
   sub(w1, v1, tv1);

   /* b = (w0,w) */
   GF2 b;

   InnerProduct(b, w0, w);

   /* initialize solution */
   div(t1, b, dvw);
   mul(x, t1, w0);

   /* variables used in loop */
   vec_GF2 v2;
   vec_GF2 w2;

   v2.SetLength(ncols);

   /* status display */
   long count = 0;
   long last_count = 0;
   long last_percent = 0;

   do
   {
      /* v2 = transpose(A)*D*A*w1 */
      clear(v2);
      for(long j = 0; j < nrows; ++j)
      {
         InnerProduct(t1, A[j], w1);
         mul(t1, D[j], t1);
         mul(ts1, A[j], t1);
         add(v2, v2, ts1);
      }

      long percent = (++count) * 100 / nrows;

      if((percent != last_percent) && (count - last_count >= 20))
      {
         std::
            cout << "Lanczos: " << count << " itterations (" << percent <<
            "%)  \r";
         std::cout.flush();
         last_percent = percent;
         last_count = count;
      }

      /* dvw = (w1,v2) */
      InnerProduct(dvw, w1, v2);
      if(IsZero(dvw))
      {
         std::cout << "Lanczos: " << count << " itterations        \n";

         if(IsZero(w1))
            break;

         std::cerr << "Lanczos() algorithm failed!\n";
         return false;
      }

      /* b = (w1,w) */
      InnerProduct(b, w1, w);

      /* add to solution */
      div(t1, b, dvw);
      mul(tv1, t1, w1);
      add(x, x, tv1);

      /* dv = |v2|^2 */
      InnerProduct(dv, v2, v2);

      /* dv0 = (v2,v1) */
      GF2 dv0;

      InnerProduct(dv0, v2, v1);

      /* dvw0 = (w0,v1) */
      GF2 dvw0;

      InnerProduct(dvw0, w0, v1);

      /* w2 = v2 - (dv/dvw)*w1 - (dv0/dvw0)*w0 */
      w2 = v2;
      div(t1, dv, dvw);
      mul(tv1, t1, w1);
      sub(w2, w2, tv1);
      div(t1, dv0, dvw0);
      mul(tv1, t1, w0);
      sub(w2, w2, tv1);

      /* shift variables */
      w0 = w1;
      w1 = w2;
      v1 = v2;

   }
   while(true);

   return true;
}
// }}}

// {{{ Legendre()
//1 if y is a quadratic residue modulo z, -1 otherwise
int Legendre(NTL::ZZ & y, NTL::ZZ & z)
{
   NTL::ZZ x, yP;

   yP = y % z;
   for(x = 0; x < z; x++)
   {
      if(((x * x) - yP) % z == 0)
      {
         return 1;
      }
   }
   return -1;
}
// }}}

// {{{ linear_algebra()
void linear_algebra(Polynomial &polynomial, Target &target, 
	FactorBase &fb, Matrix &matrix, 
    const std::vector<int> &av, const std::vector<int> &bv)
{

   NTL::ZZ aZ;
   NTL::ZZ bZ;
   NTL::ZZ pZ;
   NTL::ZZ valZ;
   NTL::ZZ numZ;
   int i;
   int j;
   int k;
   int u = polynomial.d*target.t;

   // Initialize sM
   for(j = 0; j <= matrix.row-1; j++)
   {
      // Initialize row
      for(k = 0; k < matrix.col-1; k++)
         matrix.sM[k][j] = 0;
      
      // Set the first column
      aZ = av[j];
      bZ = bv[j];
      valZ = aZ + bZ * polynomial.m;
      if(valZ < 0)
      {
         valZ *= -1;
         matrix.sM[0][j] = 1;
      }

      // Set a RFB row 
      i = 0;
      while(i < target.t && valZ != 1)
      {
         pZ = fb.RFB[i];
         if(valZ % pZ == 0)
         {
            if(matrix.sM[1+i][j]==0) matrix.sM[1+i][j]=1;
            else matrix.sM[1+i][j]=0;
            valZ = valZ / pZ;
         } 
         else i++;
      }

      // Set a AFB row 
      valZ = algebraic_norm(polynomial, av[j], bv[j]);
      if(valZ < 0)
         valZ *= -1;
      
      i = 0;
      while(i<u && valZ!=1)
      {
         pZ = fb.AFB[i];
         if(valZ % pZ == 0)
         {
            numZ = fb.AFBr[i];
            while((aZ + bZ * numZ) % pZ != 0)
            {
               pZ = fb.AFB[++i];
               numZ = fb.AFBr[i];
            }

            if(matrix.sM[1+target.t+i][j]==0) matrix.sM[1+target.t+i][j]=1;
            else matrix.sM[1+target.t+i][j]=0;

            valZ = valZ / pZ;
         }
         else i++;
         
      }

      // Set a QCB row
      for(i=0; i<target.digits; i++)
      {
         numZ = fb.QCB[i];
         valZ = fb.QCBs[i];
         valZ = aZ + bZ * valZ;
         if(Legendre(valZ, numZ) != 1)
         {
            matrix.sM[1+target.t+u+i][j]=1;
         }
      }
   }

   /*
   NTL::vec_GF2 x;
   NTL::vec_GF2 y;
   x.SetLength(matrix.sM.NumRows());
   y.SetLength(matrix.sM.NumRows());

   while(1)
   if(!block_lanczos(matrix.sM, x, y))
   {
      std::cout << "block_lanczos() error" << std::endl;
   }
   std::cout << matrix.sM << std::endl;
   exit(0);
   */

   bool stop;

   // Solve Mx = 0
   NTL::GF2 temp;

   for(i=0; i < matrix.col-1; i++)
   {
      matrix.sfreeCols[i] = 0;
      stop = false;
      j = i;
      if(matrix.sM[i][j]==0)
      {
         // find k so that M[k,j]==1
         k = i;
         while(++k < matrix.col-1 && matrix.sM[k][j]==0)
         {
         }
         // swap row i with row k
         if(k >= matrix.col-1)
         {
            matrix.sfreeCols[i] = 1;
            stop = true;
         }
         for(j = 0; j <= matrix.row-1 && !stop; j++)
         {
            temp = matrix.sM[i][j];
            matrix.sM[i][j] = matrix.sM[k][j];
            matrix.sM[k][j] = temp;
         }

      }
      for(k = 0; k < matrix.col-1 && !stop; k++)
      {
         if(matrix.sM[k][i]==1 && k != i)
         {
            for(j = 0; j <= matrix.row-1; j++)
            {
               matrix.sM[k][j] =  
               ((NTL::IsOne(matrix.sM[k][j]) ||NTL::IsOne(matrix.sM[i][j])) &&
               ((NTL::IsZero(matrix.sM[k][j])||NTL::IsZero(matrix.sM[i][j]))));
            }
         }
      }
   }

}
// }}}



