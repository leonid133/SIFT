//---------------------------------------------------------------------------

#pragma hdrstop

#include "matrc.h"

//---------------------------------------------------------------------------
//class Vektor****************************************************************
//----------------------------------------------------------------------------
const float ToGrad = 57.2957795130823;

//Перегрузка операции +  ---------------------------------------------------
vektor operator + (const vektor& u,const vektor& v)
{
	return vektor (u.z+v.z, u.x+v.x);
}

//Перегрузка операции -  ---------------------------------------------------
vektor operator - (const vektor& u,const vektor& v)
{
	return vektor (u.z-v.z, u.x-v.x);
}

//Перегрузка операции *. Умножение скаляра на вектор.------------------------
vektor operator* (const vektor& v, float f)
{
	return vektor (f*v.z, f*v.x);
}

//Перегрузка операции *. Умножение скаляра на вектор.------------------------
vektor operator* (float f, const vektor& v)
{
	return vektor (f*v.z, f*v.x);
}

//Деление вектора на скаляр.-------------------------------------------------
vektor operator / (const vektor& v, float f)
{
	return vektor (v.z/f, v.x/f);
}

//Косинус угла между векторами.---------------------------------------------
float cos(const vektor& a, const vektor& b)
{
   float ab = !a*!b;
   if(ab > 0.0000001)
   {
      float k2 = (a*b)/ab;
      if (k2 > 1) k2 = 1;
      if (k2 < -1) k2 = -1;
   	return k2;
   }
   else
      return 0;
}

//----------------------------------------------------------------------------
String vektor:: show(void)
{
  String msg;  msg.sprintf("%6.3f %6.3f ", z, x);  return msg;
}

//атангенс угла между векторами.----------------------------------------------
float atan(const vektor& a, const vektor& b)
{
   return a.alfa() - b.alfa();
}

//class Matric ***************************************************************
//----------------------------------------------------------------------------
template <class TYPE_MATRIC>
matric<TYPE_MATRIC>::matric(const long m_, const long n_)
   : m(m_)
   , n(n_)
{
   if (m == 0 && n == 0)
      ar = NULL;
   else
   {
      ar = new TYPE_MATRIC [m*n];
      for(long i = 0; i < (m*n); i++)
         ar[i] = 0;
   }
};

//----------------------------------------------------------------------------
template <class TYPE_MATRIC> matric<TYPE_MATRIC>::matric(const matric &r)
{
   m = r.m;
   n = r.n;
   ar = new TYPE_MATRIC [m*n];

   for (long j = 0; j < m*n; j++)
      ar[j] = r.ar[j];
};

//-----------------------------------------------------------------------------
template <class TYPE_MATRIC>
const matric<TYPE_MATRIC> &matric<TYPE_MATRIC>::operator=(const matric &r)
{
   if(this != &r)
   {
      delete [] ar;
      m = r.m;
      n = r.n;
      ar = new TYPE_MATRIC [m*n];

      for (long j = 0; j < m*n; j++)
         ar[j] =  r.ar[j];
   }
   return *this;
};

//Сложение--------------------------------------------------------------------
/*virtual*/template <class TYPE_MATRIC>
matric<TYPE_MATRIC> matric<TYPE_MATRIC>::operator+(const matric &r)
{
   matric temp(m, n);

   for (long j = 0; j < m*n; j++)
      temp.ar[j] =  ar[j]+r.ar[j];
   return temp;
};

//Вычитание ------------------------------------------------------------------
/*virtual*/template <class TYPE_MATRIC>
matric<TYPE_MATRIC> matric<TYPE_MATRIC>::operator-(const matric &r)
{
   matric temp(m, n);
   for (long j = 0; j < m*n; j++)
      temp.ar[j] =  ar[j]-r.ar[j];

   return temp;
};

//Умножение-------------------------------------------------------------------
/*virtual*/template <class TYPE_MATRIC>
matric<TYPE_MATRIC> matric<TYPE_MATRIC>::operator*(const matric &r)
{
   matric temp(m, n);
   for (long j = 0; j < m; j++)
      for (long k = 0; k < n; k++)
      {
         temp.ar[j*n+k] = 0;
         for (long i = 0; i < n; i++)
            temp.ar[j*n+k] += ar[j*n+i] * r.ar[i*n+k];
      }
   return temp;
};

//Умножение на скаляр----------------------------------------------------------
/*virtual*/template <class TYPE_MATRIC>
matric<TYPE_MATRIC> matric<TYPE_MATRIC>::operator*(const TYPE_MATRIC x)
{
   matric temp(m, n);
   for (int j = 0; j < m*n; j++)
      temp.ar[j] =  ar[j]*x;

   return temp;
};

//Определитель----------------------------------------------------------------
/*virtual*/template <class TYPE_MATRIC>
TYPE_MATRIC  matric<TYPE_MATRIC>::D(void)
{
   TYPE_MATRIC temp = 1, temp2 = -1, sum = 0;
   long i, i2;
   for(long kk = 1; kk < n+1; kk++)
   {
      i = kk - 1 ;
      i2 = n - kk;
      for (long j = 0; j < m; j++)
      {
         if (i >= n){i = 0;}
         temp = temp * ar[j*n+i++];
         if (i2 < 0){i2 = n - 1;}
         temp2 = temp2 * ar[j*n+i2--];
      }
      sum += temp + temp2;
      temp = 1;
      temp2 = -1;
   }
   return sum;
};

//Транспонирование матриц-----------------------------------------------------
/*virtual*/template <class TYPE_MATRIC> matric<TYPE_MATRIC> matric<TYPE_MATRIC>::T(void)
{
   matric temp(n, m);
   for(long j = 0; j < n; j++)
      for(long k = 0; k < m; k++)
         temp.ar[j*n+k] = ar[k*m+j];

   return temp;
};

//След матрицы ---------------------------------------------------------------
/*virtual*/template <class TYPE_MATRIC> TYPE_MATRIC matric<TYPE_MATRIC>::Sp(void)
{
    TYPE_MATRIC Spur=0;
    long x=m<n?m:n;
    for(long i=0; i<x; i++)
    Spur+=ar[i*m+i];
    return Spur;
};
//Вывод-----------------------------------------------------------------------
String matric<float>::Show(void)
{
   String tmp(""), txt;
   for (long j = 0; j < m; j++)
      for (long k = 0; k < n; k++)
      {
         tmp = tmp.FormatFloat("0.00     ", ar[j*n+k]);
         txt += tmp;
      }
   return txt;
};
//Вывод-----------------------------------------------------------------------
String matric<int>::Show(void)
{
   String tmp(""), txt;
   for (long j = 0; j < m; j++)
      for (long k = 0; k < n; k++)
      {
         tmp = ar[j*n+k];
         txt += tmp;
      }
   return txt;
};
//Фильтр Гаусса---------------------------------------------------------------
matric<float> FilterGauss (float sigma, int n, float factorY) // sigma(0.5..1.5) n(1..3) factorY(1..2.5)
{
   matric<float> FilterGauss(n, n);
   for (int j = 0; j < n; j++)
     for (int k = 0; k < n; k++)
        FilterGauss.ar[j*n+k] = factorY * (exp((-pow(k, 2)-pow(j, 2))/(2*pow(sigma, 2)))/(2 * M_PI * pow(sigma, 2)));
   return FilterGauss;
};

//Class Image*****************************************************************
//----------------------------------------------------------------------------
const image &image::operator=(const image &img)
{
   if(this != &img)
   {
      m = img.m;
      n = img.n;
      dm = img.dm;
      dn = img.dn;

      R = img.R;
      G = img.G;
      B = img.B;
      Y = img.Y;
      Sx = img.Sx;
      Sy = img.Sy;
      Sxy = img.Sxy;
      Lnorm = img.Lnorm;
      H = img.H;
   }
   return *this;
};

//Инициализация из потока-----------------------------------------------------
/*virtual*/image image::LoadFromFile(ifstream &file)
{
   dm=dn=0;
   unsigned char buf[1023];
   file.read(buf, 54);
   long int location = ((long int)(buf[10]))+((long int)(buf[11]) << 8)+((long int)(buf[12]) << 16)+((long int) (buf[13]) << 24);
   int n = ((long int)(buf[18])) + ((long int)(buf[19]) << 8) + ((long int)(buf[20]) << 16) + ((long int)(buf[21]) << 24);
   int m = ((long int)(buf[22]))+((long int)(buf[23]) << 8)+((long int)(buf[24]) << 16)+((long int) (buf[25]) << 24);
   long int bytes_bitmap_data = ((long int)(buf[34]))+((long int)(buf[35]) << 8)+((long int)(buf[36]) << 16)+((long int) (buf[37]) << 24);
   file.ignore((location - 54));
   long int j = (m - 1), k = 0, i = 0, mem_i = 0, s;
   B = G = R = Y = Lnorm = Sx = Sy = Sxy = H = matric<float>(m, n);
   while (mem_i <= bytes_bitmap_data)
   {
      file.read(buf, 3);
      while ((j >= 0) && (i <= 2))
      {
         B.ar[j*B.n+k] = buf[i++];
         G.ar[j*G.n+k] = buf[i++];
         R.ar[j*R.n+k++] = buf[i++];
         if (k >= n)
         {
            k = 0;
            j--;
         }
      }
      mem_i+=i;
      i = 0;
      if (file.eof())
         break;
   }
   file.close();
   Y =(B*0.11 + G*0.59 + R*0.3);
   return *this;
};

//Свертка---------------------------------------------------------------------
/*virtual*/image image :: Convolution(const matric<float> &Filter)
{
   matric<float> A(Filter.m, Filter.n);
   matric<float> temp(Y);

   for (int jm = Filter.m; jm < Y.m; jm++)
   {
      for (int kn = Filter.n; kn < Y.n; kn++)
      {
         for (int j = 0; j < Filter.m; j++)
            for (int k = 0; k < Filter.n; k++)
               A.ar[j*Filter.n+k] = temp.ar[(jm - Filter.m + j)*Y.n+(kn - Filter.n + k)];
//Умножение A = A * Filter и свертка в temp
         for (int jj = 0; jj < A.m; jj++)
            for (int kk = 0; kk < A.n; kk++)
            {
               int l = (jm - Filter.m + jj)*Y.n + (kn - Filter.n + kk);
               temp.ar[l] = 0;
               for (int i = 0; i < A.n; i++)
                  temp.ar[l] += A.ar[jj*A.n+i] * Filter.ar[i*Filter.n+kk];
            }
//-------------------------------------
      }
   }
   Y = temp;
   return *this;
};

//Эквализация-----------------------------------------------------------------
/*virtual*/ image  image :: Eqalization(void)
{
   float Ymax=0, Ymin=255, Rmax=0, Rmin=255, Gmax=0, Gmin=255, Bmax=0, Bmin=255;
   m = Y.m;
   n = Y.n;
   for (int j = 0; j < m; j++)
      for (int k = 0; k < n; k++)
      {
         if(Y.ar[j*n+k]>Ymax){Ymax=Y.ar[j*n+k];}
         if(Y.ar[j*n+k]<Ymin){Ymin=Y.ar[j*n+k];}
         if(R.ar[j*n+k]>Rmax){Rmax=R.ar[j*n+k];}
         if(R.ar[j*n+k]<Rmin){Rmin=R.ar[j*n+k];}
         if(G.ar[j*n+k]>Gmax){Gmax=G.ar[j*n+k];}
         if(G.ar[j*n+k]<Gmin){Gmin=G.ar[j*n+k];}
         if(B.ar[j*n+k]>Bmax){Bmax=B.ar[j*n+k];}
         if(B.ar[j*n+k]<Bmin){Bmin=B.ar[j*n+k];}
      }
   for (int j = 0; j < m; j++)
      for (int k = 0; k < n; k++)
      {
         Y.ar[j*n+k] = 255*(Y.ar[j*n+k]-Ymin)/(Ymax-Ymin);
         R.ar[j*n+k] = 255*(R.ar[j*n+k]-Rmin)/(Rmax-Rmin);
         G.ar[j*n+k] = 255*(G.ar[j*n+k]-Gmin)/(Gmax-Gmin);
         B.ar[j*n+k] = 255*(B.ar[j*n+k]-Bmin)/(Bmax-Bmin);
      }
   return *this;
}
//DiscretCosinusTransformation -----------------------------------------------
/*virtual*/ image image :: DiscretCosinusTransformation(void)
{
   int NN = 8;
   matric<float> Ya(NN, NN), Ra(NN, NN), Ga(NN, NN), Ba(NN, NN);
   for (int m = 0; m < (R.m); m+=NN)
      for (int n = 0; n < (R.n); n+=NN)
      {
         //Ya = DCT(Y, m, n);
         Ra = DCT(R, m, n, NN);
         Ga = DCT(G, m, n, NN);
         Ba = DCT(B, m, n, NN);
         for(int j=0; j<NN; j++)
            for(int k=0; k<NN; k++)
            {
               if((m+j)<R.m && (n+k)<R.n)
               {
                  R.ar[(m+j)*R.n + (n+k)]=Ra.ar[j*NN+k];
                  G.ar[(m+j)*G.n + (n+k)]=Ga.ar[j*NN+k];
                  B.ar[(m+j)*B.n + (n+k)]=Ba.ar[j*NN+k];
                 // Y.ar[(m+j)*Y.n + (n+k)]=Ya.ar[j*NN+k];
               }
            }
      }
   return *this;
}
//InverseDiscretCosinusTransformation -----------------------------------------------
/*virtual*/ image image :: InverseDiscretCosinusTransformation(void)
{
   int NN=8;
   matric<float> Ya(NN, NN), Ra(NN, NN), Ga(NN, NN), Ba(NN, NN);

   for (int m = 0; m < (R.m); m+=NN)
      for (int n = 0; n < (R.n); n+=NN)
      {
         //Ya = IDCT(Y, m, n);
         Ra = IDCT(R, m, n, NN);
         Ga = IDCT(G, m, n, NN);
         Ba = IDCT(B, m, n, NN);
         for(int j=0; j<NN; j++)
            for(int k=0; k<NN; k++)
            {
               if((m+j)<R.m && (n+k)<R.n)
               {
                  R.ar[(m+j)*R.n + (n+k)]=Ra.ar[j*NN+k];
                  G.ar[(m+j)*G.n + (n+k)]=Ga.ar[j*NN+k];
                  B.ar[(m+j)*B.n + (n+k)]=Ba.ar[j*NN+k];
                  //Y.ar[(m+j)*Y.n + (n+k)]=Ya.ar[j*NN+k];
               }
            }
      }
   return *this;
}
//Fast Furye Transformation --------------------------------------------------
/*virtual*/ image  image :: FastFuryeTransformation(void)
{
   matric<int> A(8, 8);
   matric<float> temp(Y);

   for (int m = 0; m < Y.m; m+=8)
      for (int n = 0; n < Y.n; n+=8)
      {
         A = FFT(Y, m, n);
         for(int j=0; j<8; j++)
            for(int k=0; k<8; k++)
            {
               Y.ar[(m+j)*Y.n + (n+k)]=A.ar[j*8+k];
            }
      }
   return *this;
}
//---------------------------------------------------------------------
/*virtual*/ image image :: HarrissLaplass(void)
{
   float Etta = 1.2, S = 0.7, h;
   matric<float> A(2, 2), GaussN(2, 2), GaussD(2, 2);
   int al, bl, cl, dl, el, fl, gl, hl, il, l, lx, ly, lxy;
   int Ym = Y.m; int Yn = Y.n;
   matric<float> L(Ym, Yn), Lx(Ym, Yn), Ly(Ym, Yn), Lxx(Ym, Yn), Lxy(Ym, Yn), Lxxyy(Ym, Yn), Lyy(Ym, Yn), LoG(Ym, Yn);

   for(int scale = 1; scale < 10; scale++)
   {
      GaussN = FilterGauss(pow(Etta, scale), 2, 1);
      GaussD = FilterGauss(S*pow(Etta, scale), 2, 1);

      for (int jm = GaussN.m+2; jm < Y.m; jm++)
      {
         for (int kn = GaussN.n+2; kn < Y.n; kn++)
         {
            for (int j = 0; j < GaussN.m; j++)
               for (int k = 0; k < GaussN.n; k++)
                  A.ar[j*GaussN.n+k] = Y.ar[(jm - GaussN.m + j)*Y.n+(kn - GaussN.n + k)];
            for (int jj = 0; jj < A.m; jj++)
               for (int kk = 0; kk < A.n; kk++)
               {
                  cl = (jm - GaussN.m + jj - 2)*Y.n + (kn - GaussN.n + kk);
                  al = cl-2;   // al = (jm - GaussN.m + jj - 2)*Y.n + (kn - GaussN.n + kk - 2);
                  bl = cl-1;   // bl = (jm - GaussN.m + jj - 2)*Y.n + (kn - GaussN.n + kk - 1);
                  dl = (jm - GaussN.m + jj - 1)*Y.n + (kn - GaussN.n + kk - 2);
                  lxy = el = (jm - GaussN.m + jj - 1)*Y.n + (kn - GaussN.n + kk - 1);
                  ly = fl = (jm - GaussN.m + jj - 1)*Y.n + (kn - GaussN.n + kk);
                  gl = (jm - GaussN.m + jj)*L.n + (kn - GaussN.n + kk - 2);
                  lx = hl = (jm - GaussN.m + jj)*Y.n + (kn - GaussN.n + kk - 1);
                  l = il = (jm - GaussN.m + jj)*Y.n + (kn - GaussN.n + kk);

                  L.ar[l] = 0;
                  //Нормальный Лаплассиан Гауссиана
                  for (int i = 0; i < A.n; i++)
                     L.ar[l] += A.ar[jj*A.n+i] * (GaussD.ar[i*GaussN.n+kk] - GaussN.ar[i*GaussN.n+kk]);
                  Lnorm.ar[l] += L.ar[l] * pow(S*pow(Etta, scale), 2);
                 //Встроенный градиент
                 /* Lx.ar[l] = (Lnorm.ar[l]-Lnorm.ar[lx]);
                  Ly.ar[l] = (Lnorm.ar[l]-Lnorm.ar[ly]);
                  Lxx.ar[l] = Lx.ar[l]-Lx.ar[lx];
                  Lyy.ar[l] = Lx.ar[l]-Lx.ar[ly];
                  Lxy.ar[l] = (Lnorm.ar[l] - Lnorm.ar[lxy]);
                  Lxxyy.ar[l] = Lxy.ar[l]-Lxy.ar[lxy];
                  if (scale == 1)
                  {
                     Sx.ar[l] = (Lnorm.ar[cl] + (2 * Lnorm.ar[fl]) + Lnorm.ar[il]) - (Lnorm.ar[al] + (2 * Lnorm.ar[dl]) + Lnorm.ar[gl]);
                     Sy.ar[l] = (Lnorm.ar[gl] + (2 * Lnorm.ar[hl]) + Lnorm.ar[il]) - (Lnorm.ar[al] + (2 * Lnorm.ar[bl]) + Lnorm.ar[cl]);
                     Sxy.ar[l] = pow(pow(Sx.ar[l], 2)+ pow(Sy.ar[l], 2), 0.5);
                  }      */
                  // Фактор Харрисона
                 /* if(pow(pow(Lxx.ar[l]+Lyy.ar[l], 2), 0.5)>LoG.ar[l] && pow(pow(Lxx.ar[l]+Lyy.ar[l], 2), 0.5)>10)
                  {
                     h = pow(Etta, scale)*(Lxx.ar[l]*Lyy.ar[l]) - Lxxyy.ar[l]*Lxxyy.ar[l] +
                     0.04*pow(pow(Etta, scale)*(Lxx.ar[l]+Lyy.ar[l]), 2);

                     if(h > H.ar[l] && h>250)
                     {H.ar[l] = h/10;}
                  }
                  LoG.ar[l] = pow(pow(Lxx.ar[l]+Lyy.ar[l], 2), 0.5);   */
               }
         }
      }
   }
   return *this;
}


//Тестовый вывод--------------------------------------------------------------
void image::Draw(void)
{
   m = R.m;
   n = R.n;
   for (int j = 0; j < m; j++)
      for (int k = 0; k < n; k++)
        Form1->Image1->Canvas->Pixels[k][j] = (int)R.ar[j*n+k] + ((int)(G.ar[j*n+k]) << 8) + ((int)(B.ar[j*n+k]) << 16);
};

void image::DrawY(void)
{
   m = Y.m;
   n = Y.n;
   for (int j = 0; j < m; j++)
      for (int k = 0; k < n; k++)
        Form1->Image1->Canvas->Pixels[k][j] = (int)Y.ar[j*n+k] + ((int)(Y.ar[j*n+k]) << 8) + ((int)(Y.ar[j*n+k]) << 16);
};

// Discret Cosinus Transformation*--------------------------------------------
matric<float> DCT(matric<float> &r, int M, int N, int NN)
{
   matric<float> DCTkvant(NN, NN);
   matric<float> kvant(NN, NN);
   float PI = 3.1415926535897932384626433832795;
   float Cu, Cv, AC;
   for (int u = 0; u < NN; u++)
   {
      for (int v = 0; v < NN; v++)
      {
         if(u != 0)
         {Cu=1;}
         else
         {Cu=1/(pow(2, 0.5));}
         if(v != 0)
         {Cv=1;}
         else
         {Cv=1/(pow(2, 0.5));}
         AC=0;
         for (int x = 0; x < NN; x++)
         {
            for (int y = 0; y < NN; y++)
            {
               if((y+M)<r.m && (x+N)<r.n) {AC += r.ar[(y+M)*r.n+x+N]*cos((2*(x)+1)*(u)*PI/(2*NN))*cos((2*(y)+1)*(v)*PI/(2*NN));}
            }
         }
         DCTkvant.ar[v*NN+u]=0.25*Cu*Cv*AC;
      }
   }
 /*  DCTkvant.ar[0]=DCTkvant.ar[0]/8;
   DCTkvant.ar[1]=DCTkvant.ar[1]/11;
   DCTkvant.ar[2]=DCTkvant.ar[2]/10;
   DCTkvant.ar[3]=DCTkvant.ar[3]/16;
   DCTkvant.ar[4]=DCTkvant.ar[4]/24;
   DCTkvant.ar[5]=DCTkvant.ar[5]/40;
   DCTkvant.ar[6]=DCTkvant.ar[6]/51;
   DCTkvant.ar[7]=DCTkvant.ar[7]/61;
   DCTkvant.ar[8]=DCTkvant.ar[8]/12;
   DCTkvant.ar[9]=DCTkvant.ar[9]/12;
   DCTkvant.ar[10]=DCTkvant.ar[10]/14;
   DCTkvant.ar[11]=DCTkvant.ar[11]/19;
   DCTkvant.ar[12]=DCTkvant.ar[12]/26;
   DCTkvant.ar[13]=DCTkvant.ar[13]/58;
   DCTkvant.ar[14]=DCTkvant.ar[14]/60;
   DCTkvant.ar[15]=DCTkvant.ar[15]/55;
   DCTkvant.ar[16]=DCTkvant.ar[16]/14;
   DCTkvant.ar[17]=DCTkvant.ar[17]/13;
   DCTkvant.ar[18]=DCTkvant.ar[18]/16;
   DCTkvant.ar[19]=DCTkvant.ar[19]/24;
   DCTkvant.ar[20]=DCTkvant.ar[20]/40;
   DCTkvant.ar[21]=DCTkvant.ar[21]/57;
   DCTkvant.ar[22]=DCTkvant.ar[22]/69;
   DCTkvant.ar[23]=DCTkvant.ar[23]/56;
   DCTkvant.ar[24]=DCTkvant.ar[24]/14;
   DCTkvant.ar[25]=DCTkvant.ar[25]/17;
   DCTkvant.ar[26]=DCTkvant.ar[26]/22;
   DCTkvant.ar[27]=DCTkvant.ar[27]/29;
   DCTkvant.ar[28]=DCTkvant.ar[28]/51;
   DCTkvant.ar[29]=DCTkvant.ar[29]/87;
   DCTkvant.ar[30]=DCTkvant.ar[30]/80;
   DCTkvant.ar[31]=DCTkvant.ar[31]/62;
   DCTkvant.ar[32]=DCTkvant.ar[32]/18;
   DCTkvant.ar[33]=DCTkvant.ar[33]/22;
   DCTkvant.ar[34]=DCTkvant.ar[34]/37;
   DCTkvant.ar[35]=DCTkvant.ar[35]/56;
   DCTkvant.ar[36]=DCTkvant.ar[36]/68;
   DCTkvant.ar[37]=DCTkvant.ar[37]/109;
   DCTkvant.ar[38]=DCTkvant.ar[38]/103;
   DCTkvant.ar[39]=DCTkvant.ar[39]/77;
   DCTkvant.ar[40]=DCTkvant.ar[40]/24;
   DCTkvant.ar[41]=DCTkvant.ar[41]/35;
   DCTkvant.ar[42]=DCTkvant.ar[42]/55;
   DCTkvant.ar[43]=DCTkvant.ar[43]/64;
   DCTkvant.ar[44]=DCTkvant.ar[44]/81;
   DCTkvant.ar[45]=DCTkvant.ar[45]/104;
   DCTkvant.ar[46]=DCTkvant.ar[46]/113;
   DCTkvant.ar[47]=DCTkvant.ar[47]/92;
   DCTkvant.ar[48]=DCTkvant.ar[48]/49;
   DCTkvant.ar[49]=DCTkvant.ar[49]/64;
   DCTkvant.ar[50]=DCTkvant.ar[50]/78;
   DCTkvant.ar[51]=DCTkvant.ar[51]/87;
   DCTkvant.ar[52]=DCTkvant.ar[52]/103;
   DCTkvant.ar[53]=DCTkvant.ar[53]/121;
   DCTkvant.ar[54]=DCTkvant.ar[54]/120;
   DCTkvant.ar[55]=DCTkvant.ar[55]/101;
   DCTkvant.ar[56]=DCTkvant.ar[56]/72;
   DCTkvant.ar[57]=DCTkvant.ar[57]/92;
   DCTkvant.ar[58]=DCTkvant.ar[58]/95;
   DCTkvant.ar[59]=DCTkvant.ar[59]/98;
   DCTkvant.ar[60]=DCTkvant.ar[60]/112;
   DCTkvant.ar[61]=DCTkvant.ar[61]/100;
   DCTkvant.ar[62]=DCTkvant.ar[62]/103;
   DCTkvant.ar[63]=DCTkvant.ar[63]/99; */

   
   return DCTkvant;
}
// Inverse Discret Cosinus Transformation*--------------------------------------------
matric<float> IDCT(matric<float> &r, int M, int N, int NN)
{
   matric<float> DCTkvant(NN, NN);
   matric<float> kvant(NN, NN);
   float PI = 3.1415926535897932384626433832795;
   float Cu, Cv, AC;
   for (int u = 0; u < NN; u++)
   {
      for (int v = 0; v < NN; v++)
      {

         AC=0;
         for (int x = 0; x < NN; x++)
         {
            for (int y = 0; y < NN; y++)
            {
               if(x != 0)
               {Cu=1;}
               else
               {Cu=1/(pow(2, 0.5));}
               if(y != 0)
               {Cv=1;}
               else
               {Cv=1/(pow(2, 0.5));}
               if((y+M)<r.m && (x+N)<r.n) {AC += Cu*Cv*r.ar[(y+M)*r.n+x+N]*cos((2*(u)+1)*(x)*PI/(2*NN))*cos((2*(v)+1)*(y)*PI/(2*NN));}
            }
         }
       /*
         if(u==7){AC=AC*4;}
         if(u==6){AC=AC*1.8;}
         if(u==5){AC=AC*1.5;}
         if(u==4){AC=AC*1.2;}
         if(u==3){AC=AC*1.1;}
         if(u==2){AC=AC*1.05;}
         if(u==1){AC=AC*1.01;}

         if(v==7){AC=AC*4;}
         if(v==6){AC=AC*1.8;}
         if(v==5){AC=AC*1.25;}
         if(v==4){AC=AC*1.1;}
         if(v==3){AC=AC*1.05;}
         if(v==2){AC=AC*1.025;}
         if(v==1){AC=AC*1.005;}  */

         DCTkvant.ar[v*NN+u]=0.25*AC;
      }
   }
 /*  DCTkvant.ar[0]=DCTkvant.ar[0]*8;
   DCTkvant.ar[1]=DCTkvant.ar[1]*11;
   DCTkvant.ar[2]=DCTkvant.ar[2]*10;
   DCTkvant.ar[3]=DCTkvant.ar[3]*16;
   DCTkvant.ar[4]=DCTkvant.ar[4]*24;
   DCTkvant.ar[5]=DCTkvant.ar[5]*40;
   DCTkvant.ar[6]=DCTkvant.ar[6]*51;
   DCTkvant.ar[7]=DCTkvant.ar[7]*61;
   DCTkvant.ar[8]=DCTkvant.ar[8]*12;
   DCTkvant.ar[9]=DCTkvant.ar[9]*12;
   DCTkvant.ar[10]=DCTkvant.ar[10]*14;
   DCTkvant.ar[11]=DCTkvant.ar[11]*19;
   DCTkvant.ar[12]=DCTkvant.ar[12]*26;
   DCTkvant.ar[13]=DCTkvant.ar[13]*58;
   DCTkvant.ar[14]=DCTkvant.ar[14]*60;
   DCTkvant.ar[15]=DCTkvant.ar[15]*55;
   DCTkvant.ar[16]=DCTkvant.ar[16]*14;
   DCTkvant.ar[17]=DCTkvant.ar[17]*13;
   DCTkvant.ar[18]=DCTkvant.ar[18]*16;
   DCTkvant.ar[19]=DCTkvant.ar[19]*24;
   DCTkvant.ar[20]=DCTkvant.ar[20]*40;
   DCTkvant.ar[21]=DCTkvant.ar[21]*57;
   DCTkvant.ar[22]=DCTkvant.ar[22]*69;
   DCTkvant.ar[23]=DCTkvant.ar[23]*56;
   DCTkvant.ar[24]=DCTkvant.ar[24]*14;
   DCTkvant.ar[25]=DCTkvant.ar[25]*17;
   DCTkvant.ar[26]=DCTkvant.ar[26]*22;
   DCTkvant.ar[27]=DCTkvant.ar[27]*29;
   DCTkvant.ar[28]=DCTkvant.ar[28]*51;
   DCTkvant.ar[29]=DCTkvant.ar[29]*87;
   DCTkvant.ar[30]=DCTkvant.ar[30]*80;
   DCTkvant.ar[31]=DCTkvant.ar[31]*62;
   DCTkvant.ar[32]=DCTkvant.ar[32]*18;
   DCTkvant.ar[33]=DCTkvant.ar[33]*22;
   DCTkvant.ar[34]=DCTkvant.ar[34]*37;
   DCTkvant.ar[35]=DCTkvant.ar[35]*56;
   DCTkvant.ar[36]=DCTkvant.ar[36]*68;
   DCTkvant.ar[37]=DCTkvant.ar[37]*109;
   DCTkvant.ar[38]=DCTkvant.ar[38]*103;
   DCTkvant.ar[39]=DCTkvant.ar[39]*77;
   DCTkvant.ar[40]=DCTkvant.ar[40]*24;
   DCTkvant.ar[41]=DCTkvant.ar[41]*35;
   DCTkvant.ar[42]=DCTkvant.ar[42]*55;
   DCTkvant.ar[43]=DCTkvant.ar[43]*64;
   DCTkvant.ar[44]=DCTkvant.ar[44]*81;
   DCTkvant.ar[45]=DCTkvant.ar[45]*104;
   DCTkvant.ar[46]=DCTkvant.ar[46]*113;
   DCTkvant.ar[47]=DCTkvant.ar[47]*92;
   DCTkvant.ar[48]=DCTkvant.ar[48]*49;
   DCTkvant.ar[49]=DCTkvant.ar[49]*64;
   DCTkvant.ar[50]=DCTkvant.ar[50]*78;
   DCTkvant.ar[51]=DCTkvant.ar[51]*87;
   DCTkvant.ar[52]=DCTkvant.ar[52]*103;
   DCTkvant.ar[53]=DCTkvant.ar[53]*121;
   DCTkvant.ar[54]=DCTkvant.ar[54]*120;
   DCTkvant.ar[55]=DCTkvant.ar[55]*101;
   DCTkvant.ar[56]=DCTkvant.ar[56]*72;
   DCTkvant.ar[57]=DCTkvant.ar[57]*92;
   DCTkvant.ar[58]=DCTkvant.ar[58]*95;
   DCTkvant.ar[59]=DCTkvant.ar[59]*98;
   DCTkvant.ar[60]=DCTkvant.ar[60]*112;
   DCTkvant.ar[61]=DCTkvant.ar[61]*100;
   DCTkvant.ar[62]=DCTkvant.ar[62]*103;
   DCTkvant.ar[63]=DCTkvant.ar[63]*99;   */
   return DCTkvant;
}
//Fast Furye Transformation --------------------------------------------------
matric<int> FFT(matric<float> &r, int M, int N)
{
   matric<int> FFTmatric(8, 8);
   float Sum1, Sum2;
   for(int u=0; u<8; u++)
   {
      for (int v=0; v<8; v++)
      {
         Sum1=0;
         for(int y=0; y<8; y++)
         {
            for(int x=0; x<8; x++)
            {
               Sum1 += exp(-2*3.14*u*x/8) * r.ar[(y+M)*8+x+N] * exp(-2*3.14*v*y/8);
            }
         }
         FFTmatric.ar[v*8+u] = (int)(0.125*0.125*Sum1);
      }
   }

   return FFTmatric;
}

//Среднее значение------------------------------------------------------------
long double X(matric<float> &r, int M, int N, float R)
{
   long double X = 0, XX = 0;
   long int  Sum = 0, Sum2 = 0;
   for (int j = (M-R/2); j < (M+R/2); j++)
   {
      for (int k = (N-R/2); k < (N+R/2); k++)
         if((pow(j-M, 2) + pow(k-N, 2) - pow(R, 2)) < 1)
         {
            X += r.ar[j*r.n+k];
            Sum++;
         }
         XX += X / Sum;
         Sum2++;
         Sum = 0;
         X = 0;
   }
   XX = XX / Sum2;
   return XX;
};

//Мат ожидание ---------------------------------------------------------------
long double Mmn(matric<float> &r, int M, int N, float R)
{
   long int Sum = 0;
   long double x, Mmn = 0;
   x = X(r, M, N, R);
   for (int j = (M-R/2); j < (M+R/2); j++)
      for (int k = (N-R/2); k < (N+R/2); k++)
         if((pow(j-M, 2) + pow(k-N, 2) - pow(R, 2)) < 1)
         {
            Mmn += r.ar[j*r.n+k]*x;
            Sum++;
         }
   Mmn = Mmn/(Sum-1);
   return Mmn;
};

//Второй момент ---------------------------------------------------------------
long double DXmn(matric<float> &r, int M, int N, float R)
{
   long double Fm, x, D = 0;
   long int Sum = 0;
   x = X(r, M, N, R);
   Fm = Mmn(r, M, N, R);
   for (int j = (M-R/2); j < (M+R/2); j++)
      for (int k = (N-R/2); k < (N+R/2); k++)
         if((pow(j-M, 2) + pow(k-N, 2) - pow(R, 2)) < 1)
         {
            D += r.ar[j*r.n+k]*Fm*x;
            Sum++;
         }
   D = D/(Sum-1);
   return D;
};

//Регрессия -------------------------------------------------------------------
long double regres(matric<float> &r1, matric<float> &r2, int r1M, int r1N, int r2M, int r2N, float R)
{
   long double Regres = 0;

   Regres = pow(pow(DXmn(r1, r1M, r1N, R)-DXmn(r2, r2M, r2N, R), 2), 0.5);

   return Regres;
}

//Плоская корреляция-----------------------------------------------------------
long double Jxy(matric<float> &r1, matric<float> &r2, int r1M, int r1N, int r2M, int r2N, float R)
{
   long double J = 0;
   for (int j = (r1M-R/2); j < (r1M+R/2); j++)
      for (int k = (r1N-R/2); k < (r1N+R/2); k++)
         if((pow(j-r1M, 2) + pow(k-r1N, 2) - pow(R, 2)) < 1)
         {
            J += r1.ar[j*r1.n+k]*(j+k) - r2.ar[(j+r2M-r1M)*r2.n+(k+r2N-r1N)]*((j+r2M-r1M)+(k+r2N-r1N));
         }
         J = pow(J*J, 0.5);
   return J;
}

//Определение коэффициентов аффинного преобразования--------------------------
matric<float> MoveOxy(image &rimg, image &rimg2, int ColorFlag)
{
   int okresnost=14, R=16, Segment=3 , alfa=7, scale=3, Points=200;
   int KeyPoint = 0, goodPoints [400], X1, Y1, X2, Y2;
   matric<float> r(rimg.Y), r2(rimg2.Y);
   String txt;
   matric<float> AffinTransformation(2, 3);
   matric<float> AffinTransformation2(2, 3);
   vektor *v, *v1, *v2, *v3;
   v = new vektor[Points];
   v1 = new vektor[Points];
   v2 = new vektor[Points];
   v3 = new vektor[Points];

   for(int i = 0; i < 400; i++)
     goodPoints[i] = 0;

   bool f = false;
   matric<float> delta(3, 3), deltaA(3, 3), deltaB(3, 3), deltaC(3, 3), deltaD(3, 3), deltaE(3, 3), deltaF(3, 3);
   matric<float> delta2(3, 3), deltaA2(3, 3), deltaB2(3, 3), deltaC2(3, 3), deltaD2(3, 3), deltaE2(3, 3), deltaF2(3, 3);

   X1 = (okresnost*2); //if((float)r.n/2-200>X1){X1 = (int)(float)r.n/2-200;}   //????? ????? ?????????? ???? ???????????? ????????? ???????????
   Y1 = (okresnost*2); //if((float)r.m/2-300>Y1){Y1 = (int)(float)r.m/2-300;}   //????? ????? ?????????? ???? ???????????? ????????? ???????????
   X2 = rimg.Y.n-(2*okresnost); //if((float)r.n/2+200<X2){X2 = (int)(float)r.n/2+200;}  //????? ????? ?????????? ???? ???????????? ????????? ???????????
   Y2 = rimg.Y.m-(2*okresnost); //if((float)r.m/2+300<Y2){Y2 = (int)(float)r.m/2+300;}  //????? ????? ?????????? ???? ???????????? ????????? ???????????

   for (int j = Y1; j < Y2+1; j++)
      for (int k = X1; k < X2+1; k++)
      {
         f = true;
        for (int u = -okresnost; u < okresnost + 1; u++)
            for (int v = -okresnost; v < okresnost + 1; v++)
               if (rimg.Lnorm.ar[(j + u)*rimg.Lnorm.n +(k + v)] > (rimg.Lnorm.ar[j*rimg.Lnorm.n+k]))
                  {f = false; u = v =  okresnost;}
         if (f == true)
         {
            v1[KeyPoint].z = j;
            v1[KeyPoint].x = k;
            KeyPoint++;
            if (KeyPoint>=Points) {j = Y2; k = X2;}
         }
      }

  float *mean, *std_dev, *skewness, *neighbor_contrast, *beta_coeff;
  mean = new float[KeyPoint];
  std_dev = new float[KeyPoint];
  skewness= new float[KeyPoint];
  neighbor_contrast= new float[KeyPoint];
  beta_coeff= new float[KeyPoint];

  float *mean2, *std_dev2, *skewness2, *neighbor_contrast2, *beta_coeff2;
  mean2 = new float[KeyPoint];
  std_dev2 = new float[KeyPoint];
  skewness2 = new float[KeyPoint];
  neighbor_contrast2 = new float[KeyPoint];
  beta_coeff2 = new float[KeyPoint];

  float meanMIN=65000, meanMAX=0;
  float std_devMIN=65000, std_devMAX=0;
  float skewnessMIN=65000, skewnessMAX=0;
  float neighbor_contrastMIN=65000, neighbor_contrastMAX=0;
  float beta_coeffMIN=65000, beta_coeffMAX=0;
  float *d;
  d = new float[KeyPoint];
  float dMIN=65000, dMAX=0, dSUM=0, dTOL;

  int incr = 1;
  float sxx=0, sxy=0, x=0, diff=0, mean_rad=0, SumN=0, col_dist=0, row_dist_sq=0, center_row=0, center_col=0;

for (int i = 0; i < KeyPoint; i++)
{
      X1 = (int)v1[i].x - Segment;
      Y1 = (int)v1[i].z - Segment;
      X2 = (int)v1[i].x + Segment;
      Y2 = (int)v1[i].z + Segment;

   SumN = ((Y2+1)*rimg.Y.n+(X2+1) - Y1*rimg.Y.n+X1);
   if(SumN)
   {
      center_row = (Y2+1 - Y1)/2 + Y1;
      center_col = (X2+1 - X1)/2 + X1;
      for (int j = Y1; j < Y2+1; j++)
      {
         row_dist_sq = (j-center_row)*(j-center_row);
         for (int k = X1; k < X2+1; k++)
         {
            col_dist = k-center_col;
            mean[i] += rimg.Y.ar[j*rimg.Y.n+k];
            mean_rad += sqrt(row_dist_sq + col_dist*col_dist);
         }
      }
      mean[i] = mean[i]/SumN;
      mean_rad /= SumN;
      for (int j = Y1; j < Y2+1; j++)
      {
         for (int k = X1; k < X2+1; k++)
            {
               x = mean_rad - sqrt(row_dist_sq + col_dist*col_dist);
               sxy += x*(rimg.Y.ar[j*rimg.Y.n+k] - mean[i]);
	            sxx += x*x;

               diff = rimg.Y.ar[j*rimg.Y.n+k]-mean[i];
               std_dev[i] += diff*diff;
               skewness[i] += diff*diff*diff;

            }
      }
      for (int j = Y1+incr; j < Y2+1; j++)
         for (int k = X1+incr; k < X2+1; k++)
            neighbor_contrast[i] += (float)abs((int)((rimg.Y.ar[j*rimg.Y.n+k] - rimg.Y.ar[(j-incr)*rimg.Y.n+k])*100))/100 + (float)abs((int)((rimg.Y.ar[j*rimg.Y.n+k] - rimg.Y.ar[j*rimg.Y.n+(k-incr)])*100))/100;

      if(std_dev[i] < 1e-10 && (-std_dev[i]) < 1e-6){std_dev[i] =0.0;}
      else{std_dev[i] = sqrt(std_dev[i]/SumN);}
      skewness[i] /= SumN;
      skewness[i] = skewness[i]/(std_dev[i]*std_dev[i]*std_dev[i]);
      neighbor_contrast[i] /= SumN;
      if (sxx < 1e-10) {beta_coeff[i] = 0.0;}
      else {beta_coeff[i] = sxy/sxx;}
   }
   if(i==0)
   {
      meanMIN = meanMAX = mean[i];
      std_devMIN = std_devMAX = std_dev[i];
      skewnessMIN = skewnessMAX = skewness[i];
      neighbor_contrastMIN = neighbor_contrastMAX = neighbor_contrast[i];
      beta_coeffMIN = beta_coeffMAX = beta_coeff[i];
   }
   if(meanMIN > mean[i])meanMIN = mean[i];
   if(meanMAX < mean[i])meanMAX = mean[i];
   if(std_devMIN > std_dev[i])std_devMIN = std_dev[i];
   if(std_devMAX < std_dev[i])std_devMAX = std_dev[i];
   if(skewnessMIN > skewness[i])skewnessMIN = skewness[i];
   if(skewnessMAX < skewness[i])skewnessMAX = skewness[i];
   if(neighbor_contrastMIN > neighbor_contrast[i])neighbor_contrastMIN = neighbor_contrast[i];
   if(neighbor_contrastMAX < neighbor_contrast[i])neighbor_contrastMAX = neighbor_contrast[i];
   if(beta_coeffMIN > beta_coeff[i])beta_coeffMIN = beta_coeff[i];
   if(beta_coeffMAX < beta_coeff[i])beta_coeffMAX = beta_coeff[i];

}
for (int i = 0; i < KeyPoint; i++)
{
  mean[i] = (mean[i]-meanMIN)/(meanMAX-meanMIN);
  std_dev[i] = (std_dev[i]-std_devMIN)/(std_devMAX-std_devMIN);
  skewness[i] = (skewness[i]-skewnessMIN)/(skewnessMAX-skewnessMIN);
  neighbor_contrast[i] = (neighbor_contrast[i]-neighbor_contrastMIN)/(neighbor_contrastMAX-neighbor_contrastMIN);
  beta_coeff[i] = (beta_coeff[i]-beta_coeffMIN)/(beta_coeffMAX-beta_coeffMIN);
}
dTOL=65532;
bool badPoint;
int  possiblyPoint[1000];

//*****************************************************************************

for (int uu = -R; uu < R+1; uu++)
{
	for (int vv = -R; vv < R+1; vv++)
	{
      possiblyPoint[0]=0;
		for (int i = 0; i < KeyPoint; i++)
		{
			sxx=sxy=x=diff=mean_rad=SumN=col_dist=row_dist_sq=center_row=center_col=0;
         badPoint = false;

         X1 = (float)v1[i].x  - Segment + vv;
			Y1 = (float)v1[i].z  - Segment + uu;
			X2 = (float)v1[i].x  + Segment + vv;
			Y2 = (float)v1[i].z  + Segment + uu;

         if((X2>rimg2.R.n)|| (Y2>rimg2.R.m)||(X1<1)||(Y1<1))
         {
            badPoint = true;
         }
         else
         {
            possiblyPoint[0]++;
            possiblyPoint[possiblyPoint[0]]=i;
         }

         if(!badPoint)
			SumN = ((Y2+1)*rimg2.Y.n+(X2+1) - Y1*rimg2.Y.n+X1);

         if(!badPoint)
			if(SumN)
			{
				center_row = (Y2+1 - Y1)/2 + Y1;
				center_col = (X2+1 - X1)/2 + X1;
				for (int j = Y1; j < Y2+1; j++)
				{
					row_dist_sq = (j-center_row)*(j-center_row);
					for (int k = X1; k < X2+1; k++)
					{
						col_dist = k-center_col;
						mean2[i] += rimg2.Y.ar[j*rimg2.Y.n+k];
						mean_rad += sqrt(row_dist_sq + col_dist*col_dist);
					}
				}
				mean2[i] = mean2[i]/SumN;
				mean_rad /= SumN;
				for (int j = Y1; j < Y2+1; j++)
				{
					for (int k = X1; k < X2+1; k++)
					{
					   x = mean_rad - sqrt(row_dist_sq + col_dist*col_dist);
					   sxy += x*(rimg2.Y.ar[j*rimg2.Y.n+k] - mean2[i]);
						sxx += x*x;

					   diff = rimg2.Y.ar[j*rimg2.Y.n+k]-mean2[i];
					   std_dev2[i] += diff*diff;
					   skewness2[i] += diff*diff*diff;
					}
				}
				for (int j = Y1+incr; j < Y2+1; j++)
					for (int k = X1+incr; k < X2+1; k++)
						neighbor_contrast[i] += (float)abs((int)((rimg2.Y.ar[j*rimg2.Y.n+k] - rimg2.Y.ar[(j-incr)*rimg2.Y.n+k])*100))/100 + (float)abs((int)((rimg2.Y.ar[j*rimg2.Y.n+k] - rimg2.Y.ar[j*rimg2.Y.n+(k-incr)])*100))/100;

				if(std_dev2[i] < 1e-10 && (-std_dev2[i]) < 1e-6){std_dev2[i] = 0.0;}
            else{std_dev2[i] = sqrt(std_dev2[i]/SumN);}
				skewness2[i] /= SumN;
				skewness2[i] = skewness2[i]/(std_dev2[i]*std_dev2[i]*std_dev2[i]);
				neighbor_contrast2[i] /= SumN;
				if (sxx < 1e-10) {beta_coeff2[i] = 0.0;}
				else {beta_coeff2[i] = sxy/sxx;}
		   }
         if(!badPoint)
		   if(i==0)
		   {
			  meanMIN = meanMAX = mean2[i];
			  std_devMIN = std_devMAX = std_dev2[i];
			  skewnessMIN = skewnessMAX = skewness2[i];
			  neighbor_contrastMIN = neighbor_contrastMAX = neighbor_contrast2[i];
			  beta_coeffMIN = beta_coeffMAX = beta_coeff2[i];
		   }
         if(!badPoint)
         {
   		   if(meanMIN > mean2[i])meanMIN = mean2[i];
   		   if(meanMAX < mean2[i])meanMAX = mean2[i];
	   	   if(std_devMIN > std_dev2[i])std_devMIN = std_dev2[i];
		      if(std_devMAX < std_dev2[i])std_devMAX = std_dev2[i];
   		   if(skewnessMIN > skewness2[i])skewnessMIN = skewness2[i];
	   	   if(skewnessMAX < skewness2[i])skewnessMAX = skewness2[i];
		      if(neighbor_contrastMIN > neighbor_contrast2[i])neighbor_contrastMIN = neighbor_contrast2[i];
		      if(neighbor_contrastMAX < neighbor_contrast2[i])neighbor_contrastMAX = neighbor_contrast2[i];
   		   if(beta_coeffMIN > beta_coeff2[i])beta_coeffMIN = beta_coeff2[i];
	   	   if(beta_coeffMAX < beta_coeff2[i])beta_coeffMAX = beta_coeff2[i];
         }

		}
      if(!badPoint)
		for (int ii = 1; ii < possiblyPoint[0]; ii++)
		{
        int i;
        i = possiblyPoint[ii];
		  mean2[i] = (mean2[i]-meanMIN)/(meanMAX-meanMIN);
		  std_dev2[i] = (std_dev2[i]-std_devMIN)/(std_devMAX-std_devMIN);
		  skewness2[i] = (skewness2[i]-skewnessMIN)/(skewnessMAX-skewnessMIN);
		  neighbor_contrast2[i] = (neighbor_contrast2[i]-neighbor_contrastMIN)/(neighbor_contrastMAX-neighbor_contrastMIN);
		  beta_coeff2[i] = (beta_coeff2[i]-beta_coeffMIN)/(beta_coeffMAX-beta_coeffMIN);
		}
		dSUM=0;
      if(!badPoint)
		for (int ii = 1; ii < possiblyPoint[0]; ii++)
		{
         int i;
         i = possiblyPoint[ii];
		   d[i] =  (float)abs((int)((mean[i]-mean2[i])*100))/100+(float)abs((int)((std_dev[i]-std_dev2[i])*100))/100+(float)abs((int)((skewness[i]-skewness2[i])*100))/100+(float)abs((int)((neighbor_contrast[i]-neighbor_contrast2[i])*100))/100+(float)abs((int)((beta_coeff[i]-beta_coeff2[i])*100))/100;
		   if(i==0)
		   {
			  dMIN=dMAX=d[i];
		   }
		   if(dMIN>d[i])dMIN=d[i];
		   if(dMAX<d[i])dMAX=d[i];
		   dSUM+=d[i];
		}
      if(possiblyPoint[0]>3){dSUM /= possiblyPoint[0]>3;}
      if(!badPoint && possiblyPoint[0]>3)
		if(dSUM<dTOL)
		{
         dTOL=dSUM;
         if (goodPoints[0]>5)
         {
            goodPoints[0]=0;
            for(int i = 1; i < 4 ; i++)
               goodPoints[i] = goodPoints[i+3];
         }

			for (int ii = 1; ii < possiblyPoint[0]; ii++)
			{
            int i;
            i = possiblyPoint[ii];
			   d[i]=(d[i]-dMIN)/(dMAX-dMIN);
			   if(d[i]<0.3)
			   {
				  v3[i].x = v1[i].x + vv;
				  v3[i].z = v1[i].z + uu;
              v2=v3;
              v[i] =(v2[i] - v1[i]);
 			     goodPoints[0]++;
			     goodPoints[goodPoints[0]] = i;
			   }
			}
		}
	}
}

for(int aa = (100-alfa); aa < (101+alfa); aa++)
{
   for(int cc = (100-alfa); cc < (101+alfa); cc++)
   {
      for(int bb = (-scale); bb < (scale+1); bb++)
      {
         for(int ga = (-scale); ga < (scale+1); ga++)
         {
            possiblyPoint[0] = 0;
            for (int i = 0; i < KeyPoint; i++)
            {
                  badPoint = false;
				      X1 = (((float)aa)/100)*(float)v3[i].x + (((float)bb)/100)*(float)v3[i].z - Segment;
				      Y1 = (((float)cc)/100)*(float)v3[i].z + (((float)ga)/100)*(float)v3[i].x - Segment;
				      X2 = (((float)aa)/100)*(float)v3[i].x + (((float)bb)/100)*(float)v3[i].z + Segment;
				      Y2 = (((float)cc)/100)*(float)v3[i].z + (((float)ga)/100)*(float)v3[i].x + Segment;

                  if((X2>rimg2.R.n)|| (Y2>rimg2.R.m)||(X1<1)||(Y1<1))
                  {
                     badPoint = true;
                  }
                  else
                  {
                     possiblyPoint[0]++;
                     possiblyPoint[possiblyPoint[0]]=i;
                  }
               if(!badPoint)
				   SumN = ((Y2+1)*rimg2.Y.n+(X2+1) - Y1*rimg2.Y.n+X1);
               if(!badPoint)
				   if(SumN)
				   {
					  center_row = (Y2+1 - Y1)/2 + Y1;
					  center_col = (X2+1 - X1)/2 + X1;
					  for (int j = Y1; j < Y2+1; j++)
					  {
						 row_dist_sq = (j-center_row)*(j-center_row);
						 for (int k = X1; k < X2+1; k++)
						 {
							 col_dist = k - center_col;
							 mean2[i] += rimg2.Y.ar[j*rimg2.Y.n+k];
						    mean_rad += sqrt(row_dist_sq + col_dist*col_dist);
						 }
					  }
					  mean2[i] = mean2[i]/SumN;
					  mean_rad /= SumN;
					  for (int j = Y1; j < Y2+1; j++)
					  {
						 for (int k = X1; k < X2+1; k++)
						 {
							x = mean_rad - sqrt(row_dist_sq + col_dist*col_dist);
							sxy += x*(rimg2.Y.ar[j*rimg2.Y.n+k] - mean2[i]);
							sxx += x*x;
							diff = rimg2.Y.ar[j*rimg2.Y.n+k]-mean2[i];
							std_dev2[i] += diff*diff;
							skewness2[i] += diff*diff*diff;
						 }
					  }
					  for (int j = Y1+incr; j < Y2+1; j++)
						 for (int k = X1+incr; k < X2+1; k++)
							neighbor_contrast[i] += (float)abs((int)((rimg2.Y.ar[j*rimg2.Y.n+k] - rimg2.Y.ar[(j-incr)*rimg2.Y.n+k])*100))/100 + (float)abs((int)((rimg2.Y.ar[j*rimg2.Y.n+k] - rimg2.Y.ar[j*rimg2.Y.n+(k-incr)])*100))/100;
                 if(std_dev2[i] < 1e-10 && (-std_dev2[i]) < 1e-6){std_dev2[i] = 0.0;}
                 else{std_dev2[i] = sqrt(std_dev2[i]/SumN);}
					  skewness2[i] /= SumN;
					  if(std_dev2[i]*std_dev2[i]*std_dev2[i] < 1e-6){skewness2[i] = 0.0;}
                 else{skewness2[i] = skewness2[i]/(std_dev2[i]*std_dev2[i]*std_dev2[i]);}
					  neighbor_contrast2[i] /= SumN;
					  if (sxx < 1e-10) {beta_coeff2[i] = 0.0;}
					  else {beta_coeff2[i] = sxy/sxx;}
				   }
				   if(i==0)
				   {
					  meanMIN = meanMAX = mean2[i];
					  std_devMIN = std_devMAX = std_dev2[i];
					  skewnessMIN = skewnessMAX = skewness2[i];
					  neighbor_contrastMIN = neighbor_contrastMAX = neighbor_contrast2[i];
					  beta_coeffMIN = beta_coeffMAX = beta_coeff2[i];
				   }
               if(!badPoint)
               {
				      if(meanMIN > mean2[i])meanMIN = mean2[i];
			   	   if(meanMAX < mean2[i])meanMAX = mean2[i];
		   		   if(std_devMIN > std_dev2[i])std_devMIN = std_dev2[i];
	   			   if(std_devMAX < std_dev2[i])std_devMAX = std_dev2[i];
   				   if(skewnessMIN > skewness2[i])skewnessMIN = skewness2[i];
				      if(skewnessMAX < skewness2[i])skewnessMAX = skewness2[i];
			   	   if(neighbor_contrastMIN > neighbor_contrast2[i])neighbor_contrastMIN = neighbor_contrast2[i];
		   		   if(neighbor_contrastMAX < neighbor_contrast2[i])neighbor_contrastMAX = neighbor_contrast2[i];
	   			   if(beta_coeffMIN > beta_coeff2[i])beta_coeffMIN = beta_coeff2[i];
   				   if(beta_coeffMAX < beta_coeff2[i])beta_coeffMAX = beta_coeff2[i];
               }
            }
            if(!badPoint)
				for (int ii = 1; ii < possiblyPoint[0]; ii++)
				{
               int i;
               i = possiblyPoint[ii];
				   mean2[i] = (mean2[i]-meanMIN)/(meanMAX-meanMIN);
				   std_dev2[i] = (std_dev2[i]-std_devMIN)/(std_devMAX-std_devMIN);
				   skewness2[i] = (skewness2[i]-skewnessMIN)/(skewnessMAX-skewnessMIN);
				   neighbor_contrast2[i] = (neighbor_contrast2[i]-neighbor_contrastMIN)/(neighbor_contrastMAX-neighbor_contrastMIN);
				   beta_coeff2[i] = (beta_coeff2[i]-beta_coeffMIN)/(beta_coeffMAX-beta_coeffMIN);
				}
				dSUM=0;
            if(!badPoint)
				for (int ii = 1; ii < possiblyPoint[0]; ii++)
				{
               int i;
               i = possiblyPoint[ii];
					d[i] =  (float)abs((int)((mean[i]-mean2[i])*100))/100+(float)abs((int)((std_dev[i]-std_dev2[i])*100))/100+(float)abs((int)((skewness[i]-skewness2[i])*100))/100+(float)abs((int)((neighbor_contrast[i]-neighbor_contrast2[i])*100))/100+(float)abs((int)((beta_coeff[i]-beta_coeff2[i])*100))/100;
					if(i==0){dMIN=dMAX=d[i];}
					if(dMIN>d[i])dMIN=d[i];
					if(dMAX<d[i])dMAX=d[i];
					dSUM+=d[i];
				}
            if(possiblyPoint[0]>3){dSUM/=possiblyPoint[0];}
            if(!badPoint && possiblyPoint[0]>3)
				if(dSUM<=dTOL)
				{
               dTOL=dSUM;
				   //goodPoints[0]=0;
					for (int ii = 1; ii < possiblyPoint[0]; ii++)
					{
                  int i;
                  i = possiblyPoint[ii];
						d[i]=(d[i]-dMIN)/(dMAX-dMIN);
						if(d[i]<0.2 && v3[i].x>1 && v3[i].z>1)
						{
						   v2[i].x = (((float)aa)/100)*(int)v3[i].x + (((float)bb)/100)*(int)v3[i].z;
						   v2[i].z = (((float)cc)/100)*(int)v3[i].z + (((float)ga)/100)*(int)v3[i].x;
						   v[i] = (v2[i] - v1[i]);
						}
      			}
				}
	      }
      }
   }
}
//******************************************************************************
   Form1->Memo1->Lines->Add("Хороших точек: " + (String)(goodPoints[0]-2));
   delta = deltaA = deltaB = deltaC = deltaD = deltaE = deltaF = delta *0;
   for (int i = 1; i < (goodPoints[0]-2); i++)
   {
      delta.ar[0*delta.n+0] += v1[goodPoints[i]].x;
      delta.ar[1*delta.n+0] += v1[goodPoints[i+1]].x;
      delta.ar[2*delta.n+0] += v1[goodPoints[i+2]].x;
      delta.ar[0*delta.n+1] += v1[goodPoints[i]].z;
      delta.ar[1*delta.n+1] += v1[goodPoints[i+1]].z;
      delta.ar[2*delta.n+1] += v1[goodPoints[i+2]].z;
      delta.ar[0*delta.n+2] = delta.ar[1*delta.n+2] = delta.ar[2*delta.n+2] += 1;

      deltaA.ar[0*deltaA.n+0] += v2[goodPoints[i]].x;
      deltaA.ar[1*deltaA.n+0] += v2[goodPoints[i+1]].x;
      deltaA.ar[2*deltaA.n+0] += v2[goodPoints[i+2]].x;
      deltaA.ar[0*deltaA.n+1] += v1[goodPoints[i]].z;
      deltaA.ar[1*deltaA.n+1] += v1[goodPoints[i+1]].z;
      deltaA.ar[2*deltaA.n+1] += v1[goodPoints[i+2]].z;
      deltaA.ar[0*deltaA.n+2] = deltaA.ar[1*deltaA.n+2] = deltaA.ar[2*deltaA.n+2] += 1;

      deltaB.ar[0*deltaB.n+0] += v1[goodPoints[i]].x;
      deltaB.ar[1*deltaB.n+0] += v1[goodPoints[i+1]].x;
      deltaB.ar[2*deltaB.n+0] += v1[goodPoints[i+2]].x;
      deltaB.ar[0*deltaB.n+1] += v2[goodPoints[i]].x;
      deltaB.ar[1*deltaB.n+1] += v2[goodPoints[i+1]].x;
      deltaB.ar[2*deltaB.n+1] += v2[goodPoints[i+2]].x;
      deltaB.ar[0*deltaB.n+2] = deltaB.ar[1*deltaB.n+2] = deltaB.ar[2*deltaB.n+2] += 1;

      deltaE.ar[0*deltaE.n+0] += v1[goodPoints[i]].x;
      deltaE.ar[1*deltaE.n+0] += v1[goodPoints[i+1]].x;
      deltaE.ar[2*deltaE.n+0] += v1[goodPoints[i+2]].x;
      deltaE.ar[0*deltaE.n+1] += v1[goodPoints[i]].z;
      deltaE.ar[1*deltaE.n+1] += v1[goodPoints[i+1]].z;
      deltaE.ar[2*deltaE.n+1] += v1[goodPoints[i+2]].z;
      deltaE.ar[0*deltaE.n+2] += v2[goodPoints[i]].x;
      deltaE.ar[1*deltaE.n+2] += v2[goodPoints[i+1]].x;
      deltaE.ar[2*deltaE.n+2] += v2[goodPoints[i+2]].x;

      deltaC.ar[0*deltaC.n+0] += v2[goodPoints[i]].z;
      deltaC.ar[1*deltaC.n+0] += v2[goodPoints[i+1]].z;
      deltaC.ar[2*deltaC.n+0] += v2[goodPoints[i+2]].z;
      deltaC.ar[0*deltaC.n+1] += v1[goodPoints[i]].z;
      deltaC.ar[1*deltaC.n+1] += v1[goodPoints[i+1]].z;
      deltaC.ar[2*deltaC.n+1] += v1[goodPoints[i+2]].z;
      deltaC.ar[0*deltaC.n+2] = deltaC.ar[1*deltaC.n+2] = deltaC.ar[2*deltaC.n+2] += 1;

      deltaD.ar[0*deltaD.n+0] += v1[goodPoints[i]].x;
      deltaD.ar[1*deltaD.n+0] += v1[goodPoints[i+1]].x;
      deltaD.ar[2*deltaD.n+0] += v1[goodPoints[i+2]].x;
      deltaD.ar[0*deltaD.n+1] += v2[goodPoints[i]].z;
      deltaD.ar[1*deltaD.n+1] += v2[goodPoints[i+1]].z;
      deltaD.ar[2*deltaD.n+1] += v2[goodPoints[i+2]].z;
      deltaD.ar[0*deltaD.n+2] = deltaD.ar[1*deltaD.n+2] = deltaD.ar[2*deltaD.n+2] += 1;

      deltaF.ar[0*deltaF.n+0] += v1[goodPoints[i]].x;
      deltaF.ar[1*deltaF.n+0] += v1[goodPoints[i+1]].x;
      deltaF.ar[2*deltaF.n+0] += v1[goodPoints[i+2]].x;
      deltaF.ar[0*deltaF.n+1] += v1[goodPoints[i]].z;
      deltaF.ar[1*deltaF.n+1] += v1[goodPoints[i+1]].z;
      deltaF.ar[2*deltaF.n+1] += v1[goodPoints[i+2]].z;
      deltaF.ar[0*deltaF.n+2] += v2[goodPoints[i]].z;
      deltaF.ar[1*deltaF.n+2] += v2[goodPoints[i+1]].z;
      deltaF.ar[2*deltaF.n+2] += v2[goodPoints[i+2]].z;

      //????? ?????????? ???? ??????
      Form1->Image2->Canvas->Refresh();
      Form1->Image2->Canvas->MoveTo(v2[goodPoints[i]].x, v2[goodPoints[i]].z);

      switch (ColorFlag){
      case 0: Form1->Image2->Canvas->Pen->Color = clRed; break;
      case 1: Form1->Image2->Canvas->Pen->Color = clGreen; break;
      case 2: Form1->Image2->Canvas->Pen->Color = clBlue;  break;
      case 3: Form1->Image2->Canvas->Pen->Color = clYellow; break;  }
      Form1->Image2->Canvas->LineTo(v1[goodPoints[i]].x, v1[goodPoints[i]].z);
    }
   if(delta.D()>(1e-15) || -delta.D()>(1e-15))
   {
      AffinTransformation.ar[0*AffinTransformation.n+0] = deltaA.D()/delta.D();
      AffinTransformation.ar[0*AffinTransformation.n+1] = deltaB.D()/delta.D();
      AffinTransformation.ar[1*AffinTransformation.n+0] = deltaC.D()/delta.D();
      AffinTransformation.ar[1*AffinTransformation.n+1] = deltaD.D()/delta.D();
      AffinTransformation.ar[0*AffinTransformation.n+2] = deltaE.D()/delta.D();
      AffinTransformation.ar[1*AffinTransformation.n+2] = deltaF.D()/delta.D();
   }
   else
   {
      AffinTransformation.ar[0*AffinTransformation.n+0] = 1.0;
      AffinTransformation.ar[0*AffinTransformation.n+1] = 0;
      AffinTransformation.ar[1*AffinTransformation.n+0] = 0;
      AffinTransformation.ar[1*AffinTransformation.n+1] = 1.0;
      AffinTransformation.ar[0*AffinTransformation.n+2] = 0;
      AffinTransformation.ar[1*AffinTransformation.n+2] = 0;
   }
   delete d, v, v1, v2, v3;
   return AffinTransformation;
}

//Аффинное преобразование-----------------------------------------------------
void AffinTransformation(image &img, matric<float> &affin)
{
   float a, b, c, d, e, f, jaff, kaff;
   a = affin.ar[0*3+0];
   b = affin.ar[0*3+1];
   e = affin.ar[0*3+2];
   c = affin.ar[1*3+0];
   d = affin.ar[1*3+1];
   f = affin.ar[1*3+2];

   int m, n, x[4], y[4], Xmin=0, Ymin=0, Xmax=0, Ymax=0;
   x[0]=(-e)/a;
   y[0]=(-f)/d;
   x[1]=(-b*(img.R.m)-e)/a;
   y[1]=((img.R.m)-f)/d;
   x[2]=((img.R.n)-e)/a;
   y[2]=(-c*(img.R.n)-f)/d;
   x[3]=((img.R.n)-b*(img.R.m)-e)/a;
   y[3]=((img.R.m)-c*(img.R.n)-f)/d;

   for(int i=0; i<4; i++)
   {
      if(Xmin>x[i]){Xmin=x[i];}
      if(Ymin>y[i]){Ymin=y[i];}
      if(Xmax<x[i]){Xmax=x[i];}
      if(Ymax<y[i]){Ymax=y[i];}
   }
   m = (int)(Ymax-Ymin);
   n = (int)(Xmax-Xmin);
   Form1->Image1->Canvas->Rectangle(0,0,800,600);
   for (int j = 0; j <800; j++)
      for (int k = 0; k < 600; k++)
         {
            jaff =c*j + d*k + f;
            kaff =a*j + b*k + e;
            if(kaff>0 && jaff>0 && kaff<img.R.n && jaff<img.R.m)
            {
               Form1->Image1->Canvas->Pixels[(int)j][(int)k] = (int)img.R.ar[(int)jaff*img.R.n + (int)kaff] + ((int)(img.G.ar[(int)jaff*img.R.n + (int)kaff]) << 8) + ((int)(img.B.ar[(int)jaff*img.R.n + (int)kaff]) << 16);
            }
         }
};
/*
//Аффинное преобразование 3D--------------------------------------------------
void AffinTransformation3D(image &img, matric<float> &affin)
{
   float a, b, c, d, e, f, g, h, k, jaff, kaff;
   a = affin.ar[0*3+0];
   b = affin.ar[0*3+1];
   e = affin.ar[0*3+2];
   c = affin.ar[1*3+0];
   d = affin.ar[1*3+1];
   f = affin.ar[1*3+2];
   g = affin.ar[2*3+0];
   h = affin.ar[2*3+1];
   k = affin.ar[2*3+2];

   int m, n, x[4], y[4], Xmin=0, Ymin=0, Xmax=0, Ymax=0;
   x[0]=(-e)/a;
   y[0]=(-f)/d;
   x[1]=(-b*(img.R.m)-e)/a;
   y[1]=((img.R.m)-f)/d;
   x[2]=((img.R.n)-e)/a;
   y[2]=(-c*(img.R.n)-f)/d;
   x[3]=((img.R.n)-b*(img.R.m)-e)/a;
   y[3]=((img.R.m)-c*(img.R.n)-f)/d;

   for(int i=0; i<4; i++)
   {
      if(Xmin>x[i]){Xmin=x[i];}
      if(Ymin>y[i]){Ymin=y[i];}
      if(Xmax<x[i]){Xmax=x[i];}
      if(Ymax<y[i]){Ymax=y[i];}
   }
   m = (int)(Ymax-Ymin);
   n = (int)(Xmax-Xmin);
   Form1->Image1->Canvas->Rectangle(0,0,800,600);
   for (int j = 0; j <800; j++)
      for (int k = 0; k < 600; k++)
         {
            jaff =c*j + d*k + f;
            kaff =a*j + b*k + e;
            if(kaff>0 && jaff>0 && kaff<img.R.n && jaff<img.R.m)
            {
               Form1->Image1->Canvas->Pixels[(int)j][(int)k] = (int)img.R.ar[(int)jaff*img.R.n + (int)kaff] + ((int)(img.G.ar[(int)jaff*img.R.n + (int)kaff]) << 8) + ((int)(img.B.ar[(int)jaff*img.R.n + (int)kaff]) << 16);
            }
         }
};    */
//Склеивание изображений-------------------------------------------------------
image Glue(image &img1, image &img2, matric<float> &affin)
{
   float a, b, c, d, e, f;
   int jaff, kaff;
   a = affin.ar[0*affin.n+0];
   b = affin.ar[0*affin.n+1];
   e = affin.ar[0*affin.n+2];
   c = affin.ar[1*affin.n+0];
   d = affin.ar[1*affin.n+1];
   f = affin.ar[1*affin.n+2];

   int m, n, x[8], y[8], Xmin=0, Ymin=0, Xmax=0, Ymax=0;

   x[0]=0;
   y[0]=0;
   x[1]=0;
   y[1]=img1.R.m;
   x[2]=img1.R.n;
   y[2]=0;
   x[3]=img1.R.n;
   y[3]=img1.R.m;
   x[4]=(-e)/a-img1.dn;
   y[4]=(-f)/d-img1.dm;
   y[5]=((img2.R.m)-f)/d-img1.dm;
   x[6]=((img2.R.n)-e)/a-img1.dn;
   y[6]=(-c*(img2.R.n)-f)/d-img1.dm;
   x[7]=((img2.R.n)-b*(img2.R.m)-e)/a -img1.dn;
   y[7]=((img2.R.m)-c*(img2.R.n)-f)/d -img1.dm;

   for(int i=0; i<8; i++)
   {
      if(Xmin>x[i]){Xmin=x[i];}
      if(Ymin>y[i]){Ymin=y[i];}
      if(Xmax<x[i]){Xmax=x[i];}
      if(Ymax<y[i]){Ymax=y[i];}
   }
   m = (int)(Ymax-Ymin);
   n = (int)(Xmax-Xmin);
   image temp, temp2;
   temp.B = temp.G = temp.R = matric<float>(m, n)*0;
   temp2.B = temp2.G = temp2.R = matric<float>(m, n)*0;
   temp.dm = img1.dm;
   temp.dn = img1.dn;
   temp.dm = temp.dm+Ymin;
   temp.dn = temp.dn+Xmin;

    for (int j = 1; j < m; j++)
      for (int k = 1; k < n; k++)
      {
        if((j<img1.R.m)&&(k<img1.R.n))
        {
         temp.R.ar[((j) - Ymin)*temp.R.n+((k)-Xmin)] = img1.R.ar[(j)*img1.R.n+(k)];
         temp.G.ar[((j) - Ymin)*temp.G.n+((k)-Xmin)] = img1.G.ar[(j)*img1.G.n+(k)];
         temp.B.ar[((j) - Ymin)*temp.B.n+((k)-Xmin)] = img1.B.ar[(j)*img1.B.n+(k)];
      //----------------------------------------------------------------------------------
        /* temp.R.ar[((j+1) - Ymin)*temp.R.n+((k)-Xmin)] = img1.R.ar[(j)*img1.R.n+(k)];
         temp.G.ar[((j+1) - Ymin)*temp.G.n+((k)-Xmin)] = img1.G.ar[(j)*img1.G.n+(k)];
         temp.B.ar[((j+1) - Ymin)*temp.B.n+((k)-Xmin)] = img1.B.ar[(j)*img1.B.n+(k)];

         temp.R.ar[((j) - Ymin)*temp.R.n+((k+1)-Xmin)] = img1.R.ar[(j)*img1.R.n+(k)];
         temp.G.ar[((j) - Ymin)*temp.G.n+((k+1)-Xmin)] = img1.G.ar[(j)*img1.G.n+(k)];
         temp.B.ar[((j) - Ymin)*temp.B.n+((k+1)-Xmin)] = img1.B.ar[(j)*img1.B.n+(k)];

         temp.R.ar[((j+1) - Ymin)*temp.R.n+((k+1)-Xmin)] = img1.R.ar[(j)*img1.R.n+(k)];
         temp.G.ar[((j+1) - Ymin)*temp.G.n+((k+1)-Xmin)] = img1.G.ar[(j)*img1.G.n+(k)];
         temp.B.ar[((j+1) - Ymin)*temp.B.n+((k+1)-Xmin)] = img1.B.ar[(j)*img1.B.n+(k)];*/
       //----------------------------------------------------------------------------------
        }
        if(j<img2.R.m && k<img2.R.n && j<img1.R.m && k<img1.R.n)
        {
            jaff =((j-temp.dm)-c*(k-temp.dn)-f)/d;
            kaff =((k-temp.dn)-b*(j-temp.dm)-e)/a;

            if(kaff>0 && jaff>0 && kaff<(n) && jaff<(m))
            {
               temp2.R.ar[jaff*temp2.R.n+kaff] = img2.R.ar[(j)*img2.R.n+(k)];
               temp2.G.ar[jaff*temp2.G.n+kaff] = img2.G.ar[(j)*img2.G.n+(k)];
               temp2.B.ar[jaff*temp2.B.n+kaff] = img2.B.ar[(j)*img2.B.n+(k)];
            }
         }

      }


    temp.DiscretCosinusTransformation();
    temp2.DiscretCosinusTransformation();
    for (int j = 0; j < (m-8); j+=8)
    {
       for (int k = 0; k < (n-8); k+=8)
       {
          if(temp.R.ar[(j)*temp.R.n+(k)]== 0 && temp.G.ar[(j)*temp.B.n+(k)] ==0 && temp.B.ar[(j)*temp.B.n+(k)] ==0)
          {
             for(int i=0; i<8; i++)
             for(int ii=0; ii<8; ii++)
             {
                 temp.R.ar[(j+i)*temp.R.n+(k+ii)] = temp2.R.ar[(j+i)*temp2.R.n+(k+ii)];
                 temp.G.ar[(j+i)*temp.G.n+(k+ii)] = temp2.G.ar[(j+i)*temp2.G.n+(k+ii)];
                 temp.B.ar[(j+i)*temp.B.n+(k+ii)] = temp2.B.ar[(j+i)*temp2.B.n+(k+ii)];
             }
          }
          else
          {
             if(temp2.R.ar[j*temp2.R.n+k] != 0 || temp2.G.ar[j*temp2.G.n+k] != 0 || temp2.B.ar[j*temp2.B.n+k] != 0)
             {
                for(int i=0; i<8; i++)
                for(int ii=0; ii<8; ii++)
                {  /*
                   if(temp.G.ar[(j+i)*temp.G.n+(k+ii)]>temp2.G.ar[(j+i)*temp2.G.n+(k+ii)])
                   {
                     temp.R.ar[(j+i)*temp.R.n+(k+ii)] = temp.R.ar[(j+i)*temp.R.n+(k+ii)]*0.9 + temp2.R.ar[(j+i)*temp2.R.n+(k+ii)]*0.1;
                     temp.G.ar[(j+i)*temp.G.n+(k+ii)] = temp.G.ar[(j+i)*temp.G.n+(k+ii)]*0.9 + temp2.G.ar[(j+i)*temp2.G.n+(k+ii)]*0.1;
                     temp.B.ar[(j+i)*temp.B.n+(k+ii)] = temp.B.ar[(j+i)*temp.B.n+(k+ii)]*0.9 + temp2.B.ar[(j+i)*temp2.B.n+(k+ii)]*0.1;
                   }
                   if(temp.G.ar[(j+i)*temp.G.n+(k+ii)]<temp2.G.ar[(j+i)*temp2.G.n+(k+ii)])
                   {
                     temp.R.ar[(j+i)*temp.R.n+(k+ii)] = temp.R.ar[(j+i)*temp.R.n+(k+ii)]*0.1 + temp2.R.ar[(j+i)*temp2.R.n+(k+ii)]*0.9;
                     temp.G.ar[(j+i)*temp.G.n+(k+ii)] = temp.G.ar[(j+i)*temp.G.n+(k+ii)]*0.1 + temp2.G.ar[(j+i)*temp2.G.n+(k+ii)]*0.9;
                     temp.B.ar[(j+i)*temp.B.n+(k+ii)] = temp.B.ar[(j+i)*temp.B.n+(k+ii)]*0.1 + temp2.B.ar[(j+i)*temp2.B.n+(k+ii)]*0.9;
                   }       */
                   temp.R.ar[(j+i)*temp.R.n+(k+ii)] = temp.R.ar[(j+i)*temp.R.n+(k+ii)]*0.5 + temp2.R.ar[(j+i)*temp2.R.n+(k+ii)]*0.5;
                   temp.G.ar[(j+i)*temp.G.n+(k+ii)] = temp.G.ar[(j+i)*temp.G.n+(k+ii)]*0.5 + temp2.G.ar[(j+i)*temp2.G.n+(k+ii)]*0.5;
                   temp.B.ar[(j+i)*temp.B.n+(k+ii)] = temp.B.ar[(j+i)*temp.B.n+(k+ii)]*0.5 + temp2.B.ar[(j+i)*temp2.B.n+(k+ii)]*0.5;
                }
             }
          }
        }
     }
      temp.InverseDiscretCosinusTransformation();
   return temp;
}


#pragma package(smart_init)
