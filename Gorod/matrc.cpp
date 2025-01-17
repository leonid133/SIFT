//---------------------------------------------------------------------------

#pragma hdrstop

#include "matrc.h"

//---------------------------------------------------------------------------
//class Vektor****************************************************************
//----------------------------------------------------------------------------
const float ToGrad = 57.2957795130823;

//���������� �������� +  ---------------------------------------------------
vektor operator + (const vektor& u,const vektor& v)
{
	return vektor (u.z+v.z, u.x+v.x);
}

//���������� �������� -  ---------------------------------------------------
vektor operator - (const vektor& u,const vektor& v)
{
	return vektor (u.z-v.z, u.x-v.x);
}

//���������� �������� *. ��������� ������� �� ������.------------------------
vektor operator* (const vektor& v, float f)
{
	return vektor (f*v.z, f*v.x);
}

//���������� �������� *. ��������� ������� �� ������.------------------------
vektor operator* (float f, const vektor& v)
{
	return vektor (f*v.z, f*v.x);
}

//������� ������� �� ������.-------------------------------------------------
vektor operator / (const vektor& v, float f)
{
	return vektor (v.z/f, v.x/f);
}

//������� ���� ����� ���������.---------------------------------------------
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

//�������� ���� ����� ���������.----------------------------------------------
float atan(const vektor& a, const vektor& b)
{
   return a.alfa() - b.alfa();
}

//class Matric ***************************************************************
//----------------------------------------------------------------------------
template <class TYPE_MATRIC> matric<TYPE_MATRIC>::matric(int m, int n)
   : m(m)
   , n(n)
{
   ar = new TYPE_MATRIC [m*n];
};

//----------------------------------------------------------------------------
template <class TYPE_MATRIC> matric<TYPE_MATRIC>::matric(const matric &r)
{
   m = r.m;
   n = r.n;
   ar = new TYPE_MATRIC [m*n];

   for (int j = 0; j < m*n; j++)
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

      for (int j = 0; j < m*n; j++)
         ar[j] =  r.ar[j];
   }
   return *this;
};

//��������--------------------------------------------------------------------
/*virtual*/template <class TYPE_MATRIC>
matric<TYPE_MATRIC> matric<TYPE_MATRIC>::operator+(const matric &r)
{
   matric temp(m, n);

   for (int j = 0; j < m*n; j++)
      temp.ar[j] =  ar[j]+r.ar[j];
   return temp;
};

//��������� ------------------------------------------------------------------
/*virtual*/template <class TYPE_MATRIC>
matric<TYPE_MATRIC> matric<TYPE_MATRIC>::operator-(const matric &r)
{
   matric temp(m, n);
   for (int j = 0; j < m*n; j++)
      temp.ar[j] =  ar[j]-r.ar[j];

   return temp;
};

//���������-------------------------------------------------------------------
/*virtual*/template <class TYPE_MATRIC>
matric<TYPE_MATRIC> matric<TYPE_MATRIC>::operator*(const matric &r)
{
   matric temp(m, n);
   for (int j = 0; j < m; j++)
      for (int k = 0; k < n; k++)
      {
         temp.ar[j*n+k] = 0;
         for (int i = 0; i < n; i++)
            temp.ar[j*n+k] += ar[j*n+i] * r.ar[i*n+k];
      }
   return temp;
};

//��������� �� ������----------------------------------------------------------
/*virtual*/template <class TYPE_MATRIC>
matric<TYPE_MATRIC> matric<TYPE_MATRIC>::operator*(const TYPE_MATRIC x)
{
   matric temp(m, n);
   for (int j = 0; j < m*n; j++)
      temp.ar[j] =  ar[j]*x;

   return temp;
};

//������������----------------------------------------------------------------
/*virtual*/template <class TYPE_MATRIC>
TYPE_MATRIC  matric<TYPE_MATRIC>::D(void)
{
   TYPE_MATRIC temp = 1, temp2 = -1, sum = 0;
   int i, i2;
   for(int kk = 1; kk < n+1; kk++)
   {
      i = kk - 1 ;
      i2 = n - kk;
      for (int j = 0; j < m; j++)
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

//���������������� ������-----------------------------------------------------
/*virtual*/template <class TYPE_MATRIC> matric<TYPE_MATRIC> matric<TYPE_MATRIC>::T(void)
{
   matric temp(n, m);
   for(int j = 0; j < n; j++)
      for(int k = 0; k < m; k++)
         temp.ar[j*n+k] = ar[k*m+j];

   return temp;
};

//���� ������� ---------------------------------------------------------------
/*virtual*/template <class TYPE_MATRIC> TYPE_MATRIC matric<TYPE_MATRIC>::Sp(void)
{
    TYPE_MATRIC Spur=0;
    int x=m<n?m:n;
    for(int i=0; i<x; i++)
    Spur+=ar[i*m+i];
    return Spur;
};
//�����-----------------------------------------------------------------------
String matric<float>::Show(void)
{
   String tmp(""), txt;
   for (int j = 0; j < m; j++)
      for (int k = 0; k < n; k++)
      {
         tmp = tmp.FormatFloat("0.00     ", ar[j*n+k]);
         txt += tmp;
      }
   return txt;
};
//�����-----------------------------------------------------------------------
String matric<int>::Show(void)
{
   String tmp(""), txt;
   for (int j = 0; j < m; j++)
      for (int k = 0; k < n; k++)
      {
         tmp = ar[j*n+k];
         txt += tmp;
      }
   return txt;
};
//������ ������---------------------------------------------------------------
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
   }
   return *this;
};

//������������� �� ������-----------------------------------------------------
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
   B = G = R = Y = matric<float>(m, n);
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

//�������� ��������������-----------------------------------------------------
/*virtual*/image image :: AffinTransformation(matric<float> &affin)
{
   m = G.m;
   n = G.n;
   image temp = *this;
   float a, b, c, d, e, f, jaff, kaff;
   a = affin.ar[0*n+0];
   b = affin.ar[0*n+1];
   e = affin.ar[0*n+2];
   c = affin.ar[1*n+0];
   d = affin.ar[1*n+1];
   f = affin.ar[1*n+2];
   for (int j = 0; j < m-1; j++)
      for (int k = 0; k < n-1; k++)
      {
         kaff = a*k + b*j + e;
         jaff = c*k + d*j + f;
         if(kaff >= n)
         {
            kaff = n-1;
            R.ar[j*R.n+(int)kaff] = 255;
            G.ar[j*G.n+(int)kaff] = 249;
            B.ar[j*B.n+(int)kaff] = 0;
         }
         if(jaff >= m)
         {
            jaff = m-1;
            R.ar[(int)jaff*R.n+k] = 255;
            G.ar[(int)jaff*G.n+k] = 249;
            B.ar[(int)jaff*B.n+k] = 0;
         }
         if(kaff <= 0)
         {
            kaff = 1;
            R.ar[j*R.n+(int)kaff] = 255;
            G.ar[j*G.n+(int)kaff] = 249;
            B.ar[j*B.n+(int)kaff] = 0;
         }
         if(jaff <= 0)
         {
            jaff = 1;
            R.ar[(int)jaff*R.n+k] = 255;
            G.ar[(int)jaff*G.n+k] = 249;
            B.ar[(int)jaff*B.n+k] = 0;
         }
         temp.R.ar[j*temp.R.n+k] = R.ar[(int)jaff*R.n+(int)kaff];
         temp.G.ar[j*temp.G.n+k] = G.ar[(int)jaff*G.n+(int)kaff];
         temp.B.ar[j*temp.B.n+k] = B.ar[(int)jaff*B.n+(int)kaff];
      }
   return temp;
};

//�������---------------------------------------------------------------------
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
//��������� A = A * Filter � ������ � temp
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

//�����������-----------------------------------------------------------------
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
   matric<int> A(8, 8);
   matric<float> temp(Y);

   for (int m = 0; m < Y.m; m+=8)
      for (int n = 0; n < Y.n; n+=8)
      {
         A = DCT(Y, m, n);
         for(int j=0; j<8; j++)
            for(int k=0; k<8; k++)
            {
               Y.ar[(m+j)*Y.n + (n+k)]=A.ar[j*8+k];
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
   float Etta=1.2, S=0.7, h;

   matric<float> A(2, 2), GaussN(2, 2), GaussD(2, 2);
   int al, bl, cl, dl, el, fl, gl, hl, il, l, lx, ly, lxy;
   matric<float> L, Lx, Ly, Lxx, Lxy, Lxxyy, Lyy, LoG;

   L=Lnorm=Lx=Ly=Lxx=Lxy=Lxxyy=Lyy=Sx=Sy=Sxy=LoG=H=Y*0;

  for(int scale=1; scale<10; scale++)
   {
      GaussN=FilterGauss(pow(Etta, scale), 2, 1);
      GaussD=FilterGauss(S*pow(Etta, scale), 2, 1);

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
                  al = (jm - GaussN.m + jj - 2)*Y.n + (kn - GaussN.n + kk - 2);
                  bl = (jm - GaussN.m + jj - 2)*Y.n + (kn - GaussN.n + kk - 1);
                  cl = (jm - GaussN.m + jj - 2)*Y.n + (kn - GaussN.n + kk);
                  dl = (jm - GaussN.m + jj - 1)*Y.n + (kn - GaussN.n + kk - 2);
                  lxy = el = (jm - GaussN.m + jj - 1)*Y.n + (kn - GaussN.n + kk - 1);
                  ly = fl = (jm - GaussN.m + jj - 1)*Y.n + (kn - GaussN.n + kk);
                  gl = (jm - GaussN.m + jj)*L.n + (kn - GaussN.n + kk - 2);
                  lx = hl = (jm - GaussN.m + jj)*Y.n + (kn - GaussN.n + kk - 1);
                  l = il = (jm - GaussN.m + jj)*Y.n + (kn - GaussN.n + kk);

                  L.ar[l] = 0;
                  for (int i = 0; i < A.n; i++)
                     L.ar[l] += A.ar[jj*A.n+i] * (GaussD.ar[i*GaussN.n+kk] - GaussN.ar[i*GaussN.n+kk]);
                  Lnorm.ar[l] += L.ar[l] * pow(S*pow(Etta, scale), 2);
                 //���������� ��������
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
                  // �������� �����
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


//�������� �����--------------------------------------------------------------
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
matric<int> DCT(matric<float> &r, int M, int N)
{
   matric<int> DCTkvant(8, 8);
   matric<int> kvant(8, 8);
   float Cu, Cv, AC;
   for (int u = 0; u < 8; u++)
   {
      for (int v = 0; v < 8; v++)
      {
         if(u == 0){Cu=1/pow(2, 0.5);}
         else{Cu = 1;}
         if(v == 0){Cv=1/pow(2, 0.5);}
         else{Cv = 1;}
         AC=0;
         for (int x = 0; x < 8; x++)
         {
            for (int y = 0; y < 8; y++)
            {
               {
                  AC += r.ar[(y+M)*8+x+N]*cos((2*x+1)*u*3.14/16)*cos((2*y+1)*v*3.14/16);
               }
            }
         }
         DCTkvant.ar[v*8+u]=0.25*Cu*Cv*AC;
      }
   }
   DCTkvant.ar[0]=DCTkvant.ar[0]/8;
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
   DCTkvant.ar[63]=DCTkvant.ar[63]/99;

   
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

//������� ��������------------------------------------------------------------
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

//��� �������� ---------------------------------------------------------------
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

//������ ������ ---------------------------------------------------------------
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


//��������� -------------------------------------------------------------------
long double regres(matric<float> &r1, matric<float> &r2, int r1M, int r1N, int r2M, int r2N, float R)
{
   long double Regres = 0;

   Regres = pow(pow(DXmn(r1, r1M, r1N, R)-DXmn(r2, r2M, r2N, R), 2), 0.5);

   return Regres;
}

//������� ����������-----------------------------------------------------------
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

//����������� ������������� ��������� ��������������--------------------------
matric<float> MoveOxy(image &rimg, image &rimg2, int ColorFlag)
{
   int okresnost=14, R=25, Segment=3 , alfa=7, scale=1, Points=200;
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

   goodPoints[0] = 0;
   bool f = false;
   matric<float> delta(3, 3), deltaA(3, 3), deltaB(3, 3), deltaC(3, 3), deltaD(3, 3), deltaE(3, 3), deltaF(3, 3);
   matric<float> delta2(3, 3), deltaA2(3, 3), deltaB2(3, 3), deltaC2(3, 3), deltaD2(3, 3), deltaE2(3, 3), deltaF2(3, 3);

   X1 = (okresnost+R+alfa); //if((float)r.n/2-200>X1){X1 = (int)(float)r.n/2-200;}   //����� ����� ���������� ���� ������������ ��������� �����������
   Y1 = (okresnost+R+alfa); //if((float)r.m/2-300>Y1){Y1 = (int)(float)r.m/2-300;}   //����� ����� ���������� ���� ������������ ��������� �����������
   X2 = rimg.Y.n-(okresnost+R+alfa); //if((float)r.n/2+200<X2){X2 = (int)(float)r.n/2+200;}  //����� ����� ���������� ���� ������������ ��������� �����������
   Y2 = rimg.Y.m-(okresnost+R+alfa); //if((float)r.m/2+300<Y2){Y2 = (int)(float)r.m/2+300;}  //����� ����� ���������� ���� ������������ ��������� �����������
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
  skewness2= new float[KeyPoint];
  neighbor_contrast2= new float[KeyPoint];
  beta_coeff2= new float[KeyPoint];

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

      std_dev[i] = sqrt(std_dev[i]/SumN);
      skewness[i] /= SumN;
      skewness[i] = skewness[i]/(std_dev[i]*std_dev[i]*std_dev[i]);
      neighbor_contrast[i] /= SumN;
      if (sxx < 1e-6) {beta_coeff[i] = 0.0;}
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
//*****************************************************************************
for (int uu = -R; uu < R+1; uu++)
{
	for (int vv = -R; vv < R+1; vv++)
	{
		for (int i = 0; i < KeyPoint; i++)
		{
			sxx=sxy=x=diff=mean_rad=SumN=col_dist=row_dist_sq=center_row=center_col=0;
			X1 = (int)v1[i].x  - Segment + vv;
			Y1 = (int)v1[i].z  - Segment + uu;
			X2 = (int)v1[i].x  + Segment + vv;
			Y2 = (int)v1[i].z  + Segment + uu;
			SumN = ((Y2+1)*rimg2.Y.n+(X2+1) - Y1*rimg2.Y.n+X1);
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

				std_dev2[i] = sqrt(std_dev2[i]/SumN);
				skewness2[i] /= SumN;
				skewness2[i] = skewness2[i]/(std_dev2[i]*std_dev2[i]*std_dev2[i]);
				neighbor_contrast2[i] /= SumN;
				if (sxx < 1e-6) {beta_coeff2[i] = 0.0;}
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
		for (int i = 0; i < KeyPoint; i++)
		{
		  mean2[i] = (mean2[i]-meanMIN)/(meanMAX-meanMIN);
		  std_dev2[i] = (std_dev2[i]-std_devMIN)/(std_devMAX-std_devMIN);
		  skewness2[i] = (skewness2[i]-skewnessMIN)/(skewnessMAX-skewnessMIN);
		  neighbor_contrast2[i] = (neighbor_contrast2[i]-neighbor_contrastMIN)/(neighbor_contrastMAX-neighbor_contrastMIN);
		  beta_coeff2[i] = (beta_coeff2[i]-beta_coeffMIN)/(beta_coeffMAX-beta_coeffMIN);
		}
		dSUM=0;
		for (int i = 0; i < KeyPoint; i++)
		{
		   d[i] =  (float)abs((int)((mean[i]-mean2[i])*100))/100+(float)abs((int)((std_dev[i]-std_dev2[i])*100))/100+(float)abs((int)((skewness[i]-skewness2[i])*100))/100+(float)abs((int)((neighbor_contrast[i]-neighbor_contrast2[i])*100))/100+(float)abs((int)((beta_coeff[i]-beta_coeff2[i])*100))/100;
		   if(i==0)
		   {
			  dMIN=dMAX=d[i];
		   }
		   if(dMIN>d[i])dMIN=d[i];
		   if(dMAX<d[i])dMAX=d[i];
		   dSUM+=d[i];
		}
		if(dSUM<dTOL)
		{
         dTOL=dSUM;
			goodPoints[0]=0;
			for (int i = 0; i < KeyPoint; i++)
			{
			   d[i]=(d[i]-dMIN)/(dMAX-dMIN);
			   if(d[i]<0.2)
			   {
				  v3[i].x = v1[i].x + vv;
				  v3[i].z = v1[i].z + uu;
              v2=v3;
              v[i] = v[i] + (v2[i] - v1[i]);
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
            for (int i = 0; i < KeyPoint; i++)
            {
				   X1 = (((float)aa)/100)*(int)v3[i].x + (((float)bb)/100)*(int)v3[i].z - Segment;
				   Y1 = (((float)cc)/100)*(int)v3[i].z + (((float)ga)/100)*(int)v3[i].x - Segment;
				   X2 = (((float)aa)/100)*(int)v3[i].x + (((float)bb)/100)*(int)v3[i].z + Segment;
				   Y2 = (((float)cc)/100)*(int)v3[i].z + (((float)ga)/100)*(int)v3[i].x + Segment;

               if(X1>rimg2.Y.n)X1=rimg2.Y.n;
               if(Y1>rimg2.Y.m)Y1=rimg2.Y.m;
               if(X1<0)X1=0;
               if(Y1<0)Y1=0;
               if(X2>rimg2.Y.n)X2=rimg2.Y.n;
               if(Y2>rimg2.Y.m)Y2=rimg2.Y.m;
               if(X2<0)X2=0;
               if(Y2<0)Y2=0;

				   SumN = ((Y2+1)*rimg2.Y.n+(X2+1) - Y1*rimg2.Y.n+X1);
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
					  std_dev2[i] = sqrt(std_dev2[i]/SumN);
					  skewness2[i] /= SumN;
					  if(std_dev2[i]*std_dev2[i]*std_dev2[i] < 1e-6){skewness2[i] = 0.0;}
                 else{skewness2[i] = skewness2[i]/(std_dev2[i]*std_dev2[i]*std_dev2[i]);}
					  neighbor_contrast2[i] /= SumN;
					  if (sxx < 1e-6) {beta_coeff2[i] = 0.0;}
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
				for (int i = 0; i < KeyPoint; i++)
				{
				  mean2[i] = (mean2[i]-meanMIN)/(meanMAX-meanMIN);
				  std_dev2[i] = (std_dev2[i]-std_devMIN)/(std_devMAX-std_devMIN);
				  skewness2[i] = (skewness2[i]-skewnessMIN)/(skewnessMAX-skewnessMIN);
				  neighbor_contrast2[i] = (neighbor_contrast2[i]-neighbor_contrastMIN)/(neighbor_contrastMAX-neighbor_contrastMIN);
				  beta_coeff2[i] = (beta_coeff2[i]-beta_coeffMIN)/(beta_coeffMAX-beta_coeffMIN);
				}
				dSUM=0;
				for (int i = 0; i < KeyPoint; i++)
				{
					d[i] =  (float)abs((int)((mean[i]-mean2[i])*100))/100+(float)abs((int)((std_dev[i]-std_dev2[i])*100))/100+(float)abs((int)((skewness[i]-skewness2[i])*100))/100+(float)abs((int)((neighbor_contrast[i]-neighbor_contrast2[i])*100))/100+(float)abs((int)((beta_coeff[i]-beta_coeff2[i])*100))/100;
					if(i==0){dMIN=dMAX=d[i];}
					if(dMIN>d[i])dMIN=d[i];
					if(dMAX<d[i])dMAX=d[i];
					dSUM+=d[i];
				}
				if(dSUM<=dTOL)
				{
               dTOL=dSUM;
				   goodPoints[0]=0;
					for (int i = 0; i < KeyPoint; i++)
					{
						d[i]=(d[i]-dMIN)/(dMAX-dMIN);
						if(d[i]<0.2 && v3[i].x>1 && v3[i].z>1)
						{
						   v2[i].x = (((float)aa)/100)*(int)v3[i].x + (((float)bb)/100)*(int)v3[i].z;
						   v2[i].z = (((float)cc)/100)*(int)v3[i].z + (((float)ga)/100)*(int)v3[i].x;
						   v[i] = v[i] + (v2[i] - v1[i]);
                     goodPoints[0]++;
						   goodPoints[goodPoints[0]] = i;

						}
      			}
				}
	      }
      }
   }
}

//******************************************************************************

  Form1->Memo1->Lines->Add((goodPoints[0]-2));
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

      //����� ���������� ���� ������
      Form1->Image2->Canvas->Refresh();
      Form1->Image2->Canvas->MoveTo(v2[goodPoints[i]].x, v2[goodPoints[i]].z);

      switch (ColorFlag){
      case 0: Form1->Image2->Canvas->Pen->Color = clRed; break;
      case 1: Form1->Image2->Canvas->Pen->Color = clGreen; break;
      case 2: Form1->Image2->Canvas->Pen->Color = clBlue;  break;
      case 3: Form1->Image2->Canvas->Pen->Color = clYellow; break;  }
      Form1->Image2->Canvas->LineTo(v1[goodPoints[i]].x, v1[goodPoints[i]].z);

   }
   AffinTransformation.ar[0*AffinTransformation.n+0] = deltaA.D()/delta.D();
   AffinTransformation.ar[0*AffinTransformation.n+1] = deltaB.D()/delta.D();
   AffinTransformation.ar[1*AffinTransformation.n+0] = deltaC.D()/delta.D();
   AffinTransformation.ar[1*AffinTransformation.n+1] = deltaD.D()/delta.D();
   AffinTransformation.ar[0*AffinTransformation.n+2] = deltaE.D()/delta.D();
   AffinTransformation.ar[1*AffinTransformation.n+2] = deltaF.D()/delta.D();
   return AffinTransformation;
}

//���������� ������-----------------------------------------------------------
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
   x[5]=(-b*(img2.R.m)-e)/a-img1.dn;
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
   temp2.B = temp2.G = temp2.R = matric<float>(2*m, 2*n)*0;
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
       
         temp2.R.ar[(j*2 - Ymin*2)*temp2.R.n+(k*2-Xmin*2)] = img1.R.ar[(j)*img1.R.n+(k)];
         temp2.G.ar[(j*2 - Ymin*2)*temp2.G.n+(k*2-Xmin*2)] = img1.G.ar[(j)*img1.G.n+(k)];
         temp2.B.ar[(j*2 - Ymin*2)*temp2.B.n+(k*2-Xmin*2)] = img1.B.ar[(j)*img1.B.n+(k)];
        }
        if(j<img2.R.m && k<img2.R.n && j<img1.R.m && k<img1.R.n)
        {
         jaff =((j-temp.dm)-c*(k-temp.dn)-f)/d;
         kaff =((k-temp.dn)-b*(j-temp.dm)-e)/a;

         if(kaff>0 && jaff>0 && kaff<n && jaff<m)
         {

          //-------------------------------------------------------------------------

               temp.R.ar[(jaff)*temp.R.n+(kaff)] = img2.R.ar[(j)*img2.R.n+(k)];
               temp.G.ar[(jaff)*temp.G.n+(kaff)] = img2.G.ar[(j)*img2.G.n+(k)];
               temp.B.ar[(jaff)*temp.B.n+(kaff)] = img2.B.ar[(j)*img2.B.n+(k)];

            if(temp2.R.ar[jaff*temp2.R.n+kaff] == 0)
            {
               temp2.R.ar[2*jaff*temp2.R.n+2*kaff] = img2.R.ar[(j)*img2.R.n+(k)];
               temp2.G.ar[2*jaff*temp2.G.n+2*kaff] = img2.G.ar[(j)*img2.G.n+(k)];
               temp2.B.ar[2*jaff*temp2.B.n+2*kaff] = img2.B.ar[(j)*img2.B.n+(k)];

               temp2.R.ar[(2*jaff-1)*temp2.R.n+2*kaff] = img2.R.ar[(j)*img2.R.n+(k)];
               temp2.G.ar[(2*jaff-1)*temp2.G.n+2*kaff] = img2.G.ar[(j)*img2.G.n+(k)];
               temp2.B.ar[(2*jaff-1)*temp2.B.n+2*kaff] = img2.B.ar[(j)*img2.B.n+(k)];

               temp2.R.ar[2*jaff*temp2.R.n+2*kaff-1] = img2.R.ar[(j)*img2.R.n+(k)];
               temp2.G.ar[2*jaff*temp2.G.n+2*kaff-1] = img2.G.ar[(j)*img2.G.n+(k)];
               temp2.B.ar[2*jaff*temp2.B.n+2*kaff-1] = img2.B.ar[(j)*img2.B.n+(k)];

               temp2.R.ar[(2*jaff-1)*temp2.R.n+2*kaff-1] = img2.R.ar[(j)*img2.R.n+(k)];
               temp2.G.ar[(2*jaff-1)*temp2.G.n+2*kaff-1] = img2.G.ar[(j)*img2.G.n+(k)];
               temp2.B.ar[(2*jaff-1)*temp2.B.n+2*kaff-1] = img2.B.ar[(j)*img2.B.n+(k)];
            }
            else
            {
             /*  if(img2.R.ar[(j)*img2.R.n+(k)]!=0 || img2.G.ar[(j)*img2.G.n+(k)]!=0 || img2.B.ar[(j)*img2.B.n+(k)])
               {
                  temp2.R.ar[2*jaff*temp2.R.n+2*kaff] = img2.R.ar[(j)*img2.R.n+(k)]*0.5 +temp2.R.ar[2*jaff*temp2.R.n+2*kaff]*0.5;
                  temp2.G.ar[2*jaff*temp2.G.n+2*kaff] = img2.G.ar[(j)*img2.G.n+(k)]*0.5+temp2.G.ar[2*jaff*temp2.G.n+2*kaff]*0.5;
                  temp2.B.ar[2*jaff*temp2.B.n+2*kaff] = img2.B.ar[(j)*img2.B.n+(k)]*0.5+temp2.B.ar[2*jaff*temp2.B.n+2*kaff]*0.5;

                  temp2.R.ar[(2*jaff-1)*temp2.R.n+2*kaff] = img2.R.ar[(j)*img2.R.n+(k)]*0.5+temp2.R.ar[2*jaff*temp2.R.n+2*kaff]*0.5;
                  temp2.G.ar[(2*jaff-1)*temp2.G.n+2*kaff] = img2.G.ar[(j)*img2.G.n+(k)]*0.5+temp2.G.ar[2*jaff*temp2.G.n+2*kaff]*0.5;
                  temp2.B.ar[(2*jaff-1)*temp2.B.n+2*kaff] = img2.B.ar[(j)*img2.B.n+(k)]*0.5+temp2.B.ar[2*jaff*temp2.B.n+2*kaff]*0.5;

                  temp2.R.ar[2*jaff*temp2.R.n+2*kaff-1] = img2.R.ar[(j)*img2.R.n+(k)]*0.5+temp2.R.ar[2*jaff*temp2.R.n+2*kaff]*0.5;
                  temp2.G.ar[2*jaff*temp2.G.n+2*kaff-1] = img2.G.ar[(j)*img2.G.n+(k)]*0.5+temp2.G.ar[2*jaff*temp2.G.n+2*kaff]*0.5;
                  temp2.B.ar[2*jaff*temp2.B.n+2*kaff-1] = img2.B.ar[(j)*img2.B.n+(k)]*0.5+temp2.B.ar[2*jaff*temp2.B.n+2*kaff]*0.5;

                  temp2.R.ar[(2*jaff-1)*temp2.R.n+2*kaff-1] = img2.R.ar[(j)*img2.R.n+(k)]*0.5+temp2.R.ar[2*jaff*temp2.R.n+2*kaff]*0.5;
                  temp2.G.ar[(2*jaff-1)*temp2.G.n+2*kaff-1] = img2.G.ar[(j)*img2.G.n+(k)]*0.5+temp2.G.ar[2*jaff*temp2.G.n+2*kaff]*0.5;
                  temp2.B.ar[(2*jaff-1)*temp2.B.n+2*kaff-1] = img2.B.ar[(j)*img2.B.n+(k)]*0.5+temp2.B.ar[2*jaff*temp2.B.n+2*kaff]*0.5;
               }
               else  */
               {
                  temp2.R.ar[2*jaff*temp2.R.n+2*kaff] = img2.R.ar[(j)*img2.R.n+(k)];
                  temp2.G.ar[2*jaff*temp2.G.n+2*kaff] = img2.G.ar[(j)*img2.G.n+(k)];
                  temp2.B.ar[2*jaff*temp2.B.n+2*kaff] = img2.B.ar[(j)*img2.B.n+(k)];

                  temp2.R.ar[(2*jaff-1)*temp2.R.n+2*kaff] = img2.R.ar[(j)*img2.R.n+(k)];
                  temp2.G.ar[(2*jaff-1)*temp2.G.n+2*kaff] = img2.G.ar[(j)*img2.G.n+(k)];
                  temp2.B.ar[(2*jaff-1)*temp2.B.n+2*kaff] = img2.B.ar[(j)*img2.B.n+(k)];

                  temp2.R.ar[2*jaff*temp2.R.n+2*kaff-1] = img2.R.ar[(j)*img2.R.n+(k)];
                  temp2.G.ar[2*jaff*temp2.G.n+2*kaff-1] = img2.G.ar[(j)*img2.G.n+(k)];
                  temp2.B.ar[2*jaff*temp2.B.n+2*kaff-1] = img2.B.ar[(j)*img2.B.n+(k)];

                  temp2.R.ar[(2*jaff-1)*temp2.R.n+2*kaff-1] = img2.R.ar[(j)*img2.R.n+(k)];
                  temp2.G.ar[(2*jaff-1)*temp2.G.n+2*kaff-1] = img2.G.ar[(j)*img2.G.n+(k)];
                  temp2.B.ar[(2*jaff-1)*temp2.B.n+2*kaff-1] = img2.B.ar[(j)*img2.B.n+(k)];
               }

            }

            //------------------------------------------------------------------

         }
        }
      }
     for (int j = 2; j < (2*m); j++)
     {
      for (int k = 2; k < (2*n); k++)
         {
            if(temp.R.ar[(j/2)*temp.R.n+(k/2)]== 0 && temp.G.ar[(j/2)*temp.B.n+(k/2)] ==0 && temp.B.ar[(j/2)*temp.B.n+(k/2)] ==0)
            {
                  temp.R.ar[(j/2)*temp.R.n+(k/2)] = temp2.R.ar[j*temp2.R.n+k];
                  temp.G.ar[(j/2)*temp.G.n+(k/2)] = temp2.G.ar[j*temp2.G.n+k];
                  temp.B.ar[(j/2)*temp.B.n+(k/2)] = temp2.B.ar[j*temp2.B.n+k];
                /*
                  temp.R.ar[(j/2)*temp.R.n+(k/2)] = temp2.R.ar[j*temp2.R.n+k]*0.25
                  +temp2.R.ar[(j-1)*temp2.R.n+k]*0.25
                  +temp2.R.ar[(j-1)*temp2.R.n+k-1]*0.25
                  +temp2.R.ar[(j)*temp2.R.n+k-1]*0.25;

                  temp.G.ar[(j/2)*temp.G.n+(k/2)] = temp2.G.ar[j*temp2.G.n+k]*0.25
                  +temp2.G.ar[(j-1)*temp2.G.n+k]*0.25
                  +temp2.G.ar[(j-1)*temp2.G.n+k-1]*0.25
                  +temp2.G.ar[(j)*temp2.G.n+k-1]*0.25;

                  temp.B.ar[(j/2)*temp.B.n+(k/2)] = temp2.B.ar[j*temp2.B.n+k]*0.25
                  +temp2.B.ar[(j-1)*temp2.B.n+k]*0.25
                  +temp2.B.ar[(j-1)*temp2.B.n+k-1]*0.25
                  +temp2.B.ar[(j)*temp2.B.n+k-1]*0.25;         */
            }
            else
            {
               if(temp2.R.ar[j*temp2.R.n+k] != 0 || temp2.G.ar[j*temp2.G.n+k] != 0 || temp2.B.ar[j*temp2.B.n+k] != 0)
               {
                  temp.R.ar[(j/2)*temp.R.n+(k/2)] = temp.R.ar[(j/2)*temp.R.n+(k/2)]*0.5 + temp2.R.ar[j*temp2.R.n+k]*0.5;
                  temp.G.ar[(j/2)*temp.G.n+(k/2)] = temp.G.ar[(j/2)*temp.G.n+(k/2)]*0.5 + temp2.G.ar[j*temp2.G.n+k]*0.5;
                  temp.B.ar[(j/2)*temp.B.n+(k/2)] = temp.B.ar[(j/2)*temp.B.n+(k/2)]*0.5 + temp2.B.ar[j*temp2.B.n+k]*0.5;
                 /*
                  temp.R.ar[(j/2)*temp.R.n+(k/2)] =temp.R.ar[(j/2)*temp.R.n+(k/2)]*0.5
                  +temp2.R.ar[j*temp2.R.n+k]*0.25/2
                  +temp2.R.ar[(j-1)*temp2.R.n+k]*0.25/2
                  +temp2.R.ar[(j-1)*temp2.R.n+k-1]*0.25/2
                  +temp2.R.ar[(j)*temp2.R.n+k-1]*0.25/2;

                  temp.G.ar[(j/2)*temp.G.n+(k/2)] =temp.G.ar[(j/2)*temp.G.n+(k/2)]*0.5
                  +temp2.G.ar[j*temp2.G.n+k]*0.25/2
                  +temp2.G.ar[(j-1)*temp2.G.n+k]*0.25/2
                  +temp2.G.ar[(j-1)*temp2.G.n+k-1]*0.25/2
                  +temp2.G.ar[(j)*temp2.G.n+k-1]*0.25/2;

                  temp.B.ar[(j/2)*temp.B.n+(k/2)] =temp.B.ar[(j/2)*temp.B.n+(k/2)]*0.5
                  +temp2.B.ar[j*temp2.B.n+k]*0.25/2
                  +temp2.B.ar[(j-1)*temp2.B.n+k]*0.25/2
                  +temp2.B.ar[(j-1)*temp2.B.n+k-1]*0.25/2
                  +temp2.B.ar[(j)*temp2.B.n+k-1]*0.25/2;   */
               }
            }
         }
      }
   return temp;
}


#pragma package(smart_init)
