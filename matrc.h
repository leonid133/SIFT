//---------------------------------------------------------------------------

#ifndef matrcH
#define matrcH

//#define TYPE_MATRIC float
#include <math.h>
#include <fstream.h>
#include <iostream.h>
#include <stdlib.h>
#include <vcl.h>

#include "Unit1.h"
 //---------------------------------------------------------------------------
class vektor{
public:
   float z, x ;                //z - ��� ������, x - ��� �������
	vektor (const vektor& v): z (v.z), x (v.x)  {};
	vektor (float vz = 0 ,float vx = 0 ): z(vz), x(vx){};

	vektor& operator = (const vektor& v) {z = v.z; x = v.x;  return *this;};
   void FromPolar(float R, float alfa_r)
   {
      z = R*sin(alfa_r);
      x = R*cos(alfa_r);
   };

	vektor  operator - () const	{ return vektor(-z,-x);};
	friend vektor operator + (const vektor& u, const vektor& v);//�������� ��������
	friend vektor operator - (const vektor&, const vektor&);//��������� ��������

	friend float operator  * (const vektor& u, const vektor& v)
			  { return u.z*v.z+u.x*v.x; };//��������� ������������.
	friend vektor operator * (const vektor&, float t);//��������� ������� ��
	friend vektor operator * (float c, const vektor&);//    ������.
	friend vektor operator / (const vektor&, float t);//������� ������� �� ������.

	float operator ! ()const { return (float)sqrt(z*z+x*x); }//������  (�����)
   float alfa(void)const
   {
      float alfa;
      if (fabs(x) < 1e-5)
      {
         if (z > 0)   alfa = M_PI/2;
         else         alfa = -M_PI/2;
      }
      else
         alfa = atan2(z, x);
      return alfa;
   };
	friend float cos(const vektor&, const vektor&);
	friend float acos(const vektor& v1, const vektor& v2)
      {
         double res = (double)cos(v1, v2);
         return (acos(res));
      };
	friend float atan(const vektor& v1, const vektor& v2);
	String show(void);
};

//----------------------------------------------------------------------------
template <class TYPE_MATRIC>
class matric{
public:
   TYPE_MATRIC *ar;
   long m, n;
   ~matric(){delete[] ar;}
   matric(const long m = 0, const long n = 0);
   matric(const matric &r);
   const matric &operator=(const matric &r);

   virtual matric operator+(const matric &r);
   virtual matric operator-(const matric &r);
   virtual matric operator*(const matric &r);
   virtual matric operator*(const TYPE_MATRIC x);
   virtual TYPE_MATRIC D(void);  // �����������
   virtual matric T(void);       // ����������������
   virtual TYPE_MATRIC Sp(void); // ���� ������� 
   String matric::Show(void);
};

//�������� ������� ������-----------------------------------------------------
matric<float> FilterGauss (float sigma, int n, float factorY);

//----------------------------------------------------------------------------
class image{
   long m, n;
public:

   matric<float> B, G, R, Y, Lnorm, Sx, Sy, Sxy, H;
   ~image(void){};
   image(void): m(0){};
   image(const fstream &file);
   const image &operator=(const image &img);

   virtual image LoadFromFile(ifstream &file);
   virtual image Convolution(const matric<float> &Filter); //�������
   virtual image Eqalization(void); //������������ ��c�������� �� ���� ��������
   virtual image DiscretCosinusTransformation(void);
   virtual image InverseDiscretCosinusTransformation(void);
   virtual image FastFuryeTransformation(void);
   virtual image HarrissLaplass(void);
   void Draw(void);
   void DrawY(void);
   int dm, dn;
};

// Discret Cosinus Transformation*--------------------------------------------
matric<float> DCT(matric<float> &r, int M, int N, int NN);

// Inverse Discret Cosinus Transformation*--------------------------------------------
matric<float> IDCT(matric<float> &r, int M, int N, int NN);

//Fast Furye Transformation --------------------------------------------------
matric<int> FFT(matric<float> &r, int M, int N);

//������� � �������-----------------------------------------------------------
long double X(matric<float> &r, int M, int N, float R);

//��� ��������----------------------------------------------------------------
long double Mmn(matric<float> &r, int M, int N, float R);

//������ ������---------------------------------------------------------------
long double DXmn(matric<float> &r, int M, int N, float R);

//���������-------------------------------------------------------------------
long double regres(matric<float> &r1, matric<float> &r2, int r1M, int r1N, int r2M, int r2N, float R);

//������� ����������----------------------------------------------------------
long double Jxy(matric<float> &r1, matric<float> &r2, int r1M, int r1N, int r2M, int r2N, float R);

//����������� �������� � �������� ������� ���������---------------------------
matric<float> MoveOxy(image &rimg, image &rimg2, int ColorFlag);

//�������� �������������� ----------------------------------------------------
void AffinTransformation(image &img, matric<float> &affin);

//�������� �������������� 3D--------------------------------------------------
//void AffinTransformation3D(image &img, matric<float> &affin);

//���������� ������-----------------------------------------------------------
image Glue(image &img1, image &img2, matric<float> &affin);

//---------------------------------------------------------------------------
#endif
