//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "Unit1.h"
#include "matrc.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;
//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
   : TForm(Owner)
{
}
//---------------------------------------------------------------------------


void __fastcall TForm1::Button1Click(TObject *Sender)
{
   Form1->Memo1->Lines->Add(GetTickCount());
   image F[6];
   ifstream img1("1.bmp", ios::binary);
   if (!img1)
   {
      throw Exception ("�� ���� ������� ���� 1");
   }
   F[1].LoadFromFile(img1);
   img1.close();

   ifstream img2("2.bmp", ios::binary);
   if (!img2)
   {
      throw Exception ("�� ���� ������� ���� 2");
   }
   F[2].LoadFromFile(img2);
   img2.close();

   ifstream img3("3.bmp", ios::binary);
   if (!img3)
   {
      throw Exception ("�� ���� ������� ���� 3");
   }
   F[3].LoadFromFile(img3);
   img3.close();

   ifstream img4("4.bmp", ios::binary);
   if (!img4)
   {
      throw Exception ("�� ���� ������� ���� 4");
   }
   F[4].LoadFromFile(img4);
   img4.close();

   ifstream img5("5.bmp", ios::binary);
   if (!img5)
   {
      throw Exception ("�� ���� ������� ���� 5");
   }
   F[5].LoadFromFile(img5);
   img5.close();


   matric<float> affin[5];
   affin[0].m=3; affin[0].n=2;
   for(int i=0; i<6 ;i++) 
      affin[0].ar[i]=0;
   int ColorFlag =0;
   for(int i=1; i<5; i++)
     {
        affin[i].m=3; affin[i].n=2;
        F[i].HarrissLaplass();
        affin[i] = MoveOxy(F[i], F[i+1], ColorFlag);
        ColorFlag++;
        if(ColorFlag>3)ColorFlag=0;
        Form1->Memo1->Lines->Add(affin[i].Show());

     }

   for(int i=1; i<5; i++)
     {
        F[1] = Glue(F[1], F[i+1], affin[1]);
        affin[1].ar[0] = affin[1].ar[0] * affin[i+1].ar[0] + affin[1].ar[1] * affin[i+1].ar[3];
        affin[1].ar[1] = affin[1].ar[0] * affin[i+1].ar[1] + affin[1].ar[1] * affin[i+1].ar[4];
        affin[1].ar[2] = affin[1].ar[0] * affin[i+1].ar[2] + affin[1].ar[1] * affin[i+1].ar[5] + affin[1].ar[2];
        affin[1].ar[3] = affin[1].ar[3] * affin[i+1].ar[0] + affin[1].ar[4] * affin[i+1].ar[3];
        affin[1].ar[4] = affin[1].ar[3] * affin[i+1].ar[1] + affin[1].ar[4] * affin[i+1].ar[4];
        affin[1].ar[5] = affin[1].ar[3] * affin[i+1].ar[2] + affin[1].ar[4] * affin[i+1].ar[5] + affin[1].ar[5];
        Form1->Memo1->Lines->Add(affin[1].Show());
     }
     F[1].Draw();
   Form1->Memo1->Lines->Add(GetTickCount());
}

//-----------------------------------------------------------------------------


