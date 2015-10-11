//---------------------------------------------------------------------------

#include "Unit1.h"
#include "matrc.h"

#pragma hdrstop


//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;
image F[100];
matric<float> affinS(3, 2);


//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
   : TForm(Owner)
{
}
//---------------------------------------------------------------------------


void __fastcall TForm1::Button1Click(TObject *Sender)
{
   affinS.ar[0*3+0] = 1;
   affinS.ar[0*3+1] = 0;
   affinS.ar[0*3+2] = 0;
   affinS.ar[1*3+0] = 0;
   affinS.ar[1*3+1] = 1;
   affinS.ar[1*3+2] = 0;
   int timeStart = GetTickCount(), file_collection = 0;

   char file[7];
   for (int i = 1; i <= 100; i++)
   {
      sprintf(file, "%d.bmp", i);
      ifstream img(file, ios::binary);
      if (img)
      {
         StatusBar1->SetTextBuf(file);
         F[i].LoadFromFile(img);
         img.close();
         file_collection++;
      }
      else
      {
         img.close();
         break;
      }
   }

   matric<float> *affin;
   affin = new matric<float>[file_collection+1];
   for (int i = 0; i <= (file_collection+1); i++)
      affin[i] = matric<float>(3, 2);

   int ColorFlag = 0;
   Memo1->Lines->Add("Прошло времени:");
   Memo1->Lines->Add((GetTickCount() - timeStart)/1000);

   for(int i = 1; i < file_collection; i++)
   {
      StatusBar1->SimpleText = "Определение функции Харриссона-Лапласса изображения № " + (String)i;
      F[i].HarrissLaplass();
      StatusBar1->SimpleText = "Определение смещения изображений № " + (String)i + " и " + (String)(i+1);
      affin[i] = MoveOxy(F[i], F[i+1], ColorFlag);
      ColorFlag++;
      if(ColorFlag > 3)
         ColorFlag = 0;
      Memo1->Lines->Add("Аффинные коэффициенты смещения изображений " + (String)i +  " и " + (String)(i+1));
      Memo1->Lines->Add(affin[i].Show());
   }
   Memo1->Lines->Add("Прошло времени:");
   Memo1->Lines->Add((GetTickCount() - timeStart)/1000);
   for(int i = 1; i < (file_collection); i++)
   {
     StatusBar1->SimpleText = "Обработка изображения № " + (String)(i+1);
     F[1] = Glue(F[1], F[i+1], affin[1]);
     Memo1->Lines->Add("Глобальные аффинные координаты смещения изображений " + (String)i +  " и " + (String)(i+1));
     Memo1->Lines->Add(affin[1].Show());
     affin[1].ar[0] = affin[1].ar[0] * affin[i+1].ar[0] + affin[1].ar[1] * affin[i+1].ar[3];
     affin[1].ar[1] = affin[1].ar[0] * affin[i+1].ar[1] + affin[1].ar[1] * affin[i+1].ar[4];
     affin[1].ar[2] = affin[1].ar[0] * affin[i+1].ar[2] + affin[1].ar[1] * affin[i+1].ar[5] + affin[1].ar[2];
     affin[1].ar[3] = affin[1].ar[3] * affin[i+1].ar[0] + affin[1].ar[4] * affin[i+1].ar[3];
     affin[1].ar[4] = affin[1].ar[3] * affin[i+1].ar[1] + affin[1].ar[4] * affin[i+1].ar[4];
     affin[1].ar[5] = affin[1].ar[3] * affin[i+1].ar[2] + affin[1].ar[4] * affin[i+1].ar[5] + affin[1].ar[5];

   }
   StatusBar1->SimpleText = "Вывод на экран склеенной матрицы.";
   F[1].Draw();

   Memo1->Lines->Add("Прошло времени:");
   Memo1->Lines->Add((GetTickCount() - timeStart)/1000);
   StatusBar1->SimpleText = "***************************";
   delete [] affin;
}

//-----------------------------------------------------------------------------

void __fastcall TForm1::SpeedButton1Click(TObject *Sender)
{
    Form1->StatusBar1->SimpleText = "Увеличение масштаба";

    affinS.ar[0*3+0] -= 0.10;
    affinS.ar[1*3+1] -= 0.10;

    AffinTransformation(F[1], affinS);

    Form1->StatusBar1->SimpleText = "***************************";
}
//---------------------------------------------------------------------------

void __fastcall TForm1::SpeedButton2Click(TObject *Sender)
{
   Form1->StatusBar1->SimpleText = "Уменьшение масштаба";

    affinS.ar[0*3+0] += 0.10;
    affinS.ar[1*3+1] += 0.10;

   AffinTransformation(F[1], affinS);
   Form1->StatusBar1->SimpleText = "***************************";
}
//---------------------------------------------------------------------------
#pragma package(smart_init)
void __fastcall TForm1::SpeedButton5Click(TObject *Sender)
{
   Form1->StatusBar1->SimpleText = "Перемещение";

    affinS.ar[0*3+2] -= 20;

   AffinTransformation(F[1], affinS);
   Form1->StatusBar1->SimpleText = "***************************";
}
//---------------------------------------------------------------------------

void __fastcall TForm1::SpeedButton6Click(TObject *Sender)
{
   Form1->StatusBar1->SimpleText = "Перемещение";

    affinS.ar[0*3+2] += 20;

   AffinTransformation(F[1], affinS);
   Form1->StatusBar1->SimpleText = "***************************";
}
//---------------------------------------------------------------------------

void __fastcall TForm1::SpeedButton4Click(TObject *Sender)
{
   Form1->StatusBar1->SimpleText = "Перемещение";

    affinS.ar[1*3+2] += 20;

   AffinTransformation(F[1], affinS);
   Form1->StatusBar1->SimpleText = "***************************";
}
//---------------------------------------------------------------------------

void __fastcall TForm1::SpeedButton3Click(TObject *Sender)
{
   Form1->StatusBar1->SimpleText = "Перемещение";

    affinS.ar[1*3+2] -= 20;

   AffinTransformation(F[1], affinS);
   Form1->StatusBar1->SimpleText = "***************************";
}
//---------------------------------------------------------------------------

void __fastcall TForm1::SpeedButton7Click(TObject *Sender)
{
   Form1->StatusBar1->SimpleText = "Сброс масштаба";
   affinS.ar[0*3+0] = 1;
   affinS.ar[0*3+1] = 0;
   affinS.ar[0*3+2] = 0;
   affinS.ar[1*3+0] = 0;
   affinS.ar[1*3+1] = 1;
   affinS.ar[1*3+2] = 0;
   AffinTransformation(F[1], affinS);
   Form1->StatusBar1->SimpleText = "***************************";   
}
//---------------------------------------------------------------------------

