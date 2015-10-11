//---------------------------------------------------------------------------

#ifndef Unit1H
#define Unit1H
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include <ComCtrls.hpp>
#include <Buttons.hpp>




//---------------------------------------------------------------------------
class TForm1 : public TForm
{
__published:	// IDE-managed Components
   TButton *Button1;
   TImage *Image2;
   TMemo *Memo1;
   TStatusBar *StatusBar1;
   TImage *Image1;
   TSpeedButton *SpeedButton1;
   TSpeedButton *SpeedButton2;
   TSpeedButton *SpeedButton3;
   TSpeedButton *SpeedButton4;
   TSpeedButton *SpeedButton5;
   TSpeedButton *SpeedButton6;
   TSpeedButton *SpeedButton7;
   void __fastcall Button1Click(TObject *Sender);
   void __fastcall SpeedButton1Click(TObject *Sender);
   void __fastcall SpeedButton2Click(TObject *Sender);
   void __fastcall SpeedButton3Click(TObject *Sender);
   void __fastcall SpeedButton4Click(TObject *Sender);
   void __fastcall SpeedButton5Click(TObject *Sender);
   void __fastcall SpeedButton6Click(TObject *Sender);
   void __fastcall SpeedButton7Click(TObject *Sender);
private:	// User declarations
public:		// User declarations
   
   __fastcall TForm1(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;
//---------------------------------------------------------------------------
#endif
