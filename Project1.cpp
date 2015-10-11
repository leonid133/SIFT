//---------------------------------------------------------------------------

#include <vcl.h>
#include "matrc.h"

#pragma hdrstop
//---------------------------------------------------------------------------
USEFORM("Unit1.cpp", Form1);
//---------------------------------------------------------------------------
const char *NamedMutex= "OneOnly";
HANDLE CheckInstance(const char *Name)
{
   HANDLE Mutex = CreateMutex(NULL, true, Name);
   int er = GetLastError();
   if (er) return 0;
   return Mutex;
}
WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{
   HANDLE Mutex = CheckInstance(NamedMutex);
   if (!Mutex)
   {
      ShowMessage("Уже работает.");
      ReleaseMutex(Mutex);
      return 1;
   }
   try
   {
       /*
      HANDLE ProcessHandle, ThreadHandle;
      DWORD ProcessID = GetCurrentProcessId();
      ProcessHandle = OpenProcess(PROCESS_ALL_ACCESS,false,ProcessID);
      SetPriorityClass(ProcessHandle,HIGH_PRIORITY_CLASS);
      ThreadHandle = GetCurrentThread();
      SetThreadPriority(ThreadHandle,THREAD_PRIORITY_TIME_CRITICAL); */

      Application->Initialize();
      Application->CreateForm(__classid(TForm1), &Form1);
      Application->Run();
   }
   catch (Exception &exception)
   {
      Application->ShowException(&exception);
   }
   return 0;
}
//---------------------------------------------------------------------------
