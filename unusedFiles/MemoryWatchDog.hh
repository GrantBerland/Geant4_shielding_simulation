

#ifndef memorywatchdog_h
#define memorywatchdog_h 1

#include "sys/types.h"
#include "sys/sysinfo.h"
#include <iostream>


class MemoryWatchDog
{
  public:
	  MemoryWatchDog();
	  ~MemoryWatchDog();
	  void GetTotalMemoryAvailable();
	  double GetMemoryUsage();
	  void KillProgram();
  private:
	  long long fSystemMemory;
	  double BtoGB;
	  double BtoMB;
	  double BtokB;
	  struct sysinfo memInfo;
};



#endif
