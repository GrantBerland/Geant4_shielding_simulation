


#include "sys/types.h"
#include "sys/sysinfo.h"
#include <stdexcept>
#include "MemoryWatchDog.hh"

// Constructor
MemoryWatchDog::MemoryWatchDog()
: BtoGB(1000000000.),
  BtoMB(1000000.),
  BtokB(1000.)
{
  sysinfo (&memInfo);

}


MemoryWatchDog::~MemoryWatchDog()
{

}

void MemoryWatchDog::GetTotalMemoryAvailable()
{
  fSystemMemory = memInfo.totalram;
  fSystemMemory += memInfo.totalswap;
  fSystemMemory *= memInfo.mem_unit;
}

double MemoryWatchDog::GetMemoryUsage()
{

 long long virtualMemUsed = memInfo.totalram - memInfo.freeram; 
 virtualMemUsed += memInfo.totalswap - memInfo.freeswap;
 virtualMemUsed *= memInfo.mem_unit;
 
 return virtualMemUsed/BtoGB;

}

void MemoryWatchDog::KillProgram(){
  std::cout 
	  << "Program exiting due to too much memory usage!" 
	  << std::endl;
  
  std::exit(0);
}

