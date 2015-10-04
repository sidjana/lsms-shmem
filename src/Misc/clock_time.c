#include <time.h>
#include "clock.h"

double clock_time_()
{
  return get_rtc();
//    return(((double)clock())/((double)CLOCKS_PER_SEC));
}
