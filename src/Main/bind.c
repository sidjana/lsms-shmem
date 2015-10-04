#include<unistd.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>

void bindmyproc()
{

  int mypid=getpid();
  char str[100];
  printf("my pid is %d",mypid);
  sprintf(str,"/bin/echo %d >> /dev/cpuset/socket0/tasks ",mypid);
  system(str);

}


void bindmyproc_()
{
 bindmyproc();
}
