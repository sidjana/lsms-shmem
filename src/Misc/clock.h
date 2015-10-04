


 
 unsigned long long get_rtc(void) {
     unsigned long long rtc;
 
     do {
 	asm volatile ("rdtsc\n\t"
 		      "salq $32, %%rdx\n\t"
 		      "orq %%rdx, %%rax\n\t"
 		      "movq %%rax, %0"
 		      : "=r" (rtc)
 		      : /* no inputs */
 		      : "%rax", "%rdx");
     } while (0);
 
     return rtc;
 }
