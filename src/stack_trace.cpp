/*
 * stack_trace.cpp
 *
 *  Created on: Jul 28, 2015
 *      Author: dmarce1
 */




/* SOURCE : http://stackoverflow.com/questions/77005/how-to-generate-a-stacktrace-when-my-gcc-c-app-crashes */


#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>


void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}


__attribute__((constructor))
void install_stack_trace(){
#ifndef NDEBUG
	 signal(SIGSEGV, handler);   // install our handler
	 signal(SIGABRT, handler);   // install our handler
	 signal(SIGFPE, handler);   // install our handler
	 signal(SIGILL, handler);   // install our handler
#endif
}
