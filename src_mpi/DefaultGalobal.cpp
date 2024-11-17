#include <stdio.h>
#include "global.h"


void DefaultGlobal() {
	int TaskType[NumberOfTaskTypes];

	for (int i=0;i<NumberOfTaskTypes; i++) {
		TaskType[i] = i;
	}

}
