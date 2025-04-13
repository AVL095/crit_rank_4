#include <string.h>
#include "stat.h"

//#####################################################################################
int main() {
 string st;
 ifstream inp("stat.inp");
 inp>>st;
 inp.close();
 cfun[st](st);
}

