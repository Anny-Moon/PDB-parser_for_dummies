/*  Copyright 2017 Anna Sinelnikova
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*/

#include "CA.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

bool pca_compare_strings(char *a, char *b)
{int i=0;
	do{
	    if(a[i]==b[i])
		i++;
	
	    else
		return false;
	}
	while(b[i]!='\0'||a[i]!='\0');
return true;
}

int pca_rounding(double number)
{
    double intpart;
    double fracpart;
    int answ;

    fracpart=modf (number , &intpart);

    if(fabs(fracpart)<0.5000)
	answ=(int)intpart;
    
    else
	answ=(int)intpart+1;
    
    return answ;
}

void pca_map_of_missings(int firstAtom, int lastAtom, int missingCounter, const int* missingFrom, const int* missingTo, int stringLength)
{
    int i,j, tmpN;
    double unit;
    
    unit = (double)stringLength/(lastAtom - firstAtom + 1);
    
    tmpN = pca_rounding((double)(missingFrom[0] - firstAtom)*unit);
    for(i=0;i<tmpN;i++)
	printf(".");
	
    for(j=0;j<missingCounter-1;j++){
	tmpN = pca_rounding((double)(missingTo[j]-missingFrom[j]+1)*unit);
	for(i=0;i<tmpN;i++)
	    printf("0");
	    
	tmpN = pca_rounding((double)(missingFrom[j+1] - missingTo[j]-1)*unit);
	for(i=0;i<tmpN;i++)
	    printf(".");
    }
    
    tmpN = pca_rounding((double)(missingTo[missingCounter-1] - missingFrom[missingCounter-1]+1)*unit);
    for(i=0;i<tmpN;i++)
	printf("0");
    
    tmpN = pca_rounding((double)(lastAtom - missingTo[missingCounter-1])*unit);
    for(i=0;i<tmpN;i++)
	printf(".");
	
    printf("\n");
}