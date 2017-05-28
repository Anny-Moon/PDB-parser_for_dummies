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

/*-------------------------------------------------------------------------------
* Parser for ProteinDataBank (.pdb) files.
* PDB: http://www.rcsb.org
* 
* The program reads pdb-file, finds C-alpha (CA) atoms in the first (or only one)
* model, extracts x-, y-, z-coordinates for them and writes the result into
* dat-file in format:
* <x_firts_atom>    <y_first_atom>    <z_first_atom>
*  ...
* <x_last_atom>     <y_last_atom>     <z_last_atom>
*
* If there are any missing atoms in the model, the program will tell about it
* and show maps for them both in percents and actual size.
*
* You will be able to rewrite dat-file with only one segment in order to not have any
* missings.
*
* If there are any CA atoms with the same number, the program automatically
* takes only the first and ignore all the rest with the same number.
* 
* (see Readme.me for more details and examples)
* 
* Anna Sinelnikova
* Uppsala, 2017
*------------------------------------------------------------------------------*/

#include "CA.h"
#include <curl/curl.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


struct Chunk{
    char* data;
    int size;
};

size_t pca_curl_write_callback(char *data, size_t size, size_t nmemb, void *buffer){
    unsigned int i;
    char* newData;
    char* tmp;
    int buff_size;
    
    newData = (char*)malloc(nmemb);
    
    for(i=0; i<nmemb; i++){
	newData[i] = data[i];
    }
    
    buff_size = ((struct Chunk*)buffer)->size;
    tmp = (char*)malloc(buff_size);
    
    for(i=0; i<buff_size; i++){
	tmp[i] =((char*)(((struct Chunk*)buffer)->data ))[i];
//	printf("%c",tmp[i]);
    }
    
    ((struct Chunk*)buffer)->data = (char*)realloc(((struct Chunk*)buffer)->data, (buff_size + nmemb));
    for(i=0; i<buff_size; i++)
	(((struct Chunk*)buffer)->data)[i] = tmp[i];

    for(i=0; i<nmemb; i++)
	(((struct Chunk*)buffer)->data)[buff_size + i] = newData[i];
    
    ((struct Chunk*)buffer)->size = buff_size+nmemb;
    free(tmp);
    free(newData);
    
    return size*nmemb;
}

int main(int np, char** p)
{
    int k;
    int i;
    int lineCounter = 0;
    int firstAtom, lastAtom;
    int N, numberOfAtom = -1;
    
    char str [100];
    
    char end []="END";
    char endmdl []="ENDMDL";
    char atom [] ="ATOM";
    char ter [] = "TER";
    char ca [] ="CA";
    
    int missingFrom[100];
    int missingTo[100];
    int missingCounter = 0;
    
    char *fname;
    char *url;
    FILE *fp, *fp1;

    if(p[1]==NULL){
	printf("\nI need the name of pdb-file without extention as an argument.\nExample: ./pdf_reader 5dn7\n");
	printf("Then I will open 'data/<name_you_gave_me>.pdb'.\n");
	printf("I will find there all CA-atoms coordinates and write them in 'results/xyz_<name_you_gave_me>.dat'.\n");
	printf("In addition I will show you maps of missing atoms if there are any.\n\n");
	
	printf("In case of missing atoms in the model:");
	printf("If you want to rewrite dat-file with only one segment of the chain\n");
	printf("between two missing parts, or between the beginning/end of the chain and\n");
	printf("the closest missing atom, then call the program again and give number\n");
	printf("of the first atom(or any missing atom before the first) of your segment as the second argument.\n");
	printf("Example: ./pdb_reader 5dn7 366\n");
	printf("The same: ./pdb_reader 5dn7 365\n");
	printf("For the first segment (from the beggining to the first missing atom)\n");
	printf("You can always pass the second argument as 0.\n\n");
	exit(1);
    }

    url = (char*) malloc(1000);
    sprintf(url,"https://files.rcsb.org/view/%s.pdb",p[1]);
    printf("%s\n",url);
    
//    char** dataFromWeb;
//    dataFromWeb =(char**)malloc(sizeof(char*));
//    dataFromWeb[0] =(char*)malloc(5);
    
    struct Chunk* dataFromWeb;
    dataFromWeb = (struct Chunk*)malloc(sizeof(struct Chunk));
    dataFromWeb->size = 1;
    dataFromWeb->data = (char*)malloc(dataFromWeb->size);
    
    CURL *easy_handle = curl_easy_init();
    curl_easy_setopt(easy_handle, CURLOPT_URL, url);
    curl_easy_setopt(easy_handle, CURLOPT_WRITEFUNCTION, pca_curl_write_callback);
    curl_easy_setopt(easy_handle, CURLOPT_WRITEDATA, dataFromWeb);
    curl_easy_perform(easy_handle);
    
    for(i=0;i<dataFromWeb->size;i++)
	printf("%c",dataFromWeb->data[i]);
    
    fname = (char*) malloc(100);
    sprintf(fname,"data/%s.pdb",p[1]);
    fp=fopen(fname,"r");
    if(fp==NULL){
	printf("Error: Can't find %s.pdb\nNote: write the name of the file without extantion.\n",p[1]);
	exit(1);
    }
    free(fname);
    
    printf("\n");
    printf("Protein: %s\n",p[1]);
    
    fname = (char*) malloc(100);
    sprintf(fname,"results/xyz_%s.dat",p[1]);
    fp1=fopen(fname, "w");
    if(fp1==NULL){
	printf("Error: Can't open '%s'\nNote: probably you don't have folde 'results'.\n",fname);
	exit(1);
    }

// reading the whole model and presenting of map of missimg atoms
    if(p[2]==NULL){
	do{
	    fscanf(fp,"%s",str);
	    if(pca_compare_strings(str,atom)){ //found ATOM
		fscanf(fp,"%s",str); //trash
		fscanf(fp,"%s",str);
	    
		if(pca_compare_strings(str,ca)){ //found CA -c-alpha atom
		    fscanf(fp,"%s",str); //trash
		    fscanf(fp,"%1s",str); //trash
		
		    fscanf(fp,"%i",&k);
		    //printf("%i\n",k);
		    if(k!=numberOfAtom+1 && lineCounter!=0){
		    
			if(k==numberOfAtom)
			    continue;
		    
			else if(k>numberOfAtom+1){
			    missingFrom[missingCounter] = numberOfAtom+1;
			    missingTo[missingCounter] = k-1;
			    printf("Missing atoms from %i to %i (%i atoms).\n",\
				missingFrom[missingCounter], missingTo[missingCounter], k-numberOfAtom-1);
			    missingCounter++;
			}
			else{
			    printf("Error: strange order of atoms numbers before atom %i\n Exit\n",k);
			    exit(1);
			}
		    
		    }
		
		    numberOfAtom = k;
		    if(lineCounter==0){
			firstAtom = numberOfAtom;
			printf("The first CA atom has number %i.\n",firstAtom);
		    }
		    
		    fscanf(fp,"%s",str);
		    fprintf(fp1,"%s\t",str);//x
		    fscanf(fp,"%s",str);
		    fprintf(fp1,"%s\t",str);//y
		    fscanf(fp,"%s",str);
		    fprintf(fp1,"%s\n",str);//z
		
		    lineCounter++;
		}
	    }
	}
	while(!pca_compare_strings(str,end) && !pca_compare_strings(str,endmdl)&& !pca_compare_strings(str,ter));//after this loop str=etalon
	lastAtom = numberOfAtom;
    
	N = lastAtom-firstAtom+1;
	
	printf("The last CA atom has number  %i.\n",lastAtom);
	
	printf("Number of CA atoms in the model: %i.\n",N);
	
	if(missingCounter!=0){
	    printf("But there is data only for %i of them.\n", lineCounter);
	    printf("\n");
	    printf("***************************************\n");
	    printf("*      Maps of missing atoms          *\n");
	    printf("*  atom: .         missing atom: 0    *\n");
	    printf("***************************************\n");
	    printf("Percentage: string length = 100 chars.\n\n");
	    pca_map_of_missings(firstAtom, lastAtom, missingCounter, missingFrom, missingTo, 100);
	    printf("\n");
	    printf("Actual: string length = number of atoms in model.\n\n");
	    pca_map_of_missings(firstAtom, lastAtom, missingCounter, missingFrom, missingTo, N);
	    printf("\n");
	    printf("*****************************************\n");
	    printf("If you want to rewrite dat-file with only one segment, call\n./pdb-reader (without arguments) for instructions.\n");
	    printf("\n");
	}
	
	else
	    printf("No missing data in pdb-file.\n\n");
	
	
    }

// reading from atom with given number to the first miss atom
    else{
	firstAtom = atoi(p[2]);
	numberOfAtom = firstAtom-1;
	do{
	    fscanf(fp,"%s",str);
	    if(pca_compare_strings(str,atom)){ //found ATOM
		fscanf(fp,"%s",str); //trash
		fscanf(fp,"%s",str);
	    
		if(pca_compare_strings(str,ca)){ //found CA -c-alpha atom
		    fscanf(fp,"%s",str); //trash
		    fscanf(fp,"%1s",str); //trash
		
		    fscanf(fp,"%i",&k);
		    //printf("%i\n",k);
		    if(k>=firstAtom){
		    
			if(k==numberOfAtom)
			    continue;
		    
			else if(k>numberOfAtom+1 && lineCounter!=0){
			    break;
			}
			
			else{
			    fscanf(fp,"%s",str);
			    fprintf(fp1,"%s\t",str);//x
			    fscanf(fp,"%s",str);
			    fprintf(fp1,"%s\t",str);//y
			    fscanf(fp,"%s",str);
			    fprintf(fp1,"%s\n",str);//z
			}
		    
		    lineCounter++;
		    }
		
		    numberOfAtom = k;
		}
	    }
	}
	while(!pca_compare_strings(str,end) && !pca_compare_strings(str,endmdl)&& !pca_compare_strings(str,ter));//after this loop str=etalon
	
	printf("Number of CA atoms: %i.\n", lineCounter);
    
    }
    
    fclose(fp1);
    fclose(fp);
    free(fname);
    return 0;
}
