// Guardará los datos generados por cada caso.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>


/* guardaVariaciones: Guarda el tamaño de spot contra número de iteraciones
 * ütil para verificar convergencia de propuesta.
 * */
void guardaVariaciones(long double * spots, int * iteraciones, int iterMax, long double epsilon1, long double epsilon2, char * prefijo, char * ruta,  _Bool termino)
{
    char fecha[12];
    time_t rawtime;
    struct tm * timeinfo;
    time(&rawtime);
    timeinfo=localtime(&rawtime);
    strftime(fecha,30,"%F_%H-%M-%S",timeinfo);
    char *bandera;
    if(termino==1)
    {
        bandera=malloc(3);
        bandera[0]='s';
        bandera[1]='i';
        bandera[2]=0;
    }
    else
    {
        bandera=malloc(3);
        bandera[0]='n';
        bandera[1]='o';
        bandera[2]=0;
    }

    int caracteres=snprintf(NULL,0,"%s/%s_E1_%.2Le_E2_%.2Le_%s_%s.csv",ruta,prefijo,epsilon1,epsilon2,fecha,bandera);
    char * buffer=malloc(caracteres+1);
    sprintf(buffer,"%s/%s_E1_%.2Le_E2_%.2Le_%s_%s.csv",ruta,prefijo,epsilon1,epsilon2,fecha,bandera);
    FILE *archivo;
    archivo = fopen(buffer,"w");
    fprintf(archivo,"Iteración,W\n");
    int i;
    for(i=0;i<iterMax;i++)
    {
        //Imprime el elemento de punto flotante
        fprintf(archivo, "%d,%.15Le\n",iteraciones[i],spots[i]);
    }
    fclose(archivo);
    free(buffer);
    free(bandera);
}

int existeArchivo(const char *nombre)
{
   FILE *fp = fopen (nombre, "r");
   if (fp!=NULL) fclose (fp);
   return (fp!=NULL);
}

/* guardaSpotIteraciones guarda el spot final y el número de iteraciones totales*/
void guardaSpotIteraciones(long double spotTan, long double spotSag, int iterMax, long double epsilon1, long double epsilon2, char * ruta, _Bool termino)
{
    int caracteres=snprintf(NULL,0,"%s",ruta);
    char *nombre=malloc(caracteres+1);
    sprintf(nombre,"%s",ruta);
    FILE *archivo;
    if(termino==0)
    {
        spotTan=0.0;
        spotSag=0.0;
    }
    if(existeArchivo(nombre)==0)
    {
        archivo=fopen(nombre,"w");
        fprintf(archivo,"epsilon1,epsilon2,Spot Tan[m], Spot Sag[m],Iteraciones,\n");
        fprintf(archivo,"%.15Le,%.15Le,%.15Le,%.15Le,%d\n",epsilon1,epsilon2,spotTan,spotSag,iterMax);
    }
    else
    {
        archivo=fopen(nombre,"a");
        fprintf(archivo,"%.15Le,%.15Le,%.15Le,%.15Le,%d\n",epsilon1,epsilon2,spotTan,spotSag,iterMax);
    }
    fclose(archivo);
    free(nombre);
}
