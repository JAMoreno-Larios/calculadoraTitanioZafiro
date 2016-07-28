#include "parser.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <cmath>
#include <QDebug>


parser::parser(std::string rutaArchivo)
{

    this->archivo.open(rutaArchivo, std::ios::out | std::ios::in);
    if(!this->archivo.is_open())
    {
        std::ofstream creaArchivo(rutaArchivo);
        creaArchivo<<std::endl;
        creaArchivo.close();
        this->archivo.open(rutaArchivo, std::ios::out | std::ios::in);
    }

}

parser::~parser()
{
    this->archivo.close();
    std::map<std::string,char*>::iterator iterador;
    for(iterador=this->mapaChar.begin(); iterador!=this->mapaChar.end();iterador++)
    {
        if(iterador->second!=NULL)
        {
            free(iterador->second);
            iterador->second=NULL;
        }
    }
    this->mapaChar.clear();
    this->mapaBool.clear();
    this->mapaDobles.clear();
    this->mapaEnteros.clear();
}

void parser::cargaVariables(void)
{
    if(this->archivo.is_open())
    {
        this->archivo.seekp(0,std::ios_base::beg);
        this->archivo.seekg(0,std::ios_base::beg);
        while(std::getline(this->archivo,linea))
        {
            std::vector<std::string> datosImportar; // Primer localidad almacena etiqueta; segunda, el valor.
            char primerDelimitador=' ';
            char segundoDelimitador='"';
            std::size_t posPrimerDelimitador=linea.find_first_of(primerDelimitador,0);
            datosImportar.push_back(linea.substr(0,posPrimerDelimitador));
            datosImportar.push_back(linea.substr(posPrimerDelimitador+1, linea.size()-posPrimerDelimitador));
            if(datosImportar[1].find_first_of(segundoDelimitador,0)==0)
            {
                std::size_t posUltimoDelimitador=datosImportar[1].find_first_of(segundoDelimitador,1);
                std::string datos=datosImportar[1].substr(1,posUltimoDelimitador-1);
                if(mapaChar[datosImportar[0]]!=NULL)
                    free(mapaChar[datosImportar[0]]);
                mapaChar[datosImportar[0]]=convierteACharPtr(datos);
               // printf("Etiqueta: %s, valor: %s\n",datosImportar[0].c_str(),mapaChar[datosImportar[0]]);
            }
            else if(datosImportar[1]=="true"|datosImportar[1]=="false")
            {
                bool valorBooleanoo = convierteABool(datosImportar[1]);
                mapaBool[datosImportar[0]]=valorBooleanoo;
                // Carga a memoria
                //std::cout<<"Etiqueta: " << datosImportar[0] << ", valor " << mapaBool[datosImportar[0]] << "\n";
            }
            else if (std::string::npos!=datosImportar[1].find('e')| datosImportar[1]=="INFINITY" | datosImportar[1]=="INF" |datosImportar[1]=="-INFINITY"|datosImportar[1]=="-INF"|datosImportar[1]=="NaN")
            {
                double valorDoble = convierteADoble(datosImportar[1]);
                mapaDobles[datosImportar[0]]=valorDoble;
                //std::cout<<"Etiqueta: " << datosImportar[0] << ", valor " << mapaDobles[datosImportar[0]] << "\n";
                // Carga a memoria.
            }
            else if (!datosImportar[1].empty() && datosImportar[1].find_first_not_of("+-0123456789") == std::string::npos)
            {
                int valorInt = convierteAEntero(datosImportar[1]);
                mapaEnteros[datosImportar[0]]=valorInt;
                //std::cout<<"Etiqueta: " << datosImportar[0] << ", valor " << mapaEnteros[datosImportar[0]] << "\n";
                // Carga a memoria.
            }
        }
    }
}

double parser::convierteADoble(const std::string& s)
{
    double x;
    if(s=="INFINITY"|s=="INF"|s=="-INFINITY"|s=="-INF"|s=="NaN")
    {
        if(s=="INFINITY"|s=="INF")
            x=INFINITY;
        if(s=="-INFINITY"|s=="-INF")
            x=-INFINITY;
        if(s=="NaN")
            x=NAN;
        return x;
    }
    else
    {
        std::vector<std::string> numero; // Primer localidad almacena etiqueta; segunda, el valor.
        char delimitador='e';
        std::size_t posDelimitador=s.find_first_of(delimitador,0);
        numero.push_back(s.substr(0,posDelimitador));
        numero.push_back(s.substr(posDelimitador+1, s.size()-posDelimitador));
        std::istringstream numeroBase(numero[0]), exponente(numero[1]);
        double parte1, x;
        int parte2;
        numeroBase>>parte1;
        exponente>>parte2;
        x=parte1*pow(10,parte2);
        return x;
    }
}

void parser::escribeDatos(std::map<std::string, char*> mapaChar, std::map<std::string, bool> mapaBool,\
                           std::map<std::string,int>mapaInt, std::map<std::string,double>mapaDobles)
                           {
                               if(this->archivo.is_open())
                               {
                                    this->archivo.clear();
                                    this->archivo.seekp(0,std::ios_base::beg);
                                    this->archivo.seekg(0,std::ios_base::beg);
                                    std::map<std::string,char*>::iterator iteradorChar;
                                    std::map<std::string,bool>::iterator iteradorBool;
                                    std::map<std::string,int>::iterator iteradorInt;
                                    std::map<std::string,double>::iterator iteradorDoble;
                                    for(iteradorChar=mapaChar.begin(); iteradorChar!=mapaChar.end();iteradorChar++)
                                    {
                                        if(iteradorChar->second!=NULL)
                                        {
                                            this->archivo << iteradorChar->first << " \"" << iteradorChar->second << "\"" << std::endl;
                                            //std::cout << iteradorChar->first << " \"" << iteradorChar->second << "\"" << std::endl;
                                        }
                                    }
                                    for(iteradorBool=mapaBool.begin(); iteradorBool!=mapaBool.end();iteradorBool++)
                                    {
                                        if(iteradorBool->second)
                                            this->archivo << iteradorBool->first << " " << "true" << std::endl;
                                        else
                                            this->archivo << iteradorBool->first << " " << "false" << std::endl;
                                        //std::cout << iteradorBool->first << " " << iteradorBool->second << std::endl;
                                    }
                                    for(iteradorInt=mapaInt.begin(); iteradorInt!=mapaInt.end();iteradorInt++)
                                    {
                                        this->archivo << iteradorInt->first << " " << iteradorInt->second << std::endl;
                                        //std::cout << iteradorInt->first << " " << iteradorInt->second << std::endl;
                                    }
                                    for(iteradorDoble=mapaDobles.begin(); iteradorDoble!=mapaDobles.end();iteradorDoble++)
                                    {
                                        this->archivo << iteradorDoble->first << " " << std::scientific << std::setprecision(15) << iteradorDoble->second << std::endl;
                                        //std::cout << iteradorDoble->first << " " << std::scientific << std::setprecision(15) << iteradorDoble->second << std::endl;
                                    }
                               }
                           }
