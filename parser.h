#ifndef PARSER_H
#define PARSER_H
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <float.h>
#include <cmath>
#include <map>

class parser
{
    public:
        explicit parser(const std::__cxx11::string rutaArchivo);
        virtual ~parser();
        void cargaVariables(void);
        void escribeDatos(std::map<std::string, char*> mapaChar, std::map<std::string, bool> mapaBool,\
                           std::map<std::string,int>mapaInt, std::map<std::string,double>mapaDobles);

        // Almacenamiento de datos

        // Rutas de archivo
        std::map<std::string, char*> mapaChar={{"RutaEM1CW",NULL},{"RutaEM2CW",NULL},{"RutaEM1ML",NULL},{"RutaEM2ML",NULL},\
                                                            {"RutaIteraciones",NULL},{"RutaIterMax",NULL},{"conjugado_largo",NULL},{"conjugado_corto",NULL}};
        // Cálculos a realizar
        std::map<std::string, bool> mapaBool={{"EM1CW",true},{"EM2CW",true},{"EM1ML",true},{"EM2ML",true},{"termico",true},\
                                                              {"guardaSpotsIter",true},{"guardaIterMax",true}};
        // Mapa de enteros
        std::map<std::string,int>mapaEnteros={{"numeroLaminas",0},{"iteraciones",0},{"NTerm",0}};
        // Mapa de dobles
        std::map<std::string,double>mapaDobles={{"L1",0},{"L2",0},{"L",0},{"n",0},{"f1",0},{"f2",0},{"lambda",0},{"PLaser",0},\
                                                       {"umbral",0},{"ladoX",0},{"ladoY",0},{"alpha",0},{"n2",0},{"nPump",0},{"kth",0},{"chi",0},{"Cp",0},{"rho",0},{"dn_dT",0},\
                                                       {"wPumpTan",0},{"wPumpSag",0},{"divPump",0},{"lambdaPump",0},{"PPump",0},{"RA",0},{"RB",0},{"tLente",0},{"nLente",0},\
                                                       {"inclinacionLente",0},{"tEspejo",0},{"nEspejo",0},{"separacionFuenteALente",0},{"separacionLenteAEspejo",0},\
                                               {"epsilon1Min",0},{"epsilon1Max",0},{"pasoEpsilon1",0},{"epsilon2Min",0},{"epsilon2Max",0},{"pasoEpsilon2",0},{"numPasoEpsilon1",0},\
                                               {"numPasoEpsilon2",0}};
        inline bool archivoVacio(void){return archivo.peek()==std::fstream::traits_type::eof();}
        inline bool archivoAbierto(void){return archivo.is_open();}

    protected:

    private:
        std::fstream archivo;
        std::string linea;
        double convierteADoble(const std::string& s);
        inline int convierteAEntero(const std::string& s)
        {
          std::istringstream i(s);
          int x;
          i>>x;
            //throw BadConversion("convertToDouble(\"" + s + "\")");
          return x;
        }
        inline bool convierteABool(const std::string& s)
        {
            if(s=="true")
                return true;
            else if (s=="false")
                return false;
        }
        inline char * convierteACharPtr(const std::string& s)
        {
            int caracteres=snprintf(NULL,0,s.c_str());
            char * ptr =(char*)malloc(caracteres+1);
            sprintf(ptr,s.c_str());
            return ptr;
        }

};

#endif // PARSER_H
