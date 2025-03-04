#ifndef FUSEE
#define FUSEE

//----------------------------------------------------------------------
// Bibliothèques utilisées (assert.h représentait un intérêt pour le débogage)
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>   // Nécessite le rajout de -lm en option de compilation
#include <time.h>
#include <string.h>

//----------------------------------------------------------------------
// Constantes utiles (macros)

#define _PI 3.14159
#define _R 8.315
#define _RHO_EAU 998

//----------------------------------------------------------------------
// Structures utilisées

struct rocket_t
{
    double empty_volume;
    double water_volume;
    double radius;
    double nozzle_radius;
    double empty_mass;
    double drag_coeff;
    long pressure;
    double surface;

    // Volumes in m³ ; Mass in kg ; Radius in m ; Pressure in Pa

};
typedef struct rocket_t fusee;

struct bouteille_t
{
    double empty_volume;
    double radius;
    double nozzle_radius;
    double empty_mass;

    // Volumes in m³ ; Mass in kg ; Radius in m

};
typedef struct bouteille_t bouteille;

struct environement_t
{
    double g;
    double rho;
    double T_ext;
    long P_ext;
    double gamma;

    // g in m.s^-2 ; rho in kg.m^-3 ; T_ext in K ; P_ext in Pa

};
typedef struct environement_t env;

//----------------------------------------------------------------------
// Fonctions "interfichiers"

double imp(fusee r, env e);

void expulsion_eau_v1(double T[3], fusee r, env e);

void methode1(double T[3], fusee r, env e);

double methode2(fusee r, env e);

double calcul_vt(fusee r, env e);

double calcul_v_i(fusee r, env e);

double calcul_t_top(double v_i, double v_t, env e, double t_i);

double calcul_z_top(double z_i,double v_t, double v_i, env e);

double calcul_z_top_fusee(fusee r, env e);

void moy_methode(double T[3], fusee r, env e);

fusee* retrouve_fusee_z(env e, double z_max, bouteille b);

double calcul_t_end_fusee(fusee r, env e);

fusee* fusee_max_t(env e, bouteille b);


#endif