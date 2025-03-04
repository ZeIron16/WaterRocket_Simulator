#include "_fusee.h"

// Fonctions pour retrouver les valeurs :

//----------------------------------------------------------------------
// Première "épreuve" de compétition

bool correct_fusee(fusee r, env e, double z_max){
    double z_top = calcul_z_top_fusee(r, e);
    if(z_top - 2 <= z_max && z_top + 2 >= z_max){return true;}
    return false;
}

void random_rocket(fusee* r){
    r->pressure = (rand() % 500)*1000 + 300000;
    r->water_volume = (rand() % 13) / 10000.0 + 0.0002;
}

fusee* z_top_inverse(env e, double z_max, double vo, double r, double r_n,double m){
    fusee* ro = (fusee*) malloc(sizeof(fusee));
    ro->drag_coeff = 0.5;
    ro->empty_mass = m;
    ro->empty_volume = vo;
    ro->nozzle_radius = r_n;
    ro->pressure = 0;
    ro->radius = r;
    ro->surface = 0;
    ro->water_volume = 0;

    srand(time(NULL));
     while (!correct_fusee(*ro,e, z_max)){
        random_rocket(ro);
    }
    printf("\nres = %lf\n", calcul_z_top_fusee(*ro,e));


    return ro;
}

fusee* retrouve_fusee_z(env e, double z_max, bouteille b){
    fusee* r = z_top_inverse(e, z_max, b.empty_volume, b.radius, b.nozzle_radius, b.empty_mass);
    return r;
}

//----------------------------------------------------------------------
// Deuxième "épreuve" de compétition


fusee* temps_max_inverse(env e, double vo, double r, double r_n,double m){
    fusee* ro = (fusee*) malloc(sizeof(fusee));
    ro->drag_coeff = 0.5;
    ro->empty_mass = m;
    ro->empty_volume = vo;
    ro->nozzle_radius = r_n;
    ro->pressure = 0;
    ro->radius = r;
    ro->surface = ro->radius*_PI*ro->radius;
    ro->water_volume = 0;

    int i = 300, j = 2; // i correspond à 300000 Pa et j à 0.0002 m³ d'eau
    double temps = 0.0;
    fusee* new_ro = malloc(sizeof(fusee));
    
    new_ro->drag_coeff = 0.5;
    new_ro->empty_mass = m;
    new_ro->empty_volume = vo;
    new_ro->nozzle_radius = r_n;
    new_ro->pressure = 0;
    new_ro->radius = r;
    new_ro->surface = new_ro->radius*_PI*new_ro->radius;
    new_ro->water_volume = 0;

    while (i!=801)
    {
        j = 2;
        while (j!=16)
        {
            new_ro->pressure = (long)i*1000;
            new_ro->water_volume = (double) j/10000;

            if (calcul_v_i(*new_ro,e) > calcul_v_i(*ro,e))
            {
                ro->pressure = new_ro->pressure;
                ro->water_volume = new_ro->water_volume;
            }
            j+=1;
        }
        i+=1;
    }
    free(new_ro);
    
    temps = calcul_t_end_fusee(*ro,e);
    printf("\nres = %lf\n", temps);

    return ro;
}

fusee* fusee_max_t(env e, bouteille b){
    fusee* r = temps_max_inverse(e, b.empty_volume, b.radius, b.nozzle_radius, b.empty_mass);
    return r;
}