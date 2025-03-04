#include "_fusee.h"

int main(){
    //----------------------------------------------------------------------
    // Structures exemples pour les tests

    fusee r;
    r.drag_coeff = 0.5;
    r.empty_mass = 0.1;
    r.radius = 0.05;
    r.nozzle_radius = 0.011;
    r.pressure = 700000;
    r.water_volume = 0.0005;
    r.surface = 0;
    r.empty_volume = 0.002;

    bouteille b;
    b.empty_mass = 0.1;
    b.empty_volume = 0.002;
    b.nozzle_radius = 0.011;
    b.radius = 0.05;

    env e;
    e.g = 9.81;
    e.rho = 1.2;
    e.T_ext = 298.15;
    e.gamma = 1.4;
    e.P_ext = 100000;

    //----------------------------------------------------------------------
    // Tests

    double v_t = calcul_vt(r,e);
    printf("\nvt = %lf\n", v_t);
    printf("\nmax = %lf\n", calcul_z_top(3.28,20.4046,49.9,e));

    double T[3] = {0.0, 0.0, 0.0};
    expulsion_eau_v1(T, r, e);

    printf("\nv_plus_eau = %lf\n", T[1]);
    methode1(T,r,e);
    printf("\nv_plus_gas_m1 = %lf\n", T[1]);

    moy_methode(T,r,e);
 
    fusee* end = retrouve_fusee_z(e,60, b);
    printf("\n\n Volume d'eau à mettre (en m³) : %lf\n\n", end->water_volume);
    printf("\n\n Pression dans la bouteille (en Pa) : %ld\n\n", end->pressure);

    double test[3];
    methode1(test, *end, e);
    printf("\test = %lf\n", calcul_t_end_fusee(r,e));

    end = fusee_max_t(e,b);
    printf("\n\n TEMPS \n\n Volume d'eau à mettre (en m³) : %lf\n\n", end->water_volume);
    printf("\n\n Pression dans la bouteille (en Pa) : %ld\n\n", end->pressure);

    free(end);


    return EXIT_SUCCESS;
}