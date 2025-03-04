#include "_fusee.h"

double  calcul_vt(fusee r, env e){
    double res;
    res = (double)2*(r.empty_mass)*e.g;
    double a = (e.rho)*(r.drag_coeff)*r.radius*_PI*r.radius;
    res = res/a;
    res = sqrt(res);
    return res;
}

double calcul_t_top(double v_i, double v_t, env e, double t_i){ // Renvoie le temps mis par la fusee pour atteintre le sommet de sa trajectoire
    double at = atan(v_i/v_t);
    double mul = v_t/e.g;
    return at*mul+t_i;
}


double calcul_z_top(double z_i,double v_t, double v_i, env e){ // correct à environ 1 m près d'après test, renvoie l'altitude max atteinte par la fusee
    double v = v_i/v_t;
    double l = log(1+(v*v));
    double v_t2 = (v_t*v_t);
    double mul = v_t2/(2*e.g);
    return z_i+mul*l;
}

double calcul_v_i(fusee r, env e){ // Renvoie la vitesse de la fusee a la fin de sa phase de propulsion
    double T[3] = {0.0,0.0,0.0};
    methode1(T,r,e);
    return T[1];
}

double calcul_v_end(double vt, double z_max, env e){ // Renvoie la vitesse de la fusee a la fin de son mouvement
    double res = (double)2*e.g*z_max/pow(vt,2);
    res = 0-res;
    res = 1-exp(res);
    if(res <= 0){res = 0;}
    else{res = vt*sqrt(res);}
    return 0-res;
}

double calcul_t_end_fusee(fusee r, env e){ // Renvoie le temps mis par la fusee pour effectuer son mouvement en entier
    double T[3];
    methode1(T,r,e);
    double vt = calcul_vt(r,e);
    double t_top = calcul_t_top(T[1],vt,e,T[0]);
    double z_max = calcul_z_top_fusee(r,e);
    double v_end = calcul_v_end(vt,z_max,e);
    double res = (double) vt/e.g;
    res = res*atanh(0-(v_end/vt));
    return res+t_top;
}

double calcul_z_top_fusee(fusee r, env e){  // A le même effet que "calcul_z_top" mais prend en entrée un variable de type fusee
    if(r.pressure == 0 || r.water_volume == 0){return 0;}
    double vt = calcul_vt(r,e);
    double T[3];
    methode1 (T,r,e);
    double res = calcul_z_top(T[2], vt, T[1], e);
    return res;
}

//----------------------------------------------------------------------
// On essaye de prendre un angle de tire ce qui nous ramène à un problème en 2 dimensions
// On tente ici de considérer teta, étant un fonction du temps, comme une constante valant un peu moins que sa moitié (au vu de la trajectoire théorique de la fusee)
// Par rapport aux résultats d'expérimentations réelles, c'est un échec ...

double  calcul_vt_z(fusee r, env e, double teta){
    double res;
    res = (double)2*(r.empty_mass)*e.g;
    double a = (e.rho)*(r.drag_coeff)*r.radius*_PI*r.radius*sin(teta);
    res = res/a;
    res = sqrt(res);
    return res;
}

double  calcul_vt_x(fusee r, env e, double teta){
    double res;
    res = (double)2*(r.empty_mass)*e.g;
    double a = (e.rho)*(r.drag_coeff)*r.radius*_PI*r.radius*cos(teta);
    res = res/a;
    res = sqrt(res);
    return res;
}

// calcul_t_top reste la même fonction mais avec vt = vt_z ; de même pour v_end et t_end

double calcul_x_end(fusee r, env e, double teta){ // teta appartient à [0 ; pi/2]
    double T[3];
    methode1(T,r,e);
    double vt_z = calcul_vt_z(r,e,teta*0.4);
    double vt_x = calcul_vt_x(r,e,teta*0.4);

    double t_ap = calcul_t_top(T[1], vt_z,e,T[0]);
    double z_max = calcul_z_top(T[2], vt_z, T[1], e);
    double v_end = calcul_v_end(vt_z, z_max, e);
    double t_end = (double) vt_z/e.g;
    t_end = t_end*atanh(0-(v_end/vt_z));
    t_end += t_ap;

    double v_x_0 = T[1]*cos(teta);
    double terme1 = pow(v_x_0/vt_x,2); terme1 += 1;
    terme1 = log(terme1); terme1 = terme1*(pow(vt_x,2)/(2*e.g)); 
    
    double x_ap = v_x_0*T[0] + terme1;
    double terme2 = (t_ap-t_end)/vt_x;
    terme2 = cosh(terme2); terme2 = log(terme2);
    terme2 = terme2*(pow(vt_x,2)/e.g);


    return x_ap - terme2;

    
}

